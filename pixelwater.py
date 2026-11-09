"""🟦 pixelwater — ビット盤（bitboard）で落ちるドット絵の水を作り、PNG・GIF・絵文字で描く。

盤の 1 行は Python の int 1 つ。ビット x が列 x。壁（wall）・水（water）・向き（dir、1 = 右へ進む）を
行ごとのビットマスクで持ち、1 行まとめてビット演算で動かす。標準ライブラリだけ。

1 ティックの動き:

1. 落下: 下のセルが空いていれば 1 行落ちる。落ちられなければ斜め下の空きへ滑る
   （下の行から順に。1 ティックで 1 行まで。重力なので流体によらず毎ティック）
2. 横移動: 支えのある水が向き（dir）の方へ 1 マス進む。ふさがれていたら向きを反転する
   （同じ空きセルを左右から取り合ったら右から来た方を通す）。
   横移動の回数は流体ごとの整数スケジュール（Bresenham 型の積算）で決める

水のドット数は移動元を消して移動先に立てるだけなので、毎ティック popcount（int.bit_count）で完全に一致する。

⚠️ 流体の違いは「見た目のたとえ」: 横移動の頻度 = base_rate × ν_water / ν（動粘度が大きいほど横に広がる
回数が減る）。落下は流体によらず毎ティック 1 行。粘性流れの物理モデルではない。
"""

from __future__ import annotations

import hashlib
import json
import math
import random
import struct
import zlib
from pathlib import Path

from fixedq import ONE, to_q
from fluids import get_fluid

WALL, WATER, EMPTY = "#", "~", "."

# 色は Ironmate（myon-bioinformatics/Ironmate の pixel_grid.py）の色 ID と派生規則に揃える。
# 色番号: 0 空き（Ironmate では透明）/ 1 壁 / 2 水 / 3 落ちている水
IRONMATE_BASE = {  # pixel_grid.PALETTE_BASE の一部（main, shadow）
    "C01": ("#5BBD3E", "#418A2D"),  # grass
    "C02": ("#63B3FF", "#3F75CC"),  # water
    "C03": ("#FF6A00", "#CC3D00"),  # lava
    "C04": ("#A84DFF", "#7520B8"),  # poison
    "C09": ("#FFD23F", "#C99A17"),  # gold
    "C07": ("#9EA4AB", "#6D737B"),  # stone
}
# 色は自由。🟩 は Namekian Water、🟪 は Poison Water の仮置き（fluids の placeholder）
FLUID_COLOR_ID = {"water": "C02", "ethanol": "C09", "methanol": "C03", "namekian-water": "C01", "poison-water": "C04"}
EMOJI = {"water": "🟦", "ethanol": "🟨", "methanol": "🟧", "namekian-water": "🟩", "poison-water": "🟪"}


def _rgb(value: str) -> tuple[int, int, int]:
    return tuple(int(value[i:i + 2], 16) for i in (1, 3, 5))


def derive_palette(main: str, shadow: str) -> tuple:
    """Ironmate の pixel_grid.derive_palette と同じ: (透明, main, shadow, highlight, outline)。

    highlight = round((3·main + 2·255) / 5)、outline = round(shadow / 2)（四捨五入は整数で）。
    """
    m, sh = _rgb(main), _rgb(shadow)
    highlight = tuple((2 * (3 * c + 510) + 5) // 10 for c in m)
    outline = tuple((c + 1) // 2 for c in sh)
    return ((0, 0, 0, 0), m + (255,), sh + (255,), highlight + (255,), outline + (255,))


def palette(fluid: str) -> list[tuple[int, int, int, int]]:
    """色番号 0..3 の RGBA。空き = 白（GIF・PNG）、壁 = 石の outline、水 = main、落下中 = highlight。"""
    body = derive_palette(*IRONMATE_BASE[FLUID_COLOR_ID.get(fluid, "C02")])
    wall = derive_palette(*IRONMATE_BASE["C07"])[4]
    return [(255, 255, 255, 255), wall, body[1], body[3]]


# --- 場面（テキストの地図） --------------------------------------------------

def _box(width: int, height: int) -> list[list[str]]:
    g = [[EMPTY] * width for _ in range(height)]
    for y in range(height):
        g[y][0] = g[y][width - 1] = WALL
    g[height - 1] = [WALL] * width
    return g


def scene_map(name: str) -> list[str]:
    """場面の地図（'#' 壁 / '~' 水 / '.' 空き）。自作の地図は simulate(map_lines=...) に渡せる。"""
    if name == "dam-break":
        w, h = 32, 24
        g = _box(w, h)
        for y in range(4, h - 1):
            for x in range(1, 9):
                g[y][x] = WATER
    elif name == "moniac-drain":
        w, h = 32, 32
        g = _box(w, h)
        for y in range(2, 15):
            g[y][7] = g[y][24] = WALL
        for x in range(7, 25):
            if x not in (15, 16):  # 穴
                g[14][x] = WALL
        for y in range(3, 14):
            for x in range(8, 24):
                g[y][x] = WATER
    else:
        raise KeyError(f"unknown scene {name!r}; choose from {sorted(SCENES)}")
    return ["".join(row) for row in g]


SCENES = ("dam-break", "moniac-drain")


class Grid:
    """行ごとのビットマスク。y = 0 が上、ビット x が列 x。"""

    def __init__(self, width: int, height: int, wall: list[int], water: list[int], dir_: list[int]):
        self.width, self.height = width, height
        self.full = (1 << width) - 1
        self.wall, self.water, self.dir = wall, water, dir_

    def count(self) -> int:
        return sum(r.bit_count() for r in self.water)

    def overlap(self) -> int:
        return sum((w & s).bit_count() for w, s in zip(self.water, self.wall))

    def free(self, y: int) -> int:
        """y 行の空きセル。盤の下（y = height）は床としてふさがっている。"""
        if y >= self.height:
            return 0
        return self.full & ~(self.water[y] | self.wall[y])

    def falling(self) -> list[int]:
        return [self.water[y] & self.free(y + 1) for y in range(self.height)]

    def column_counts(self) -> list[int]:
        """列ごとの水のドット数（立っているビットだけをたどる）。"""
        cols = [0] * self.width
        for r in self.water:
            while r:
                low = r & -r
                cols[low.bit_length() - 1] += 1
                r ^= low
        return cols

    def snapshot(self) -> tuple[tuple[int, ...], tuple[int, ...]]:
        return tuple(self.water), tuple(self.falling())


def parse_map(lines, seed: int | None = 0) -> Grid:
    rows = [ln.rstrip("\n") for ln in lines if ln.strip("\n")]
    if not rows:
        raise ValueError("empty map")
    width = len(rows[0])
    if any(len(r) != width for r in rows):
        raise ValueError("map rows must have equal width")
    bad = set("".join(rows)) - {WALL, WATER, EMPTY}
    if bad:
        raise ValueError(f"unknown map characters: {sorted(bad)}")
    wall = [sum(1 << x for x, c in enumerate(r) if c == WALL) for r in rows]
    water = [sum(1 << x for x, c in enumerate(r) if c == WATER) for r in rows]
    if seed is None:  # 交互（偶数列は右、奇数列は左）
        pattern = sum(1 << x for x in range(0, width, 2))
        dir_ = [w & pattern for w in water]
    else:
        rng = random.Random(seed)
        dir_ = [w & rng.getrandbits(width) for w in water]
    return Grid(width, len(rows), wall, water, dir_)


# --- 1 ティック --------------------------------------------------------------

def fall(g: Grid) -> int:
    """下の行から順に、空きの上にある水を 1 行落とす。落ちられない水は斜め下の空きへ滑る（向きの側を先に）。

    向きも一緒に運ぶ。斜めは両側を見るので、すれ違う穴（下の行の空き）と上の水が必ず出会う。
    動いたドット数を返す。
    """
    moved = 0
    full = g.full
    for y in range(g.height - 2, -1, -1):
        below = g.free(y + 1)
        cand = g.water[y] & below
        if cand:
            moved += cand.bit_count()
            g.water[y] &= ~cand
            g.water[y + 1] |= cand
            carried = g.dir[y] & cand
            g.dir[y] &= ~cand
            g.dir[y + 1] |= carried
            below = g.free(y + 1)
        w = g.water[y]
        if not w:
            continue
        same = full & ~(w | g.wall[y])
        supported = w & ~below
        r_ok = (below >> 1) & (same >> 1)  # (x+1, y+1) と (x+1, y) が空き（角をすり抜けない）
        l_ok = (below << 1) & (same << 1) & full
        d = g.dir[y]
        go_r = supported & ((d & r_ok) | (~d & ~l_ok & r_ok))
        go_l = supported & ((~d & l_ok) | (d & ~r_ok & l_ok)) & full
        clash = (go_r << 1) & (go_l >> 1)
        go_l &= ~(clash << 1)
        if go_r | go_l:
            moved += go_r.bit_count() + go_l.bit_count()
            g.water[y] &= ~(go_r | go_l)
            g.dir[y] &= g.water[y]
            g.water[y + 1] |= ((go_r << 1) & full) | (go_l >> 1)
            g.dir[y + 1] |= (go_r << 1) & full
    return moved


def spread(g: Grid) -> int:
    """支えのある水を向きの方へ 1 マス。ふさがれたら向きを反転。動いたドット数を返す。"""
    moved = 0
    full = g.full
    for y in range(g.height):
        w = g.water[y]
        if not w:
            continue
        supported = w & ~g.free(y + 1)
        free = full & ~(w | g.wall[y])
        right = w & g.dir[y] & supported
        left = w & ~g.dir[y] & supported
        go_r = right & (free >> 1)
        go_l = left & ((free << 1) & full)
        clash = (go_r << 1) & (go_l >> 1)  # 同じ空きを左右から取り合う → 右から来た方を通す
        go_l &= ~(clash << 1)
        stay = w & ~go_r & ~go_l
        flip_r = right & ~go_r  # 右へ行けなかった → 左向きに
        flip_l = left & ~go_l  # 左へ行けなかった → 右向きに
        g.water[y] = stay | ((go_r << 1) & full) | (go_l >> 1)
        g.dir[y] = ((g.dir[y] & stay & ~flip_r) | flip_l | ((go_r << 1) & full)) & g.water[y]
        moved += go_r.bit_count() + go_l.bit_count()
    return moved


def lateral_rate(fluid: str, base_rate: float = 0.5, mu: float | None = None, mu_scale: float | None = None) -> int:
    """横移動の頻度（Q）。base_rate × ν_water / ν（🍯 mu・mu_scale で粘度を変えた流体にも使える）。"""
    nu_w = get_fluid("water").nu
    return to_q(base_rate * nu_w / get_fluid(fluid, mu=mu, mu_scale=mu_scale).nu)


# --- 回す --------------------------------------------------------------------

def run_pixels(scene: str = "dam-break", fluid: str = "water", *, ticks: int | None = None, seed: int | None = 0,
               base_rate: float = 0.5, frame_every: int = 2, map_lines=None, mu: float | None = None,
               mu_scale: float | None = None) -> dict:
    """シミュレーション本体。frames（各フレームの水・落下中の行マスク）を含む Python の dict を返す。"""
    if frame_every < 1:
        raise ValueError("frame_every must be >= 1")
    lines = map_lines if map_lines is not None else scene_map(scene)
    g = parse_map(lines, seed)
    rate = lateral_rate(fluid, base_rate, mu, mu_scale)
    if ticks is None:
        ticks = {"dam-break": 800, "moniac-drain": 600}.get(scene, 800)
    n0 = g.count()
    acc = 0
    counts_worst = 0
    lateral_passes = 0
    frames = [g.snapshot()]
    heights = []
    first_floor = None
    for t in range(1, ticks + 1):
        fall(g)
        acc += rate
        while acc >= ONE:
            acc -= ONE
            spread(g)
            lateral_passes += 1
        counts_worst = max(counts_worst, abs(g.count() - n0), g.overlap())
        cols = g.column_counts()
        inner = [c for c, s in zip(cols, _reachable_columns(g)) if s]
        heights.append(max(inner) - min(inner) if inner else 0)
        if first_floor is None and g.water[g.height - 2]:
            first_floor = t
        if t % frame_every == 0 or t == ticks:
            frames.append(g.snapshot())
    settle = _settle_tick(heights)
    digest = hashlib.sha256()
    for water, _f in frames:
        for row in water:
            digest.update(row.to_bytes((g.width + 7) // 8, "little"))
    return {
        "scene": scene, "fluid": fluid, "seed": seed, "ticks": ticks, "frame_every": frame_every,
        "width": g.width, "height": g.height,
        "pixels": n0, "count_error": counts_worst, "lateral_rate": rate / ONE, "lateral_passes": lateral_passes,
        "settle_tick": settle, "flatness": heights[-1] if heights else 0, "first_floor_tick": first_floor,
        "frames": frames, "wall": tuple(g.wall), "hash": digest.hexdigest(), "grid": g,
    }


def _reachable_columns(g: Grid) -> list[bool]:
    """底の行で壁でない列（水が溜まりうる列）。"""
    floor_row = g.height - 2
    return [not (g.wall[floor_row] >> x) & 1 for x in range(g.width)]


def _settle_tick(heights: list[int]) -> int | None:
    """列の高さの差（最大 − 最小）が 1 以下になり、その後ずっと 1 以下のままの最初のティック。"""
    if not heights or heights[-1] > 1:
        return None
    t = len(heights)
    while t > 0 and heights[t - 1] <= 1:
        t -= 1
    return t + 1


def simulate_pixels(scene: str = "dam-break", fluid: str = "water", ticks: int | None = None, seed: int | None = 0,
                    base_rate: float = 0.5) -> dict:
    """JSON にそのまま出せる要約（Node から ``--call simulate_pixels`` で呼べる）。最後の盤は絵文字の文字列で返す。"""
    res = run_pixels(scene, fluid, ticks=ticks, seed=seed, base_rate=base_rate)
    out = {k: v for k, v in res.items() if k not in ("frames", "wall", "grid")}
    out["frames"] = len(res["frames"])
    out["final"] = render_text(res, -1)
    return out


# --- 描画 -------------------------------------------------------------------

def indices(res: dict, frame: int) -> list[list[int]]:
    """フレームを色番号の 2 次元配列へ（0 空き / 1 壁 / 2 水 / 3 落ちている水）。"""
    water, falling = res["frames"][frame]
    w = res["width"]
    out = []
    for y, wall in enumerate(res["wall"]):
        row = []
        for x in range(w):
            bit = 1 << x
            row.append(1 if wall & bit else 3 if falling[y] & bit else 2 if water[y] & bit else 0)
        out.append(row)
    return out


def render_text(res: dict, frame: int = -1, ascii_only: bool = False, fill: str | None = None) -> str:
    """盤を絵文字で（⬜ 空き / ⬛ 壁 / 流体のブロック）。fill で好きなブロックに替えられる（色は自由）。"""
    if ascii_only:
        chars = [EMPTY, WALL, WATER, "|"]
    else:
        drop = fill or EMOJI.get(res["fluid"], "🟦")
        chars = ["⬜", "⬛", drop, drop]
    return "\n".join("".join(chars[c] for c in row) for row in indices(res, frame))


def to_ironmate_artifact(res: dict, frame: int = -1) -> dict:
    """フレームを Ironmate の ``ironmate-pixel-grid/1`` 形式へ（色番号 0 は透明）。

    Ironmate の grid_tools / grid_png（to_png・encode・emoji など）がそのまま読める。
    """
    colors = ["#%02X%02X%02X%02X" % c for c in palette(res["fluid"])]
    colors[0] = "#00000000"
    return {"schema": "ironmate-pixel-grid/1", "mask": res["scene"], "palette": FLUID_COLOR_ID.get(res["fluid"], "C02"),
            "effect": "normal", "width": res["width"], "height": res["height"], "colors": colors,
            "pixels": indices(res, frame)}


def _scale(grid: list[list[int]], scale: int) -> list[bytes]:
    """最近傍で拡大（1 マス = scale × scale ピクセル）。"""
    rows = []
    for row in grid:
        line = bytes(c for c in row for _ in range(scale))
        rows.extend([line] * scale)
    return rows


def png_bytes(res: dict, frame: int = -1, scale: int = 8) -> bytes:
    """パレット PNG（色番号 8 bit、フィルタ 0、zlib 9 で同じ入力なら同じバイト列）。"""
    if scale < 1:
        raise ValueError("scale must be >= 1")
    rows = _scale(indices(res, frame), scale)
    height, width = len(rows), len(rows[0])
    raw = b"".join(b"\x00" + r for r in rows)
    pal = [c[:3] for c in palette(res["fluid"])]

    def chunk(kind: bytes, data: bytes) -> bytes:
        return struct.pack(">I", len(data)) + kind + data + struct.pack(">I", zlib.crc32(kind + data) & 0xFFFFFFFF)

    return (b"\x89PNG\r\n\x1a\n"
            + chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 3, 0, 0, 0))
            + chunk(b"PLTE", b"".join(bytes(c) for c in pal))
            + chunk(b"IDAT", zlib.compress(raw, 9))
            + chunk(b"IEND", b""))


def _lzw(data: bytes, min_size: int) -> bytes:
    """GIF の LZW 符号化（可変長コード、LSB から詰める、表が満ちたら clear）。"""
    clear, eoi = 1 << min_size, (1 << min_size) + 1
    size = min_size + 1
    table: dict[tuple[int, int], int] = {}
    nxt = eoi + 1
    acc = nbits = 0
    out = bytearray()

    def emit(code: int) -> None:
        nonlocal acc, nbits
        acc |= code << nbits
        nbits += size
        while nbits >= 8:
            out.append(acc & 0xFF)
            acc >>= 8
            nbits -= 8

    emit(clear)
    prefix = data[0]
    for k in data[1:]:
        key = (prefix, k)
        code = table.get(key)
        if code is not None:
            prefix = code
            continue
        emit(prefix)
        table[key] = nxt
        nxt += 1
        if nxt > (1 << size) and size < 12:
            size += 1
        if nxt == 4096:
            emit(clear)
            table.clear()
            nxt, size = eoi + 1, min_size + 1
        prefix = k
    emit(prefix)
    emit(eoi)
    if nbits:
        out.append(acc & 0xFF)
    return bytes(out)


def _blocks(data: bytes) -> bytes:
    return b"".join(bytes([len(data[i:i + 255])]) + data[i:i + 255] for i in range(0, len(data), 255)) + b"\x00"


def gif_bytes(res: dict, scale: int = 8, delay_cs: int = 4, frames=None) -> bytes:
    """アニメーション GIF。2 枚目からは前のフレームと違う矩形だけを重ねる（容量を抑える）。"""
    if scale < 1:
        raise ValueError("scale must be >= 1")
    picks = range(len(res["frames"])) if frames is None else frames
    grids = [indices(res, i) for i in picks]
    h, w = len(grids[0]) * scale, len(grids[0][0]) * scale
    pal = [c[:3] for c in palette(res["fluid"])]
    out = bytearray(b"GIF89a" + struct.pack("<HHBBB", w, h, 0xF1, 0, 0))  # 全体の色表 4 色
    out += b"".join(bytes(c) for c in pal)
    out += b"\x21\xff\x0bNETSCAPE2.0\x03\x01\x00\x00\x00"  # 無限ループ
    prev = None
    for grid in grids:
        if prev is None:
            x0, y0, x1, y1 = 0, 0, len(grid[0]), len(grid)
        else:
            diff = [(x, y) for y, (a, b) in enumerate(zip(grid, prev)) for x in range(len(a)) if a[x] != b[x]]
            if not diff:
                diff = [(0, 0)]
            xs, ys = [d[0] for d in diff], [d[1] for d in diff]
            x0, y0, x1, y1 = min(xs), min(ys), max(xs) + 1, max(ys) + 1
        sub = [row[x0:x1] for row in grid[y0:y1]]
        pix = b"".join(_scale(sub, scale))
        out += b"\x21\xf9\x04" + struct.pack("<BHBB", 0x04, delay_cs, 0, 0)  # 処分しない（重ねる）
        out += b"\x2c" + struct.pack("<HHHHB", x0 * scale, y0 * scale, (x1 - x0) * scale, (y1 - y0) * scale, 0)
        out += b"\x02" + _blocks(_lzw(pix, 2))
        prev = grid
    out += b"\x3b"
    return bytes(out)


# --- 読み戻し（検証用の最小デコーダ） ----------------------------------------

def read_png(data: bytes) -> dict:
    """署名・各チャンクの CRC・IHDR・展開後の長さを確かめ、色番号の行を返す。"""
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError("bad PNG signature")
    pos, chunks = 8, []
    while pos < len(data):
        n = struct.unpack(">I", data[pos:pos + 4])[0]
        kind, body = data[pos + 4:pos + 8], data[pos + 8:pos + 8 + n]
        crc = struct.unpack(">I", data[pos + 8 + n:pos + 12 + n])[0]
        if zlib.crc32(kind + body) & 0xFFFFFFFF != crc:
            raise ValueError(f"bad CRC in {kind!r}")
        chunks.append((kind, body))
        pos += 12 + n
    kinds = [k for k, _ in chunks]
    if kinds[0] != b"IHDR" or kinds[-1] != b"IEND":
        raise ValueError("PNG must start with IHDR and end with IEND")
    width, height = struct.unpack(">II", chunks[0][1][:8])
    raw = zlib.decompress(b"".join(b for k, b in chunks if k == b"IDAT"))
    if len(raw) != height * (width + 1):
        raise ValueError("PNG data length mismatch")
    rows = [raw[y * (width + 1) + 1:(y + 1) * (width + 1)] for y in range(height)]
    return {"width": width, "height": height, "chunks": [k.decode() for k in kinds], "rows": rows}


def _unlzw(data: bytes, min_size: int, count: int) -> bytes:
    clear, eoi = 1 << min_size, (1 << min_size) + 1
    size, pos, out = min_size + 1, 0, bytearray()
    table: list[bytes] = []
    prev = None
    bits = int.from_bytes(data, "little")
    total = len(data) * 8

    def reset():
        nonlocal table, size, prev
        table = [bytes([i]) for i in range(clear)] + [b"", b""]
        size, prev = min_size + 1, None

    reset()
    while pos + size <= total and len(out) < count:
        code = (bits >> pos) & ((1 << size) - 1)
        pos += size
        if code == clear:
            reset()
            continue
        if code == eoi:
            break
        if code < len(table):
            entry = table[code]
            if prev is not None:
                table.append(prev + entry[:1])
        elif prev is not None and code == len(table):
            entry = prev + prev[:1]
            table.append(entry)
        else:
            raise ValueError("bad LZW code")
        out += entry
        prev = entry
        if len(table) == (1 << size) and size < 12:
            size += 1
    return bytes(out)


def read_gif(data: bytes) -> dict:
    """GIF を読み、各フレームを重ねた後の色番号（全画面）を返す。"""
    if data[:6] != b"GIF89a":
        raise ValueError("bad GIF signature")
    w, h, flags = struct.unpack("<HHB", data[6:11])
    pos = 13 + (3 * (2 << (flags & 7)) if flags & 0x80 else 0)
    canvas = bytearray(w * h)
    frames = []
    while pos < len(data):
        b = data[pos]
        if b == 0x3B:
            break
        if b == 0x21:  # 拡張
            pos += 2
            while data[pos]:
                pos += data[pos] + 1
            pos += 1
            continue
        if b != 0x2C:
            raise ValueError(f"unexpected block 0x{b:02x}")
        x, y, fw, fh, _ = struct.unpack("<HHHHB", data[pos + 1:pos + 10])
        min_size = data[pos + 10]
        pos += 11
        body = bytearray()
        while data[pos]:
            body += data[pos + 1:pos + 1 + data[pos]]
            pos += data[pos] + 1
        pos += 1
        pix = _unlzw(bytes(body), min_size, fw * fh)
        if len(pix) != fw * fh:
            raise ValueError("GIF frame size mismatch")
        for r in range(fh):
            canvas[(y + r) * w + x:(y + r) * w + x + fw] = pix[r * fw:(r + 1) * fw]
        frames.append(bytes(canvas))
    return {"width": w, "height": h, "frames": frames}


# --- ファイルへ ---------------------------------------------------------------

RECEIPT_SCHEMA = "fluidfullspeedsimulation-pixel-water/1"


def file_stem(scene: str, fluid: str, seed, ticks: int, out_dir: str | Path = ".", file_prefix: str = "",
              file_suffix: str = "") -> Path:
    """意味のある名前で採番する: <prefix><場面>_<流体>_s<seed>_t<ティック>_r<NN><suffix>（拡張子なし）。

    並びは Ironmate の ``{mask}_{colorID}_{effect}_{size}`` と同じく条件を ``_`` でつなぐ。
    r<NN> は同じ条件で何回目の出力か（GIF・PNG・受領記録 JSON のどれかが既にあれば次の番号）。
    """
    base = f"{file_prefix}{scene}_{fluid}_s{'alt' if seed is None else seed}_t{ticks:04d}"
    out = Path(out_dir)
    n = 1
    while any((out / f"{base}_r{n:02d}{file_suffix}.{ext}").exists() for ext in ("gif", "png", "json")):
        n += 1
    return out / f"{base}_r{n:02d}{file_suffix}"


def draw(scene: str = "dam-break", fluid: str = "water", out_dir: str = "build/pixel", *, ticks: int | None = None,
         seed: int | None = 0, scale: int = 8, frame_every: int = 2, formats=("gif", "png"), file_prefix: str = "",
         file_suffix: str = "", base_rate: float = 0.5) -> dict:
    """GIF（動き）と PNG（最後のフレーム）を書き、隣に受領記録 <name>.json を置く。

    受領記録は BrowserTestKit の evidence_gif と同じ考え方: 何をどう描いたか（plan）と、出来たもの（bytes・sha256・
    frames・canvas）を残す。plan の値で draw を呼べば同じバイト列になる。
    """
    res = run_pixels(scene, fluid, ticks=ticks, seed=seed, frame_every=frame_every, base_rate=base_rate)
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    stem = file_stem(scene, fluid, seed, res["ticks"], out_dir, file_prefix, file_suffix)
    paths, files = {}, {}
    blobs = {}
    if "gif" in formats:
        blobs["gif"] = gif_bytes(res, scale=scale)
    if "png" in formats:
        blobs["png"] = png_bytes(res, -1, scale=scale)
    for kind, data in blobs.items():
        p = stem.with_name(f"{stem.name}.{kind}")
        p.write_bytes(data)
        paths[kind] = str(p)
        files[kind] = {"name": p.name, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
    receipt = {
        "schema": RECEIPT_SCHEMA, "files": files, "frames": len(res["frames"]),
        "canvas": f"{res['width'] * scale}x{res['height'] * scale}",
        "plan": {"scene": scene, "fluid": fluid, "ticks": res["ticks"], "seed": seed, "scale": scale,
                 "frame_every": frame_every, "base_rate": base_rate},
        "pixels": res["pixels"], "count_error": res["count_error"], "settle_tick": res["settle_tick"],
        "hash": res["hash"],
    }
    p = stem.with_name(f"{stem.name}.json")
    p.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    paths["json"] = str(p)
    return {"paths": paths, **{k: receipt[k] for k in ("frames", "canvas", "pixels", "count_error", "settle_tick",
                                                       "hash")}}


# --- 🟦 を積む断面図（行 = 高さ、列 = 区間・水槽）で、溜まる・流れる・遅れるを見る ---------

SECTION_ROWS = 8  # 縦の段数（8 の倍数を既定に。引数で変えられる）
EMPTY_BLOCK, WALL_BLOCK = "⬜", "⬛"
CARRIED, IN_PIPE = "🍜", "⏳"  # 運ばれる物 / 管の中
SECTION_CASES = ("flume", "flume-hopper")
CAP_MODES = ("run-max", "rim")


def fill_block(fluid: str) -> str:
    """流体ごとの既定のブロック（色は自由。--fill で好きな絵文字に替えてよい）。"""
    return EMOJI.get(fluid, "🟦")


def blocks(h: float, cap: float, rows: int = SECTION_ROWS) -> int:
    """水深 h を段数へ: round(h / cap × rows)。0..rows に収める（cap を超えたら rows）。"""
    if cap <= 0 or h <= 0:
        return 0
    return max(0, min(rows, int(h / cap * rows + 0.5)))


def readback(n: int, cap: float, rows: int = SECTION_ROWS) -> float:
    """段数から水深を読み戻す（1 段 = cap / rows。ずれは ±½ 段まで）。"""
    return n * cap / rows


def _columns(res: dict, columns) -> list[tuple[str, str, object]]:
    """列の指定を展開: 樋の prefix → 各セル、水槽名 → その水位、'src>dst' → 管（水深なし）。"""
    row0 = res["series"][0] if res["series"] else {}
    out = []
    for spec in columns:
        if ">" in spec:
            out.append(("pipe", spec, tuple(spec.split(">", 1))))
        elif f"h:{spec}" in row0 and isinstance(row0[f"h:{spec}"], list):
            out += [("cell", f"{spec}:{i}", (spec, i)) for i in range(len(row0[f"h:{spec}"]))]
        elif f"h:{spec}" in row0:
            out.append(("tank", spec, spec))
        else:
            raise ValueError(f"unknown section column {spec!r}")
    return out


def tracer_intervals(res: dict) -> list[list[tuple[float, float, str]]]:
    """各トレーサーの居場所の区間 [(始まり, 終わり, 場所)]。遅れのある管の中は 'src>dst'。"""
    delays = {(ln["src"], ln["dst"]): ln["delay"] for ln in res["links"]}
    out = []
    for tr in res["tracers"]:
        spans = []
        path = tr["path"]
        for a, b in zip(path, path[1:]):
            d = delays.get((a["node"], b["node"]), 0.0)
            spans.append((a["t"], b["t"] - d, a["node"]))
            if d:
                spans.append((b["t"] - d, b["t"], f"{a['node']}>{b['node']}"))
        if path:
            # 取られた（箸）・はじき出された麺はそこで消える
            end = tr.get("until")
            spans.append((path[-1]["t"], math.inf if end is None else end, path[-1]["node"]))
        out.append(spans)
    return out


def _caps(res: dict, cols, depths, cap_mode: str, caps: dict | None) -> list[float]:
    """列ごとの「縦いっぱい」の水深。caps（名前 → 水深）が最優先。run-max は走らせた中の最大
    （樋はセルどうしを比べられるよう樋全体で 1 つ）、rim は樋の縁（水槽は run-max と同じ）。"""
    if cap_mode not in CAP_MODES:
        raise ValueError(f"cap mode must be one of {CAP_MODES}")
    caps = caps or {}
    group_max: dict[str, float] = {}
    for c, (kind, name, ref) in enumerate(cols):
        key = ref[0] if kind == "cell" else name
        group_max[key] = max(group_max.get(key, 0.0), max((row[c] or 0.0) for row in depths) if depths else 0.0)
    out = []
    for kind, name, ref in cols:
        key = ref[0] if kind == "cell" else name
        if key in caps:
            out.append(float(caps[key]))
        elif kind == "cell" and cap_mode == "rim":
            out.append(res["channels"][key]["depth_cap"])
        else:
            out.append(group_max.get(key, 0.0))
    return out


def _pick(n: int, count: int | None) -> list[int]:
    if count is None or count >= n:
        return list(range(n))
    if count < 2:
        return [n - 1]
    return sorted({round(i * (n - 1) / (count - 1)) for i in range(count)})


def network_frames(res: dict, columns, *, frames: int | None = 6, rows: int = SECTION_ROWS,
                   cap_mode: str = "run-max", caps: dict | None = None, fill: str | None = None,
                   empty: str = EMPTY_BLOCK, carried: str = CARRIED, in_pipe: str = IN_PIPE) -> dict:
    """network.Network.run の結果を、時刻ごとの断面図（🟦 を下から積む）にする。

    series の各行を時刻とみなす（細かく見るなら run の sample_every を小さく）。frames は並べる時刻の数。
    arrivals は列ごとに最初の 1 段が積まれた時刻（水の先端が届いた時刻）、readback_error は全時刻での
    |段数 × 1 段 − 水深| の最大（縦いっぱいを超えていない分だけ）。
    """
    if rows < 1:
        raise ValueError("rows must be >= 1")
    fill = fill or fill_block(res.get("fluid", "water"))
    cols = _columns(res, columns)
    names = [name for _k, name, _x in cols]
    series = res["series"]
    times = [row["t"] for row in series]
    depths = [[None if kind == "pipe" else row[f"h:{ref[0]}"][ref[1]] if kind == "cell" else row[f"h:{ref}"]
               for kind, _n, ref in cols] for row in series]
    cap = _caps(res, cols, depths, cap_mode, caps)
    counts = [[0 if d is None else blocks(d, cap[c], rows) for c, d in enumerate(line)] for line in depths]
    worst = 0.0
    for line, cnt in zip(depths, counts):
        for c, d in enumerate(line):
            if d is not None and cap[c] > 0 and d <= cap[c]:
                worst = max(worst, abs(readback(cnt[c], cap[c], rows) - d))
    arrivals = [next((t for t, cnt in zip(times, counts) if cnt[c] > 0), None) if cols[c][0] != "pipe" else None
                for c in range(len(cols))]
    spans = tracer_intervals(res)
    where = []
    for t in times:
        here = []
        for tr_spans in spans:
            loc = next((place for a, b, place in tr_spans if a <= t < b), None)
            if loc in names:
                here.append((names.index(loc), ">" in loc))
            elif loc and ">" in loc and loc.split(">", 1)[1] in names:  # 管の列が無ければ行き先の列の底に ⏳
                here.append((names.index(loc.split(">", 1)[1]), True))
            else:
                here.append(None)
        where.append(here)
    shown = []
    for r in _pick(len(times), frames):
        grid = [[empty] * len(cols) for _ in range(rows)]
        for c, n in enumerate(counts[r]):
            for k in range(n):
                grid[rows - 1 - k][c] = fill
        spots = [s_ for s_ in where[r] if s_ is not None]
        for body in series[r].get("bodies", {}).values():  # 長さのある麺は体の区間すべてに 🍜
            spots += [(names.index(cell), False) for cell in body if cell in names]
        for c, piped in spots:
            k = max(counts[r][c], 1) - 1  # 水面の段（空なら底）
            grid[rows - 1 - k][c] = in_pipe if piped else carried
        shown.append({"t": times[r], "lines": ["".join(row) for row in grid], "counts": counts[r]})
    return {"columns": names, "caps": cap, "rows": rows, "fill": fill, "times": times, "counts": counts,
            "frames": shown, "arrivals": arrivals, "readback_error": worst,
            "half_block": max((c / (2 * rows) for c in cap if c > 0), default=0.0),
            "tracer_columns": [[None if s is None else s[0] for s in here] for here in where],
            "tracer_spans": [[(a, None if b == math.inf else b, p) for a, b, p in sp] for sp in spans]}


def pixel_arrivals(res: dict) -> list[int | None]:
    """ドット絵の盤で、列ごとに最初の水のドットが来たティック（壁で届かない列は除く）。"""
    keep = [x for x, ok in enumerate(_reachable_columns(res["grid"])) if ok]
    every = res["frame_every"]
    last = len(res["frames"]) - 1
    out = {}
    for i, (water, _f) in enumerate(res["frames"]):
        tick = res["ticks"] if i == last else i * every
        seen = 0
        for r in water:
            seen |= r
        for x in keep:
            if x not in out and seen >> x & 1:
                out[x] = tick
    return [out.get(x) for x in keep]


def pixel_frames(res: dict, *, frames: int | None = 6, fill: str | None = None) -> dict:
    """ドット絵の水を同じ絵文字で時刻ごとに並べる（盤そのものが断面図）。"""
    picks = _pick(len(res["frames"]), frames)
    last = len(res["frames"]) - 1
    out = []
    for i in picks:
        text = render_text(res, i, fill=fill)
        out.append({"t": float(res["ticks"] if i == last else i * res["frame_every"]), "lines": text.splitlines()})
    return {"frames": out, "arrivals": pixel_arrivals(res)}


def section_case(name: str = "flume", fluid: str = "water", *, dt: float = 0.05, mu: float | None = None,
                 mu_scale: float | None = None) -> tuple[dict, list[str]]:
    """断面図を試す回路（network の部品だけで組む）。戻り値は (run の結果, 列)。

    - flume: 空の樋（長方形、層流、20 区間）へポンプで流し込み、10 秒後にトレーサー 🍜 を放す
    - flume-hopper: ポンプ → じょうご（小さな水槽）→ 遅れ 3 秒の細管 → 樋。じょうごに溜まってから流れる
    """
    from network import Network  # 遅延 import（pixelwater 単体はセル・オートマトンだけで使える）

    net = Network(dt, fluid, mu=mu, mu_scale=mu_scale).tank("sea", 1.0, volume=1.0)
    net.channel("flume", "sea", length=2.0, cells=20, slope=0.05, section={"shape": "rect", "size": 0.05},
                law="laminar", max_depth=0.005)
    if name == "flume":
        net.pump("sea", "flume:0", rate=2e-6)
        net.tracer("somen", "flume:0", release=10.0, stop="sea")
        return net.run(60.0, sample_every=1), ["flume"]
    if name == "flume-hopper":
        net.tank("hopper", 1e-4)
        net.pump("sea", "hopper", rate=2e-6)
        net.pipe("hopper", "flume:0", radius=1e-3, length=0.5, delay=3.0)
        net.tracer("somen", "hopper", release=20.0, stop="sea")
        return net.run(90.0, sample_every=1), ["hopper", "hopper>flume:0", "flume"]
    raise KeyError(f"unknown section case {name!r}; choose from {SECTION_CASES}")


def render_frames(view: dict) -> str:
    """断面図を時刻の順に縦へ並べる（各図の頭に ⏱️ 時刻）。"""
    out = []
    for fr in view["frames"]:
        out.append(f"⏱️ {fr['t']:.2f}")
        out += fr["lines"]
    return "\n".join(out)
