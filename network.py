"""🌊 network — 水槽・樋（ストック）と管（フロー）をつないだ汎用の水力回路。

MONIAC は「水が流れ落ちる様子」の手本。部品は物理のまま持ち、組み合わせで流しそうめんの台のような
別の題材にも広げられるようにする。

節点（ノード）:

- ``tank``: 体積を Q 整数で持つ水槽。水位 h = V / A。流入と流出の差を積み上げる積分器そのもの
- ``channel``: 傾いた樋を ``cells`` 個の区間セル（``<prefix>:0`` …）に分けたもの。断面は長方形・半円（竹）
- ``flume``: 竹を継いだ流しそうめんの台。継ぎ目は段差（自由落下）か縁（たまってから越える）、曲がりはあふれうる

管（リンク）の法則は ``LAWS`` 登録簿にあり、``register_law`` で増やせる:

- ``pipe``（Hagen–Poiseuille、層流）/ ``orifice``（Torricelli）/ ``weir_rect`` ``weir_v``（堰）
- ``cam``（水位 → 流量の表）/ ``rating``（セル体積 → 流量の表。樋の区間をつなぐ）
- ``lip``（継ぎ目の縁を越える堰）/ ``spill``（曲がりの外側から横へあふれる）
- ``pump``（一定流量）/ ``float_valve``（浮き弁。測った水位が遅れの段を通って弁を開く）
- ``valve``（弁つきの出口。開き具合は一定か、別の水槽の浮き → カムで決まる。MONIAC の弁）
- ``tip``（🎍 鹿威し: 溜まると傾いて一気に落とす）/ ``siphon``（🍷 ピタゴラスの杯: 頂に届くと吸い切る）
  / ``mariotte``（🫗 マリオットびん: 空気管の下端まで流量が一定）
- ``tip``（🎍 鹿威し: 溜まると傾いて一気に落とす）/ ``siphon``（🍷 ピタゴラスの杯: 頂に届くと吸い切る）
  / ``mariotte``（🫗 マリオットびん: 空気管の下端まで流量が一定）
- どの管にも ``delay``（輸送の遅れ）を付けられる。管の中の水も総量に数える

運ばれる物 ``tracer``（そうめん・葉っぱ・浮き）: 入った節点に先にいた水が全部出たら、いちばん多く
流れた管で次へ進む（押し出し流れ・先入れ先出し）。水槽では平均滞留時間 V/Q、樋では L/v に当たる。
太さ ``size`` を持たせると、継ぎ目の縁の越流水深がそれより浅いあいだ引っかかる（``PASSABLE``）。

1 ステップでは全部の管の流量を同じ状態から求め（同時更新）、出す側の水が足りなければ整数のまま按分する。
移動は引いた整数をそのまま足すだけなので、総量は 1 ulp も変わらない。``build(spec)`` で JSON から組める。
"""

from __future__ import annotations

import bisect
import math
import random
from typing import Callable

from fractions import Fraction

from fixedq import ONE, Cam, DelayLine, clamp, exp_neg_q, from_q, lag_gain, lag_step, qdiv, qmul, qsqrt, to_q
from fluids import G, RE_LAMINAR, Section, channel_flow, fall_time, get_fluid, lip_flow, plan_width, reynolds
from fluids import side_spill
from fluids import viscous_orifice
from fluids import superelevation

LAWS: dict[str, Callable] = {}
OUTLET_MODES = ("exact", "independent")
INFLOW_MODES = ("end", "uniform")
ROUNDING_MODES = ("q48", "carry64")
CARRY_BITS = 64  # 🧮 carry64 の係数の刻み（2^-64）。量は Q（2^-48）のまま
PASSABLE: dict[str, Callable] = {}  # 運ばれる物がその管を通れるか（net, link, tracer）→ bool。無ければ常に通れる
EJECT: dict[str, Callable] = {}  # 運ばれる物がその管へはじき出されるか（net, link, tracer）→ bool

# 麺のたまり（縁に積もった麺）の詰め方: 平行な円柱の最密充填 π/(2√3)（六方。いちばん低く積もる側の仮定）
PACKING = math.pi / (2.0 * math.sqrt(3.0))


def register_law(kind: str):
    """管の法則を登録する。関数は (net, link, n) を受け、このステップに送りたい量（Q）を返す。"""

    def deco(fn):
        LAWS[kind] = fn
        return fn

    return deco


@register_law("pipe")
def _law_pipe(net, ln, n):
    h = net.head_q(ln["src"])
    ln["max_head"] = max(ln["max_head"], h)
    nd = net.nodes[ln["src"]]
    if nd["profile"] is None:
        return qmul(nd["v"], ln["k"])  # 断面一定: 厳密離散化
    return qmul(qmul(ln["c_q"], h), net.dt_q)  # 断面が変わる器: Q = c·h を前進 Euler で


ORIFICE_RE_STEADY = 1.0e4  # これより低い Re では鋭いふちの小孔の Cd が Re で大きく変わる（Sci. Rep. 2023 ほか）


@register_law("orifice")
def _law_orifice(net, ln, n):
    h = net.head_q(ln["src"])
    if h <= 0:
        return 0
    ln["max_head"] = max(ln["max_head"], h)
    if not ln["viscous"]:
        return qmul(ln["coef"], qsqrt(h))
    hf = from_q(h)
    cd, re = viscous_orifice(hf, ln["diameter"], net.fluid.nu, ln["cd"], ln["delta"], net.g)
    ln["min_cd"] = cd if ln["min_cd"] is None else min(ln["min_cd"], cd)
    return qmul(to_q(cd * ln["orifice_area"] * math.sqrt(2.0 * net.g * hf)), net.dt_q)


@register_law("weir_rect")
def _law_weir_rect(net, ln, n):
    h = net.head_q(ln["src"])
    return qmul(ln["coef"], qmul(h, qsqrt(h))) if h > 0 else 0


@register_law("weir_v")
def _law_weir_v(net, ln, n):
    h = net.head_q(ln["src"])
    return qmul(ln["coef"], qmul(qmul(h, h), qsqrt(h))) if h > 0 else 0


@register_law("cam")
def _law_cam(net, ln, n):
    h = net.head_q(ln["src"])
    return qmul(ln["cam"](h), net.dt_q) if h > 0 else 0


@register_law("rating")
def _law_rating(net, ln, n):
    v = net.nodes[ln["src"]]["v"]
    ln["max_v"] = max(ln["max_v"], v)
    return qmul(ln["cam"](v), net.dt_q)


@register_law("pump")
def _law_pump(net, ln, n):
    if n < ln["start_step"] or (ln["until_step"] is not None and n >= ln["until_step"]):
        return 0
    return ln["amount"]


@register_law("lip")
def _law_lip(net, ln, n):
    """継ぎ目の段（次の竹の縁が lip だけ高い）。セルの水深が縁を越えた分だけ流れる。
    縁に麺がたまっていれば、そのぶん縁が高くなる（crest = lip + たまりの高さ）。"""
    v = net.nodes[ln["src"]]["v"]
    ln["max_v"] = max(ln["max_v"], v)
    h = from_q(net.head_q(ln["src"]))
    return qmul(to_q(lip_flow(ln["section"], h, crest(net, ln), ln["cd"], net.g, width=ln["width"])), net.dt_q)


def pile_members(net, ln) -> list[dict]:
    """その縁に引っかかって積もっている麺（長さのある麺だけ。長さ 0 の点は重さも体積も持たない）。"""
    return [tr for tr in net.tracers if tr["held"] and tr["held_link"] is ln and tr["length"] > 0]


def pile_height(net, ln) -> float:
    """縁に積もった麺の高さ Δ = Σ(本数·πd²/4) / (PACKING·w)。w は縁の高さでの水面幅（たまりの間は一定とみなす）。"""
    area = sum(tr["strands"] * math.pi * tr["size"] ** 2 / 4.0 for tr in pile_members(net, ln))
    return area / (PACKING * ln["width"]) if area else 0.0


def crest(net, ln) -> float:
    return ln["lip"] + pile_height(net, ln)


def _lip_passable(net, ln, tr):
    """越流水深（水深 − 縁。麺が積もっていればその上）が運ばれる物の太さ以上なら越えられる。
    押し流された（flushed）麺はその縁をそのまま越える。"""
    if id(ln) in tr["flushed"]:
        return True
    return from_q(net.head_q(ln["src"])) - crest(net, ln) >= tr["size"]


PASSABLE["lip"] = _lip_passable


@register_law("spill")
def _law_spill(net, ln, n):
    """樋の縁からあふれる。曲がり（radius あり）では外側の水面 h + Δh/2（Δh = v²B/(gR)）が縁を越えた分が
    横越流堰として流れ、どのセルでも縁の高さの容量を超えた体積はそのまま外へ落ちる。"""
    nd = net.nodes[ln["src"]]
    sec = nd["section"]
    h = from_q(net.head_q(ln["src"]))
    a = sec.area(h)
    amount = 0
    ln["last_outer"] = h
    if a > 0 and ln["radius"] is not None:
        v = from_q(qdiv(net.last_flow[ln["main"]], net.dt_q)) / a  # そのセルの本流（前のステップ）の流速
        h_out = h + superelevation(v, sec.top_width(h), ln["radius"], net.g) / 2.0
        ln["max_outer"] = max(ln["max_outer"], h_out)
        ln["last_outer"] = h_out
        amount = qmul(to_q(side_spill(nd["length"], h_out, ln["rim"], ln["cd"], net.g)), net.dt_q)
    return max(amount, nd["v"] - ln["v_cap"])


def _spill_eject(net, ln, tr):
    """曲がりの外側の水面に乗った麺は、その上面（水面 + 太さ）が縁に届くと縁の外へはじき出される。
    （麺は水面に浮いて外側の壁へ寄せられる、とみなす。水そのものがあふれるより太さの分だけ早い）"""
    return ln["radius"] is not None and ln.get("last_outer", 0.0) + tr["size"] >= ln["rim"]


EJECT["spill"] = _spill_eject


@register_law("valve")
def _law_valve(net, ln, n):
    """弁つきの出口: 全開で時定数 tau の線形流出に、開き具合（0..1）を掛ける。
    開き具合は一定か、別の水槽の浮き（sensor の水位）→ カム（表）で決まる。浮きと弁の間に遅れの段も持てる。"""
    if ln["cam"] is None:
        o = ln["opening"]
    else:
        u = net.head_q(ln["sensor"])
        for i in range(ln["stages"]):
            ln["z"][i] = lag_step(ln["z"][i], u, ln["k_link"])
            u = ln["z"][i]
        o = clamp(ln["cam"](u), 0, ONE)
    ln["last_opening"] = o
    return qmul(qmul(net.nodes[ln["src"]]["v"], ln["k"]), o)


@register_law("boil")
def _law_boil(net, ln, n):
    """🔥 単式蒸留の釜: 熱 power [W] で 1 ステップに P·dt / λ(y) mol の蒸気をつくり、気液平衡の y で
    エタノールと水に分けて dst（受け器）へ送る。成分モード（ethanol・water、モル）専用。"""
    from fluids import latent_heat, vle_y
    parts = net.nodes[ln["src"]]["parts"]
    ne, nw = parts["ethanol"], parts["water"]
    if ne + nw <= 0 or n < ln["start_step"] or (ln["until_step"] is not None and n >= ln["until_step"]):
        ln["_want"] = {}
        return 0
    x = ne / (ne + nw)
    y = vle_y(x)
    vapor = to_q(ln["power"] * net.dt / latent_heat(y))
    e = qmul(vapor, to_q(y))
    ln["last_x"], ln["last_y"] = x, y
    if ln["cuts"]:
        _cut(net, ln, n, y)
    ln["_want"] = {"ethanol": e, "water": vapor - e}
    return vapor


def _cut(net, ln, n, y):
    """✂️ 蒸留の切り替え: 流れ出る留液の度数（running ABV）か、その受け器にたまった体積で、次の受け器へ移る。"""
    from fluids import abv, mix_volume
    running = abv(y, 1.0 - y)
    ln["last_abv"] = running
    while ln["cut_index"] < len(ln["cuts"]) - 1:
        cut = ln["cuts"][ln["cut_index"]]
        rec = net.nodes[cut["to"]]["parts"]
        done = False
        if "abv_below" in cut and running < cut["abv_below"]:
            done = True
        if "volume" in cut and mix_volume(from_q(rec["ethanol"]), from_q(rec["water"])) >= cut["volume"]:
            done = True
        if not done:
            break
        cut.setdefault("ended", n * net.dt)
        ln["cut_index"] += 1
        ln["cuts"][ln["cut_index"]].setdefault("started", n * net.dt)
    ln["dst"] = ln["cuts"][ln["cut_index"]]["to"]


@register_law("evaporate")
def _law_evaporate(net, ln, n):
    """👼 成分ごとの 1 次の抜け: 成分 s は 1 ステップに parts_s·(1 − e^(−k_s·dt)) だけ出ていく（厳密離散化）。
    rates は {成分: k [1/s]}。成分モード専用。"""
    parts = net.nodes[ln["src"]]["parts"]
    out = {sp: qmul(parts[sp], g) for sp, g in ln["gains"].items()}
    ln["_want"] = out
    return sum(out.values())


@register_law("float_valve")
def _law_float_valve(net, ln, n):
    u = net.head_q(ln["sensor"])
    for i in range(ln["stages"]):
        ln["z"][i] = lag_step(ln["z"][i], u, ln["k"])
        u = ln["z"][i]
    q = clamp(qmul(ln["gain"], ln["setpoint"] - u), 0, ln["qmax"])
    return qmul(q, net.dt_q)


@register_law("tip")
def _law_tip(net, ln, n):
    """🎍 鹿威し（添水）: 筒の体積が v_tip に届いた瞬間に傾き、dump_steps ステップのあいだ中身を落とす
    （dump_amount が None なら一度に全部、あれば 1 ステップにその量まで）。そのあいだに注ぎ込んだ分も一緒に落ちる。
    傾いている時間が過ぎたら戻る（戻りは力学を解かず、時間で決める）。"""
    v = net.nodes[ln["src"]]["v"]
    if ln["tipped_left"] == 0:
        if ln["dump_from"] is not None:  # 戻った直後: 前の傾きで落ちた量を締める
            ln["tips"][-1]["dumped_q"] = ln["flow_sum"] - ln["dump_from"]
            ln["tips"][-1]["left_q"] = v
            ln["dump_from"] = None
        if v < ln["v_tip"]:
            return 0
        ln["tips"].append({"t": n * net.dt, "volume_q": v, "dumped_q": None, "left_q": None})
        ln["tipped_left"] = ln["dump_steps"]
        ln["dump_from"] = ln["flow_sum"]
    ln["tipped_left"] -= 1
    return v if ln["dump_amount"] is None else min(v, ln["dump_amount"])


@register_law("siphon")
def _law_siphon(net, ln, n):
    """🍷 サイフォン（ピタゴラスの杯）: 水面が管の頂（h_crest）に届くと流れ始め、水面と出口（底より h_drop 下）の
    高さの差で Torricelli の流量。水面が吸い口（h_intake）まで下がると空気が入って止まり、また頂まで溜まるのを待つ。
    止まる判定は「このステップで吸い口まで吸い切る」こと（要求が残りを超えたら、残りだけ流して止める）。"""
    v = net.nodes[ln["src"]]["v"]
    h = net.head_q(ln["src"])
    if not ln["primed"]:
        if h < ln["crest"]:
            return 0
        ln["primed"] = True
        ln["primes"].append(n * net.dt)
    avail = max(v - ln["intake_v"], 0)
    out = qmul(ln["coef"], qsqrt(h + ln["drop"]))
    if out >= avail:
        ln["primed"] = False
        ln["breaks"].append({"t": (n + 1) * net.dt, "left_q": v - avail})
        return avail
    return out


@register_law("mariotte")
def _law_mariotte(net, ln, n):
    """🫗 マリオットびん: 空気管の下端（出口より h_air 上）が水面より下にあるあいだ、出口を押す水頭は h_air で一定
    （管の下端が大気圧になるため）。水面が下端より下がると、ふつうの Torricelli（水頭 = 水面）になる。"""
    v = net.nodes[ln["src"]]["v"]
    h = net.head_q(ln["src"])
    if h <= 0:
        return 0
    if h >= ln["h_air"]:
        out = ln["const_out"]
        ln["const_steps"] += 1
    else:
        if ln["falling_from"] is None:
            ln["falling_from"] = n * net.dt
        out = qmul(ln["coef"], qsqrt(h))
    return min(out, v)


class Network:
    def __init__(self, dt: float, fluid: str = "water", g: float = G, temp_c: float | None = None,
                 mu: float | None = None, mu_scale: float | None = None, species=None, volume_of=None,
                 outlets: str = "exact", inflows: str = "uniform", rounding: str = "q48"):
        """species を渡すと成分モード（🍷）: 各水槽が成分ごとの量（Q 整数、単位は法則に任せる。蒸留ではモル）を持ち、
        v はその合計。普通の法則が返す量は、出す側の成分の比で整数のまま分ける。成分ごとの量を返す法則
        （boil・evaporate）もある。volume_of(parts) は水位（断面一定の水槽）に使う体積 [m³]。"""
        if dt <= 0:
            raise ValueError("dt must be positive")
        if outlets not in OUTLET_MODES:
            raise ValueError(f"outlets must be one of {OUTLET_MODES}")
        # 🚿 同じ槽から弁・ポンプが何本か出るときの取り方。exact = 槽ごとに 1 ステップを厳密に解く（既定）、
        # independent = どの弁もステップ最初の量を見て別々に取る（ずれを測るために残してある）
        self.outlets = outlets
        if inflows not in INFLOW_MODES:
            raise ValueError(f"inflows must be one of {INFLOW_MODES}")
        # 🚰 槽に入ってくる量の扱い。end = 入った量はステップの終わりに足す（そのステップのあいだ減らない）、
        # uniform = 1 ステップのあいだ一定の速さで入ると見なし、そのぶん弁から多く出す
        self.inflows = inflows
        if rounding not in ROUNDING_MODES:
            raise ValueError(f"rounding must be one of {ROUNDING_MODES}")
        # 🧮 丸め方。q48 = 係数も量も 2^-48 刻み、掛け算は切り捨て（既定）。
        # carry64 = 弁だけが出ている槽（開き具合が一定）で、係数を 2^-64 刻みにし、掛け算の切り捨ての残りを次のステップへ繰り越す
        self.rounding = rounding
        self._carry: dict = {}  # 掛け算の場所ごとの、切り捨ての残り
        self._carry_plan: dict = {}  # 槽ごとの係数（管の本数が変わったら作り直す）
        self.dt = dt
        self.dt_q = to_q(dt)
        # 🧊 temp_c を渡すと温度表で μ・ρ が動く。🍯 mu（Pa·s）か mu_scale（倍率）で粘度を直接変えられる
        self.fluid = get_fluid(fluid, temp_c, mu=mu, mu_scale=mu_scale)
        self.g = g
        self.nodes: dict[str, dict] = {}
        self.links: list[dict] = []
        self.channels: dict[str, dict] = {}
        self.tracers: list[dict] = []
        self.catchers: list[dict] = []
        self._groups: tuple[int, dict] | None = None
        self.last_flow: list[int] = []  # 前のステップに各管で動いた量（Q）
        self.species = tuple(species) if species else None
        if self.species is not None and len(set(self.species)) != len(self.species):
            raise ValueError("species must be unique")
        self.volume_of = volume_of
        self.last_parts: list[dict] = []

    # --- 節点 -----------------------------------------------------------
    def _node(self, name: str, **data) -> dict:
        if name in self.nodes:
            raise ValueError(f"duplicate node {name!r}")
        self.nodes[name] = data
        return data

    def tank(self, name: str, area: float | None = None, head: float = 0.0, volume: float | None = None,
             profile=None, samples: int = 512, parts: dict | None = None, outlets: str | None = None) -> "Network":
        """水槽。area 一定か、profile で高さにより断面が変わる器（🏺 水時計のような形）。
        outlets はこの槽だけ出口の取り方を変える（"exact"・"independent"。既定は Network(outlets=...) のまま）。

        profile は ``[[h, A(h)], ...]``（間は直線）か、Python の関数 ``A(h)``（spec からは表だけ）。
        0 から h_max まで samples 個に刻んで体積 V(h) を台形で積み上げ、V → h の表（Q 整数の Cam）で
        水位を引く。h_max を超えた分は一番上の断面で伸ばす。表の刻みの粗さは結果の profile_samples で分かる。
        """
        if outlets is not None and outlets not in OUTLET_MODES:
            raise ValueError(f"outlets must be one of {OUTLET_MODES}")
        if self.species is not None:  # 🍷 成分モード: parts（成分 → 量）で中身を決める
            if profile is not None or head or volume is not None:
                raise ValueError("in species mode give parts=, not head/volume/profile")
            if area is None or area <= 0:
                raise ValueError("area must be positive")
            parts = parts or {}
            unknown = set(parts) - set(self.species)
            if unknown:
                raise ValueError(f"unknown species: {sorted(unknown)}")
            pq = {sp: to_q(parts.get(sp, 0.0)) for sp in self.species}
            if any(v < 0 for v in pq.values()):
                raise ValueError("parts must be >= 0")
            self._node(name, kind="tank", area=area, area_q=to_q(area), v=sum(pq.values()), profile=None, parts=pq)
            return self
        if parts is not None:
            raise ValueError("parts= needs Network(species=...)")
        if profile is None:
            if area is None or area <= 0:
                raise ValueError("area must be positive")
            v = head * area if volume is None else volume
            self._node(name, kind="tank", area=area, area_q=to_q(area), v=to_q(v), profile=None)
            self.nodes[name]["outlets"] = outlets
            return self
        if samples < 2:
            raise ValueError("samples must be >= 2")
        if callable(profile):
            if area is None or area <= 0:
                raise ValueError("a callable profile needs area=h_max (the height the table covers)")
            h_max, area_at = float(area), profile
            knots = []
        else:
            pts = sorted((float(h), float(a)) for h, a in profile)
            if len(pts) < 2 or pts[0][0] != 0.0 or any(a <= 0 for _h, a in pts):
                raise ValueError("profile table must start at h = 0 and have positive areas")
            h_max = pts[-1][0]
            knots = [h for h, _a in pts]

            def area_at(h: float, _pts=pts) -> float:
                if h <= _pts[0][0]:
                    return _pts[0][1]
                if h >= _pts[-1][0]:
                    return _pts[-1][1]
                i = bisect.bisect_left([x for x, _ in _pts], h)
                (x0, a0), (x1, a1) = _pts[i - 1], _pts[i]
                return a0 + (a1 - a0) * (h - x0) / (x1 - x0)
        hs = sorted(set(knots) | {h_max * i / (samples - 1) for i in range(samples)})
        table, v_acc = [(0.0, 0.0)], 0.0
        for h0, h1 in zip(hs, hs[1:]):
            v_acc += 0.5 * (area_at(h0) + area_at(h1)) * (h1 - h0)
            table.append((v_acc, h1))
        a_top = area_at(h_max)
        if head and volume is None:
            if head > h_max:
                volume = v_acc + (head - h_max) * a_top
            else:
                volume = from_q(Cam([(h, v) for v, h in table])(to_q(head)))
        self._node(name, kind="tank", area=None, area_q=None, v=to_q(volume or 0.0), profile="table",
                   vh=Cam(table), h_max=h_max, v_max=v_acc, a_top=a_top, area_at=area_at, samples=len(hs))
        return self

    def channel(self, prefix: str, dst: str, *, length: float, cells: int, slope: float, section,
                law: str = "manning", n_manning: float = 0.010, beta: float | None = None,
                max_depth: float | None = None, table_points: int = 256, delay: float = 0.0,
                lip: float = 0.0, bend_radius: float | None = None, spill: str | None = None,
                cd: float = 0.62, cracks: list[dict] | None = None) -> dict:
        """傾いた樋。区間セル <prefix>:0 … <prefix>:cells-1 を作り、末端を dst へつなぐ。

        各セルの流出は「セル体積 → 流量」の表（rating カム）で引く。表は等流の式から一度だけ作り、
        水深の浅い側ほど点を密にする（(i/N)² 刻み）。表の上端（半円なら縁、長方形なら max_depth）を
        超えた分は流量が頭打ちになるので、結果の overflow で知らせる。

        - lip > 0: 末端は「次の竹の縁が lip だけ高い」継ぎ目。末端セルに水がたまり、越えた分だけ流れる
        - bend_radius: 曲がり（半径 R）。外側の水面が縁を越えた分は各セルから spill の節点へあふれる
        - delay: 末端の管の遅れ（段差の自由落下などに使う）
        - cracks: 竹の割れ目 ``[{"cell": i, "area": m², "cd": 0.61}]``。そのセルの水深 h で
          Torricelli（Cd·a·√(2gh)）だけ spill の節点へ漏れる
        """
        if cells < 1 or length <= 0 or slope <= 0:
            raise ValueError("channel needs cells >= 1, length > 0, slope > 0")
        if dst not in self.nodes:
            raise ValueError(f"unknown node {dst!r}")
        if bend_radius is not None and spill is None:
            raise ValueError("a bend needs a spill node")
        if cracks and spill is None:
            raise ValueError("cracks need a spill node")
        if spill is not None and spill not in self.nodes:
            raise ValueError(f"unknown spill node {spill!r}")
        sec = section if isinstance(section, Section) else Section(section["shape"], section["size"])
        lc = length / cells
        h_cap = sec.size if sec.shape == "semicircle" else (0.05 if max_depth is None else max_depth)
        a_max = sec.max_area(h_cap)
        table = []
        for i in range(table_points + 1):
            a = a_max * (i / table_points) ** 2
            q = channel_flow(sec, sec.depth(a), slope, law, nu=self.fluid.nu, n_manning=n_manning, beta=beta,
                             g=self.g)
            table.append((a * lc, q))
        cam = Cam(table)
        names = [f"{prefix}:{i}" for i in range(cells)]
        for name in names:
            self._node(name, kind="reach", length=lc, section=sec, v=0, channel=prefix)
        for i, name in enumerate(names):
            nxt = names[i + 1] if i + 1 < cells else dst
            if nxt == dst and lip > 0:
                self._link("lip", name, nxt, delay, section=sec, lip=lip, cd=cd, max_v=0, channel=prefix,
                           width=sec.top_width(lip), releases=[])
            else:
                self._link("rating", name, nxt, delay if nxt == dst else 0.0, cam=cam, max_v=0, channel=prefix)
            if bend_radius is not None or (spill is not None and nxt == dst and lip > 0):
                # 曲がりの全セルと、縁の手前でたまるセルには縁からあふれる口を付ける
                self._link("spill", name, spill, 0.0, main=len(self.links) - 1, radius=bend_radius, rim=h_cap,
                           v_cap=to_q(a_max * lc), cd=cd, max_outer=0.0, channel=prefix)
        for crack in cracks or ():
            unknown = set(crack) - {"cell", "area", "cd"}
            if unknown or not 0 <= crack["cell"] < cells:
                raise ValueError(f"bad crack {crack!r}")
            ln = self.orifice(names[crack["cell"]], spill, area=crack["area"], cd=crack.get("cd", 0.61))
            ln.update(channel=prefix, crack=True)
        info = {"prefix": prefix, "cells": names, "length": length, "cell_length": lc, "slope": slope,
                "section": sec.to_dict(), "law": law, "n_manning": n_manning, "beta": beta, "depth_cap": h_cap,
                "volume_cap": a_max * lc, "table": table, "lip": lip, "bend_radius": bend_radius,
                "cracks": [dict(c) for c in cracks or ()]}
        self.channels[prefix] = info
        return info

    def flume(self, prefix: str, dst: str, segments: list[dict], spill: str | None = None) -> list[dict]:
        """竹をつないだ流しそうめんの台。segments を上流から順に並べ、<prefix>/<i>:<cell> を作る。

        各 segment は channel() の引数（length・cells・slope・section・law・cracks …）に加えて
        ``joint``（次の竹との継ぎ目）を持てる: ``{"drop": 段差 m}``（自由落下の遅れ）か
        ``{"lip": 縁の高さ m}``（たまってから越える）。``bend_radius`` の竹は spill へあふれうる。
        ``node_remnant``（節の削り残し）``{"at": 0〜1 の位置, "height": m}`` は竹の途中の低い縁で、
        竹を <prefix>/<i>（節の手前、末端に縁）と <prefix>/<i>b（節の後ろ）の 2 つに分けて作る。
        cracks の cell は分ける前の竹の通し番号。
        """
        if not segments:
            raise ValueError("flume needs at least one segment")
        infos = [None] * len(segments)
        nxt = dst
        for i in range(len(segments) - 1, -1, -1):  # 下流から作る（行き先の節点を先に用意する）
            seg = dict(segments[i])
            joint = seg.pop("joint", None) or {}
            node = seg.pop("node_remnant", None)
            unknown = set(joint) - {"drop", "lip"}
            if unknown:
                raise ValueError(f"unknown joint keys: {sorted(unknown)}")
            end = {"delay": fall_time(joint.get("drop", 0.0), self.g), "lip": joint.get("lip", 0.0)}
            if node is None:
                infos[i] = self.channel(f"{prefix}/{i}", nxt, spill=spill, **end, **seg)
            else:
                if set(node) - {"at", "height"} or not 0 < node["at"] < 1 or node["height"] <= 0:
                    raise ValueError(f"bad node_remnant {node!r}")
                cells, length = seg.pop("cells"), seg.pop("length")
                cracks = seg.pop("cracks", None) or []
                cut = max(1, min(cells - 1, round(cells * node["at"])))
                after = self.channel(f"{prefix}/{i}b", nxt, spill=spill, cells=cells - cut,
                                     length=length * (cells - cut) / cells, **end,
                                     cracks=[{**c, "cell": c["cell"] - cut} for c in cracks if c["cell"] >= cut],
                                     **seg)
                infos[i] = self.channel(f"{prefix}/{i}", after["cells"][0], spill=spill, cells=cut,
                                        length=length * cut / cells, lip=node["height"],
                                        cracks=[c for c in cracks if c["cell"] < cut], **seg)
                infos[i]["node_remnant"] = dict(node)
                infos[i]["after_node"] = after
            infos[i]["joint"] = joint
            nxt = infos[i]["cells"][0]
        return infos

    # --- 管 -------------------------------------------------------------
    def _link(self, kind: str, src: str, dst: str, delay: float, **extra) -> dict:
        if kind not in LAWS:
            raise ValueError(f"unknown law {kind!r}; registered: {sorted(LAWS)}")
        for t in (src, dst):
            if t not in self.nodes:
                raise ValueError(f"unknown node {t!r}")
        if src == dst:
            raise ValueError("src and dst must differ")
        steps = int(round(delay / self.dt))
        if steps and self.species is not None:
            raise ValueError("delays are not supported in species mode yet")
        link = {
            "kind": kind, "src": src, "dst": dst, "delay": steps * self.dt, "delay_steps": steps,
            "line": DelayLine(steps - 1) if steps else None,
            "pending": 0, "flow_sum": 0, "max_head": 0, "first_arrival_step": None, "first_flow_step": None,
            **extra,
        }
        self.links.append(link)
        return link

    def link(self, kind: str, src: str, dst: str, *, delay: float = 0.0, **params) -> dict:
        """登録した任意の法則で管を足す（params はそのまま link に入り、法則の関数が読む）。"""
        return self._link(kind, src, dst, delay, **params)

    def pipe(self, src: str, dst: str, *, radius: float, length: float, delay: float = 0.0) -> dict:
        """細管（層流）。src の水位 h だけで決まる自由流出。"""
        if self.nodes[src]["kind"] != "tank":
            raise ValueError("pipe needs a tank as src")
        c = math.pi * radius ** 4 * self.g / (8.0 * self.fluid.nu * length)  # [m²/s]
        area = self.nodes[src]["area"]
        tau = None if area is None else area / c
        return self._link("pipe", src, dst, delay, radius=radius, length=length, conductance=c, c_q=to_q(c),
                          tau=tau, k=None if tau is None else lag_gain(self.dt, tau))

    def orifice(self, src: str, dst: str, *, area: float, cd: float = 0.61, delay: float = 0.0,
                viscous: bool = False, delta: float | None = None) -> dict:
        """小孔（Torricelli）。viscous=True で、粘度で下がる流量係数（fluids.orifice_cd）を毎ステップ解く。
        省くと Cd 一定のまま（粘度によらない）。どちらでも結果に Re と、Re の低さの印 viscous_regime を残す。"""
        from fluids import ORIFICE_DELTA
        coef = to_q(cd * area * math.sqrt(2.0 * self.g) * self.dt)
        return self._link("orifice", src, dst, delay, coef=coef, orifice_area=area, cd=cd, viscous=viscous,
                          delta=ORIFICE_DELTA if delta is None else delta, diameter=math.sqrt(4.0 * area / math.pi),
                          min_cd=None)

    def weir(self, src: str, dst: str, *, shape: str = "rect", width: float = 0.02, angle_deg: float = 90.0,
             cd: float | None = None, delay: float = 0.0) -> dict:
        s2g = math.sqrt(2.0 * self.g)
        if shape == "rect":
            cd = 0.62 if cd is None else cd
            k = cd * (2.0 / 3.0) * width * s2g
        elif shape == "v":
            cd = 0.58 if cd is None else cd
            k = cd * (8.0 / 15.0) * math.tan(math.radians(angle_deg) / 2.0) * s2g
        else:
            raise ValueError("weir shape must be 'rect' or 'v'")
        return self._link("weir_" + shape, src, dst, delay, k=k, coef=to_q(k * self.dt), cd=cd)

    def cam(self, src: str, dst: str, *, table, delay: float = 0.0) -> dict:
        """[[水位, 流量], ...] の表で決まる流出（特殊形状のスリット）。"""
        return self._link("cam", src, dst, delay, cam=Cam(table), table=[list(map(float, p)) for p in table])

    def rain(self, prefix: str, intensity: float, src: str, *, start: float = 0.0, until: float | None = None) -> list[dict]:
        """🌧️ 樋 prefix の全セルに一様な雨 intensity [m/s] を降らせる（src の水槽から各セルへのポンプ）。

        1 セルの受ける流量 = intensity × 平面の幅（fluids.plan_width）× セル長。ポンプなので start〜until の間だけ降り、
        src が空になれば按分で止まる。戻り値は作ったポンプの並び（link に rain=True・channel=prefix を付ける）。
        """
        if prefix not in self.channels:
            raise ValueError(f"unknown channel {prefix!r}")
        if intensity < 0:
            raise ValueError("intensity must be >= 0")
        ch = self.channels[prefix]
        sec = self.nodes[ch["cells"][0]]["section"]
        rate = intensity * plan_width(sec) * ch["cell_length"]
        out = []
        for cell in ch["cells"]:
            ln = self.pump(src, cell, rate=rate, start=start, until=until)
            ln.update(rain=True, channel=prefix, intensity=intensity)
            out.append(ln)
        return out

    def pump(self, src: str, dst: str, *, rate: float, start: float = 0.0, until: float | None = None,
             delay: float = 0.0) -> dict:
        """一定流量。start から until まで動く。"""
        return self._link("pump", src, dst, delay, amount=qmul(to_q(rate), self.dt_q), rate=rate,
                          start_step=int(round(start / self.dt)),
                          until_step=None if until is None else int(round(until / self.dt)))

    def valve(self, src: str, dst: str, *, tau: float, opening: float = 1.0, sensor: str | None = None,
              cam=None, stages: int = 0, link_tau: float = 1.0, delay: float = 0.0) -> dict:
        """弁つきの出口（MONIAC の弁）。全開の時定数 tau、開き具合は opening（0..1）か、
        sensor の水位を cam（[[水位, 開き具合], ...]）で読み替えたもの。stages 段の遅れ（各 link_tau）は
        浮き → ワイヤー → 弁のつながりの遅れ。"""
        if self.nodes[src]["kind"] != "tank" or self.nodes[src]["profile"] is not None:
            raise ValueError("valve needs a constant-area tank as src")
        if not 0.0 <= opening <= 1.0:
            raise ValueError("opening must be within 0..1")
        if (cam is None) != (sensor is None):
            raise ValueError("give both sensor and cam, or neither")
        if sensor is not None and sensor not in self.nodes:
            raise ValueError(f"unknown sensor node {sensor!r}")
        h0 = self.head_q(sensor) if sensor else 0
        return self._link("valve", src, dst, delay, tau=tau, k=lag_gain(self.dt, tau), opening=to_q(opening),
                          sensor=sensor, cam=None if cam is None else Cam(cam),
                          table=None if cam is None else [list(map(float, p)) for p in cam],
                          stages=stages, k_link=lag_gain(self.dt, link_tau), z=[h0] * stages, last_opening=0)

    def boil(self, src: str, dst: str | None = None, *, power: float, start: float = 0.0, until: float | None = None,
             cuts=None) -> dict:
        """🔥 釜で熱 power [W] を加えて蒸気を dst へ（凝縮して受け器にたまる）。成分 ethanol・water（モル）が要る。

        cuts で蒸留の切り替え（✂️ 前留・本留・後留）: ``[{"to": 受け器, "volume": m³}, {"to": …, "abv_below": %}, {"to": …}]``。
        今の受け器の体積が volume に届くか、流れ出る留液の度数が abv_below を切ったら次へ移る。最後の受け器で止まる。
        """
        if self.species is None or not {"ethanol", "water"} <= set(self.species):
            raise ValueError("boil needs Network(species=('ethanol', 'water'))")
        if power <= 0:
            raise ValueError("power must be positive")
        cuts = [dict(c) for c in (cuts or [])]
        for c in cuts:
            if c.get("to") not in self.nodes:
                raise ValueError(f"unknown cut receiver {c.get('to')!r}")
            if set(c) - {"to", "volume", "abv_below"}:
                raise ValueError("cut keys are to, volume, abv_below")
        if cuts:
            cuts[0]["started"] = start
            dst = cuts[0]["to"]
        if dst is None:
            raise ValueError("give dst or cuts")
        return self._link("boil", src, dst, 0.0, power=power, start_step=int(round(start / self.dt)),
                          until_step=None if until is None else int(round(until / self.dt)), last_x=None, last_y=None,
                          last_abv=None, cuts=cuts, cut_index=0)

    def evaporate(self, src: str, dst: str, *, rates: dict) -> dict:
        """👼 成分ごとの 1 次の抜け（rates = {成分: k [1/s]}）。樽の木を通って空へ逃げる天使の取り分など。"""
        if self.species is None:
            raise ValueError("evaporate needs Network(species=...)")
        unknown = set(rates) - set(self.species)
        if unknown or any(k < 0 for k in rates.values()):
            raise ValueError("rates must be known species with k >= 0")
        gains = {sp: (lag_gain(self.dt, 1.0 / k) if k > 0 else 0) for sp, k in rates.items()}
        return self._link("evaporate", src, dst, 0.0, rates=dict(rates), gains=gains)

    def float_valve(self, src: str, dst: str, *, sensor: str, setpoint: float, gain: float, qmax: float,
                    stages: int = 1, tau: float = 1.0, delay: float = 0.0) -> dict:
        if sensor not in self.nodes:
            raise ValueError(f"unknown sensor node {sensor!r}")
        if stages < 0:
            raise ValueError("stages must be >= 0")
        h0 = self.head_q(sensor)
        return self._link("float_valve", src, dst, delay, sensor=sensor, setpoint=to_q(setpoint),
                          gain=to_q(gain), qmax=to_q(qmax), k=lag_gain(self.dt, tau), z=[h0] * stages,
                          stages=stages, tau=tau)

    def tip(self, src: str, dst: str, *, v_tip: float, t_dump: float | None = None, dump_rate: float | None = None,
            delay: float = 0.0) -> dict:
        """🎍 鹿威しの筒。v_tip [m³] 溜まると傾き、t_dump 秒（省くと 1 ステップ）のあいだ dump_rate [m³/s]
        （省くと一度に全部）で落とす。"""
        if v_tip <= 0:
            raise ValueError("v_tip must be positive")
        dump_steps = 1 if t_dump is None else max(1, int(round(t_dump / self.dt)))
        if dump_rate is not None and dump_rate <= 0:
            raise ValueError("dump_rate must be positive")
        return self._link("tip", src, dst, delay, v_tip=to_q(v_tip), dump_steps=dump_steps, t_dump=dump_steps * self.dt,
                          dump_amount=None if dump_rate is None else qmul(to_q(dump_rate), self.dt_q),
                          dump_rate=dump_rate, tipped_left=0, dump_from=None, tips=[])

    def siphon(self, src: str, dst: str, *, h_crest: float, h_drop: float = 0.0, area: float, cd: float = 0.6,
               h_intake: float = 0.0, delay: float = 0.0) -> dict:
        """🍷 サイフォン。頂 h_crest、出口は底より h_drop 下、管の断面 area、Cd は損失をまとめた係数（既定 0.6）、
        吸い口の高さ h_intake（既定 0 = 底まで吸い切る）。src は水槽。"""
        if self.nodes[src]["kind"] != "tank":
            raise ValueError("siphon needs a tank as src")
        if h_crest <= h_intake or h_intake < 0 or h_drop < 0 or area <= 0:
            raise ValueError("siphon needs h_crest > h_intake >= 0, h_drop >= 0, area > 0")
        return self._link("siphon", src, dst, delay, crest=to_q(h_crest), drop=to_q(h_drop),
                          intake_v=qmul(to_q(h_intake), self.nodes[src]["area_q"]), h_crest=h_crest, h_drop=h_drop,
                          h_intake=h_intake, tube_area=area, cd=cd,
                          coef=to_q(cd * area * math.sqrt(2.0 * self.g) * self.dt), primed=False, primes=[], breaks=[])

    def mariotte(self, src: str, dst: str, *, area: float, h_air: float, cd: float = 0.61, delay: float = 0.0) -> dict:
        """🫗 マリオットびんの出口。小孔の断面 area、空気管の下端は出口より h_air 上。src は水槽。"""
        if self.nodes[src]["kind"] != "tank":
            raise ValueError("mariotte needs a tank as src")
        if h_air <= 0 or area <= 0:
            raise ValueError("mariotte needs h_air > 0 and area > 0")
        coef = to_q(cd * area * math.sqrt(2.0 * self.g) * self.dt)
        return self._link("mariotte", src, dst, delay, coef=coef, orifice_area=area, cd=cd, h_air=to_q(h_air),
                          h_air_m=h_air, const_out=qmul(coef, qsqrt(to_q(h_air))), const_steps=0, falling_from=None)

    # --- 運ばれる物 -------------------------------------------------------
    def tracer(self, name: str, node: str, release: float = 0.0, stop=(), size: float = 0.0, length: float = 0.0,
               strands: int = 1, density: float = 1100.0, friction: float = 0.3) -> dict:
        """運ばれる物を node に release 時刻で放す。stop（節点名か、その並び）に入るか、出口のない節点で止まる。

        - size: 太さ [m]。継ぎ目の段（lip）では越流水深がこれより浅いと引っかかって待つ
        - length: 長さ [m]。頭の後ろに、通ってきた区間を合計 length まで体として引く。体のどこかが
          越えられない縁に掛かっていれば全体が止まる。曲がりではじき出される判定も体のどのセルでも受ける
        - strands: 束の本数（一口分）。縁に積もる体積に効く
        - density・friction: 茹でた麺の密度 [kg/m³] と竹との摩擦係数（どちらも仮定の値）。積もった麺を
          押し流す力のつり合いに使う
        """
        stops = (stop,) if isinstance(stop, str) else tuple(stop)
        for nd in (node, *stops):
            if nd not in self.nodes:
                raise ValueError(f"unknown node {nd!r}")
        if length < 0 or strands < 1:
            raise ValueError("length must be >= 0 and strands >= 1")
        tr = {"name": name, "node": node, "release_step": int(round(release / self.dt)), "state": "waiting",
              "credit": 0, "pipe_left": 0, "pipe_dst": None, "pipe_t": None, "pipe_f": 0.0, "path": [],
              "arrival": None, "stop": stops, "size": size, "stuck_steps": 0, "stuck_at": None,
              "held": False, "held_link": None, "length": length, "strands": strands, "density": density,
              "friction": friction, "cells": [], "flushed": set(), "until": None, "end": None,
              "held_since": None}
        self.tracers.append(tr)
        return tr

    def catcher(self, name: str, node: str, *, t_grab: float | None = None, grab=(0.15, 0.6),
                cooldown: float = 0.0, seed: int = 0) -> dict:
        """🥢 取り手（箸）。node（区間セル）に入った麺を、入ってから t_grab 秒たったところで取る。

        - 麺がその前にセルを出てしまえば取り逃し（fast）。そのとき手がふさがっていても取り逃し（busy）
        - t_grab（手を伸ばして挟むまでの時間）は麺ごとに grab = (最短, 最長) の一様乱数（seed で再現）。
          t_grab を渡すと全部その値。取ったら cooldown 秒は次を取れない（口へ運ぶ時間）
        - 取った麺は水から外れるだけで、水の量は変えない
        """
        if node not in self.nodes:
            raise ValueError(f"unknown node {node!r}")
        lo, hi = grab
        if t_grab is None and not 0 <= lo <= hi:
            raise ValueError("grab must be (min, max) with 0 <= min <= max")
        rng = random.Random(seed)

        def draw(tr):
            return t_grab if t_grab is not None else rng.uniform(lo, hi)

        c = {"name": name, "node": node, "t_grab": t_grab, "grab": (lo, hi), "cooldown": cooldown, "seed": seed,
             "draw": draw, "free_at": -math.inf, "caught": [], "missed": []}
        self.catchers.append(c)
        return c

    # --- 状態 -----------------------------------------------------------
    def head_q(self, name: str) -> int:
        nd = self.nodes[name]
        if nd["kind"] == "tank":
            if "parts" in nd:  # 🍷 成分モード: 体積は volume_of（なければ量そのものを体積とみなす）
                if self.volume_of is None:
                    return qdiv(nd["v"], nd["area_q"])
                return to_q(self.volume_of({k: from_q(v) for k, v in nd["parts"].items()}) / nd["area"])
            if nd["profile"] is None:
                return qdiv(nd["v"], nd["area_q"])
            v = nd["v"]
            if v <= 0:
                return 0
            if from_q(v) > nd["v_max"]:  # 表の上は一番上の断面で伸ばす
                return to_q(nd["h_max"] + (from_q(v) - nd["v_max"]) / nd["a_top"])
            return nd["vh"](v)
        return to_q(nd["section"].depth(from_q(nd["v"]) / nd["length"]))

    def head(self, name: str) -> float:
        return from_q(self.head_q(name))

    def total_q(self) -> int:
        inside = sum(nd["v"] for nd in self.nodes.values())
        transit = sum((ln["line"].in_transit if ln["line"] else 0) + ln["pending"] for ln in self.links)
        return inside + transit

    def _by_src(self) -> dict[str, list[int]]:
        if self._groups is None or self._groups[0] != len(self.links):
            groups: dict[str, list[int]] = {}
            for i, ln in enumerate(self.links):
                groups.setdefault(ln["src"], []).append(i)
            self._groups = (len(self.links), groups)
        return self._groups[1]

    def _cmul(self, site, v: int, coef: int) -> int:
        """🧮 量 × 係数（2^-64 刻み）。切り捨ての残りを場所ごとに持ち、次のステップの同じ掛け算に足す。"""
        acc = self._carry.get(site, 0) + v * coef
        self._carry[site] = acc & ((1 << CARRY_BITS) - 1)
        return acc >> CARRY_BITS

    def _carry_tank(self, src: str, idx: list[int]):
        """🧮 carry64 が効く槽なら (弁の番号, Φ, 分け前, 入った量のうち同じステップに出る割合) を返す。効かなければ None。

        効くのは、出ている管が弁だけで、開き具合が一定（浮きで動かない）で、槽ごとに厳密に解く設定の槽。
        L = Σ 開き具合·dt/τ は分数で厳密に、Φ = 1 − e^(−L) は整数の級数で出す（float を通さない）。"""
        if self.rounding != "carry64" or self.species is not None:
            return None
        key = (src, len(self.links))
        if key not in self._carry_plan:
            links = [self.links[i] for i in idx]
            plan = None
            if (all(ln["kind"] == "valve" and ln["cam"] is None for ln in links)
                    and (self.nodes[src].get("outlets") or self.outlets) == "exact"):
                weights = [Fraction(ln["opening"], ONE) * Fraction(self.dt) / Fraction(ln["tau"]) for ln in links]
                scale = sum(weights)
                if scale > 0:
                    one = 1 << CARRY_BITS
                    phi = one - exp_neg_q(scale, CARRY_BITS)
                    plan = (list(idx), phi, [round(w / scale * one) for w in weights[:-1]],
                            round((1 - Fraction(phi, one) / scale) * one))
            self._carry_plan[key] = plan
        return self._carry_plan[key]

    def _carry_outlets(self, req: list[int]) -> None:
        """🧮 carry64 が効く槽の、今までどおりの量（V·Φ）を出し直す。弁どうしは速さの比で分け、端数は最後の弁へ。"""
        for src, idx in self._by_src().items():
            plan = self._carry_tank(src, idx)
            if plan is None:
                continue
            valves, phi, shares, _ = plan
            have = self.nodes[src]["v"]
            rest = total = min(max(self._cmul(("base", src), have, phi), 0), have)
            for i, share in zip(valves[:-1], shares):
                req[i] = min(self._cmul(("share", i), total, share), rest)
                rest -= req[i]
            req[valves[-1]] = rest

    def _exact_outlets(self, req: list[int]) -> None:
        """🚿 弁とポンプだけが出ている槽を、1 ステップぶん厳密に解く（dV/dt = −K·V − q）。

        弁はどれもステップ最初の量 V の k_i·開き具合 を取りたがる。別々に取ると、合計は厳密な 1 − e^(−K·dt) より
        少し多い（n 本なら刻みの 1 次、(1 − 1/n) に比例して速く空く）。ここで弁の合計を (V + P/L)·Φ − P に置き直す:
        L = K·dt = Σ 開き具合_i·dt/τ_i、Φ = 1 − e^(−L)、P = このステップのポンプの量（q·dt）。ポンプは P のまま。
        弁どうしは速さ（開き具合/τ）の比で分け、端数は最後の弁へ（合計は整数で一致）。
        弁が 1 本でポンプも無い槽は法則のまま（元から厳密）。槽ごとに tank(outlets=...) で切り替えられる。
        槽と槽のあいだ（上の槽から入る量を同じステップの中で知らない）のずれは、ここでは直らない。
        """
        for src, idx in self._by_src().items():
            node = self.nodes[src]
            if (node.get("outlets") or self.outlets) != "exact" or self._carry_tank(src, idx) is not None:
                continue
            if any(self.links[i]["kind"] not in ("valve", "pump") for i in idx):
                continue
            valves = [i for i in idx if self.links[i]["kind"] == "valve" and req[i] > 0]
            pumped = sum(req[i] for i in idx if self.links[i]["kind"] == "pump")
            if not valves or (len(valves) == 1 and not pumped):
                continue
            # 速さ × dt = 開き具合 × dt/τ（dt/τ = −ln(1 − k)。k は弁を作ったときの厳密離散化ゲイン）
            weights = [from_q(self.links[i]["last_opening"]) * -math.log1p(-from_q(self.links[i]["k"])) for i in valves]
            scale = sum(weights)
            if node["v"] <= 0 or scale <= 0.0:
                continue
            phi = to_q(-math.expm1(-scale))
            if pumped:
                total = qmul(node["v"] + to_q(from_q(pumped) / scale), phi) - pumped
            else:
                total = qmul(node["v"], phi)
            rest = total = max(total, 0)
            for i, w in zip(valves[:-1], weights):
                req[i] = min(qmul(total, to_q(w / scale)), rest)
                rest -= req[i]
            req[valves[-1]] = rest

    def _uniform_inflows(self, req: list[int]) -> None:
        """🚰 このステップに届く量 I が、1 ステップのあいだ一定の速さで入ると見なす（dV/dt = −K·V + (I − P)/dt）。

        入った量をステップの終わりに足すだけだと、そのステップのあいだ減らないので、下の槽は dt/(2·τ) ほど多くなる。
        ここで、弁とポンプだけが出ている槽の弁の合計を V·Φ + (I − P)·(1 − Φ/L) まで増やす
        （L = K·dt、Φ = 1 − e^(−L)、P = 出るポンプの量）。req はもう全部決まっている（出る量はステップ最初の量まで）ので、
        I は実際に届く量。増やす分は I までなので槽は負にならず、総量は整数のまま変わらない。
        増やした分が下の槽に届いても、下の槽の I には数えない（刻みの 2 次の差）。
        上がポンプ（本当に一定の速さ）なら厳密、上が弁（減りながら入る）なら刻みの 2 次まで小さくなる。
        """
        arrivals: dict[str, int] = {}
        for i, ln in enumerate(self.links):
            arrive = req[i] if ln["line"] is None else ln["pending"]
            if arrive:
                arrivals[ln["dst"]] = arrivals.get(ln["dst"], 0) + arrive
        base = list(req)
        for src, idx in self._by_src().items():
            inflow = arrivals.get(src, 0)
            if inflow <= 0 or any(self.links[i]["kind"] not in ("valve", "pump") for i in idx):
                continue
            plan = self._carry_tank(src, idx)
            if plan is not None:  # 🧮 増やす分も 2^-64 刻みの係数と繰り越しで出す
                valves, _, shares, through = plan
                rest = extra = min(max(self._cmul(("extra", src), inflow, through), 0), inflow)
                for i, share in zip(valves[:-1], shares):
                    part = min(self._cmul(("xshare", i), extra, share), rest)
                    req[i] += part
                    rest -= part
                req[valves[-1]] += rest
                continue
            valves = [i for i in idx if self.links[i]["kind"] == "valve"]
            weights = [from_q(self.links[i]["last_opening"]) * -math.log1p(-from_q(self.links[i]["k"])) for i in valves]
            scale = sum(weights)
            if scale <= 0.0:
                continue
            node = self.nodes[src]
            through = to_q(1.0 + math.expm1(-scale) / scale)  # 入った量のうち、同じステップの中で弁から出る割合
            pumped = sum(base[i] for i in idx if self.links[i]["kind"] == "pump")
            if pumped and (node.get("outlets") or self.outlets) == "exact":
                target = qmul(node["v"], to_q(-math.expm1(-scale))) + qmul(inflow - pumped, through)
                extra = max(target, 0) - sum(base[i] for i in valves)
            else:
                extra = qmul(inflow, through)
            rest = extra = min(max(extra, 0), inflow)
            for i, w in zip(valves[:-1], weights):
                part = min(qmul(extra, to_q(w / scale)), rest)
                req[i] += part
                rest -= part
            req[valves[-1]] += rest

    # --- 1 ステップ -----------------------------------------------------
    def step(self, n: int) -> list[int]:
        if len(self.last_flow) != len(self.links):
            self.last_flow = [0] * len(self.links)
        req = [max(LAWS[ln["kind"]](self, ln, n), 0) for ln in self.links]
        if self.rounding == "carry64":
            self._carry_outlets(req)
        self._exact_outlets(req)
        # 出す側の水が足りなければ、その節点から出る管で整数のまま按分する
        for src, idx in self._by_src().items():
            want = sum(req[i] for i in idx)
            have = self.nodes[src]["v"]
            if want > have:
                for i in idx:
                    req[i] = req[i] * have // want
                # 床で余った端数は最大の管へ（足りない水槽は空になるまで出し切る）
                req[max(idx, key=lambda i: req[i])] += have - sum(req[i] for i in idx)
        if self.species is not None:
            return self._step_species(n, req)
        if self.inflows == "uniform":
            self._uniform_inflows(req)
        for i, ln in enumerate(self.links):
            self.nodes[ln["src"]]["v"] -= req[i]
            ln["flow_sum"] += req[i]
            if req[i] and ln["first_flow_step"] is None:
                ln["first_flow_step"] = n
        for i, ln in enumerate(self.links):
            if ln["line"] is None:
                arrive = req[i]
            else:  # 今届く分を足し、今出た分は管の中へ（delay ステップ後に届く）
                arrive = ln["pending"]
                ln["pending"] = ln["line"].push(req[i])
            self.nodes[ln["dst"]]["v"] += arrive
            if arrive and ln["first_arrival_step"] is None:
                ln["first_arrival_step"] = n
        if self.tracers:
            self._move_tracers(n, req)
        self.last_flow = req
        return req

    def _step_species(self, n: int, req: list[int]) -> list[int]:
        """🍷 成分モードの 1 ステップ。法則が成分ごとの量（ln["_want"]）を置いていればそれを、なければ
        合計 req を出す側の成分の比で床で分ける。成分ごとに足りなければ整数のまま按分する。
        移した成分をそのまま足すので、成分ごとの総量が 1 ulp も変わらない。"""
        want = []
        for i, ln in enumerate(self.links):
            explicit = ln.pop("_want", None)
            if explicit is not None:
                want.append({sp: max(int(explicit.get(sp, 0)), 0) for sp in self.species})
                continue
            src = self.nodes[ln["src"]]
            v = src["v"]
            want.append({sp: (req[i] * q // v if v > 0 else 0) for sp, q in src["parts"].items()})
        for src, idx in self._by_src().items():
            parts = self.nodes[src]["parts"]
            for sp in self.species:
                need = sum(want[i][sp] for i in idx)
                if need > parts[sp]:
                    for i in idx:
                        want[i][sp] = want[i][sp] * parts[sp] // need
        moved = []
        for i, ln in enumerate(self.links):
            src, dst = self.nodes[ln["src"]], self.nodes[ln["dst"]]
            total = 0
            for sp, q in want[i].items():
                src["parts"][sp] -= q
                dst["parts"][sp] += q
                total += q
            src["v"] -= total
            dst["v"] += total
            ln["flow_sum"] += total
            ln.setdefault("parts_sum", {sp: 0 for sp in self.species})
            for sp, q in want[i].items():
                ln["parts_sum"][sp] += q
            if total and ln["first_flow_step"] is None:
                ln["first_flow_step"] = n
                ln["first_arrival_step"] = n
            moved.append(total)
        self.last_flow = moved
        self.last_parts = want
        return moved

    def parts_total(self) -> dict[str, int]:
        """成分ごとの総量（Q）。"""
        return {sp: sum(nd["parts"][sp] for nd in self.nodes.values()) for sp in self.species}

    def _move_tracers(self, n: int, req: list[int]) -> None:
        """先入れ先出し。先にいた水（credit）を流出で減らし、0 を切ったステップ内の位置 f まで時刻を補間する。

        1. 縁に積もった麺が水に押し流されるか（力のつり合い）を先に見る
        2. 体のどこかが曲がりの縁の外へはじき出されるか
        3. 先にいた水が出たら、頭の行き先の管と、体がまたいでいる管を全部越えられるか（越えられなければ待つ）
        4. 取り手（箸）の区間なら、入ってから t_grab たったところで取れるか
        """
        by_src = self._by_src()
        now = (n + 1) * self.dt
        self._flush_piles(n)
        for tr in self.tracers:
            if tr["state"] == "waiting":
                if n < tr["release_step"]:
                    continue
                # 放した時点でその節点にいた水（このステップの流入を含む）が「先にいた水」
                tr["state"] = "node"
                tr["credit"] = self.nodes[tr["node"]]["v"] + sum(req[i] for i in by_src.get(tr["node"], []))
                tr["path"].append({"node": tr["node"], "t": n * self.dt})
                self._note_cell(tr, tr["node"], n * self.dt)
            if tr["state"] == "pipe":
                tr["pipe_left"] -= 1
                if tr["pipe_left"] <= 0:
                    self._enter(tr, tr["pipe_dst"], tr["pipe_t"], tr["pipe_f"], req)
                continue
            if tr["state"] != "node":
                continue
            if self._eject(tr, now):
                continue
            outs = by_src.get(tr["node"], [])
            if not outs:
                tr["state"], tr["arrival"] = "done", tr["path"][-1]["t"]
                continue
            out_now = sum(req[i] for i in outs)
            tr["credit"] -= out_now
            if tr["credit"] <= 0 and out_now > 0:
                ln = self.links[max(outs, key=lambda i: req[i])]
                blocked = self._blocking(tr, ln)
                if blocked is not None:
                    tr["stuck_steps"] += 1  # 引っかかり: 先にいた水は出たが、どこかの縁を越えられない
                    if not tr["held"]:
                        tr["held_since"] = now
                    tr["stuck_at"], tr["held"], tr["held_link"] = blocked["src"], True, blocked
                    continue
                # このステップのうち、越えた後に残る割合（引っかかりから外れた直後は 0）
                f = 0.0 if tr["held"] else -tr["credit"] / out_now
                tr["held"], tr["held_link"] = False, None
                t_cross = (n + 1 - f) * self.dt
                if self._catch_on_leave(tr, t_cross):
                    continue
                if ln["delay_steps"]:
                    tr.update(state="pipe", pipe_left=ln["delay_steps"], pipe_dst=ln["dst"],
                              pipe_t=t_cross + ln["delay"], pipe_f=f)
                else:
                    self._enter(tr, ln["dst"], t_cross, f, req)
            self._catch_waiting(tr, now)

    # --- 体（長さ）・縁・はじき出し・箸 -------------------------------------
    def _note_cell(self, tr: dict, node: str, t: float) -> None:
        if self.nodes[node]["kind"] == "reach":
            tr["cells"].append(node)
        for c in self.catchers:
            if c["node"] == node:
                tr["station"] = {"catcher": c, "t_in": t, "t_grab": c["draw"](tr), "done": False}

    def body(self, tr: dict) -> list[str]:
        """頭から後ろへ、通ってきた区間を合計 length になるまで（頭の区間も 1 つ分に数える）。"""
        head = tr["node"]
        if self.nodes[head]["kind"] != "reach" or not tr["cells"] or tr["cells"][-1] != head:
            return [head]
        out, total = [], 0.0
        for c in reversed(tr["cells"]):
            out.append(c)
            total += self.nodes[c]["length"]
            if total >= tr["length"] - 1e-12:
                break
        return out

    def _pair_link(self, src: str, dst: str):
        for i in self._by_src().get(src, []):
            if self.links[i]["dst"] == dst and self.links[i]["kind"] in PASSABLE:
                return self.links[i]
        return None

    def _blocking(self, tr: dict, head_link: dict):
        """越えられない管（頭の行き先、または体がまたいでいる縁）。全部越えられれば None。"""
        candidates = [head_link]
        body = self.body(tr)
        for down, up in zip(body, body[1:]):
            ln = self._pair_link(up, down)
            if ln is not None:
                candidates.append(ln)
        for ln in candidates:
            ok = PASSABLE.get(ln["kind"])
            if ok is not None and not ok(self, ln, tr):
                return ln
        return None

    def _eject(self, tr: dict, now: float) -> bool:
        by_src = self._by_src()
        for cell in self.body(tr):
            for i in by_src.get(cell, []):
                ln = self.links[i]
                rule = EJECT.get(ln["kind"])
                if rule is not None and rule(self, ln, tr):
                    tr.update(state="spilled", node=ln["dst"], until=now, end="spilled", held=False, held_link=None)
                    tr["spilled_at"] = cell
                    tr["path"].append({"node": ln["dst"], "t": now})
                    return True
        return False

    def _flush_piles(self, n: int) -> None:
        """縁に積もった麺（たまり）を、水の押す力が摩擦を超えたら全部まとめて押し流す。

        押す力 F = ρ·g·w·Δ·(H − lip − Δ/2)（たまりの上流面にかかる静水圧。H はそのセルの水深）
        支える力 R = Σ μ·(ρ_麺 − ρ)·g·(本数·πd²/4)·L（水中の重さ × 摩擦）。F ≥ R で流れ出す。
        """
        for ln in self.links:
            if ln["kind"] != "lip":
                continue
            members = [tr for tr in pile_members(self, ln) if id(ln) not in tr["flushed"]]
            if not members:
                continue
            delta = pile_height(self, ln)
            h = from_q(self.head_q(ln["src"]))
            rho = self.fluid.rho
            push = rho * self.g * ln["width"] * delta * (h - ln["lip"] - delta / 2.0)
            hold = sum(tr["friction"] * (tr["density"] - rho) * self.g * tr["strands"] * math.pi * tr["size"] ** 2
                       / 4.0 * tr["length"] for tr in members)
            if push >= hold:
                for tr in members:
                    tr["flushed"].add(id(ln))
                ln["releases"].append({"t": (n + 1) * self.dt, "members": [tr["name"] for tr in members],
                                       "depth": h, "pile": delta, "push": push, "hold": hold})

    def _catch_on_leave(self, tr: dict, t_cross: float) -> bool:
        st = tr.get("station")
        if not st or st["done"] or st["catcher"]["node"] != tr["node"]:
            return False
        if t_cross - st["t_in"] >= st["t_grab"]:
            return self._try_catch(tr, st)
        st["done"] = True
        st["catcher"]["missed"].append({"tracer": tr["name"], "t_in": st["t_in"], "residence": t_cross - st["t_in"],
                                        "t_grab": st["t_grab"], "why": "fast"})
        return False

    def _catch_waiting(self, tr: dict, now: float) -> None:
        st = tr.get("station")
        if tr["state"] != "node" or not st or st["done"] or st["catcher"]["node"] != tr["node"]:
            return
        if now - st["t_in"] >= st["t_grab"]:
            self._try_catch(tr, st)

    def _try_catch(self, tr: dict, st: dict) -> bool:
        c, t_c = st["catcher"], st["t_in"] + st["t_grab"]
        st["done"] = True
        if c["free_at"] <= t_c + 1e-12:
            c["free_at"] = t_c + c["cooldown"]
            c["caught"].append({"tracer": tr["name"], "t_in": st["t_in"], "t_catch": t_c, "t_grab": st["t_grab"]})
            tr.update(state="caught", until=t_c, end="caught", held=False, held_link=None)
            return True
        c["missed"].append({"tracer": tr["name"], "t_in": st["t_in"], "t_grab": st["t_grab"], "why": "busy",
                            "free_at": c["free_at"]})
        return False

    def _enter(self, tr: dict, node: str, t: float, f: float, req: list[int]) -> None:
        tr["node"], tr["state"] = node, "node"
        outs = self._by_src().get(node, [])
        # 入った後のステップの残り f の間に、すでに先へ出ていった分を差し引く
        tr["credit"] = self.nodes[node]["v"] - int(f * sum(req[i] for i in outs))
        tr["path"].append({"node": node, "t": t})
        self._note_cell(tr, node, t)
        if node in tr["stop"] or not outs:
            tr["state"], tr["arrival"] = "done", t

    def run(self, t_end: float, *, sample_every: int | None = None, trace=()) -> dict:
        """t_end まで流す。series はプロッター記録（抽出）、trace は指定した節点の毎ステップ水位。"""
        steps = int(round(t_end / self.dt))
        if sample_every is None:
            sample_every = max(1, steps // 200)
        total0 = self.total_q()
        parts0 = self.parts_total() if self.species is not None else None
        parts_worst = {sp: 0 for sp in self.species} if self.species is not None else None
        worst = 0
        series = []
        traces = {name: [] for name in trace}
        tank_names = [k for k, nd in self.nodes.items() if nd["kind"] == "tank"]
        for n in range(steps):
            last = self.step(n)
            worst = max(worst, abs(self.total_q() - total0))
            if parts0 is not None:
                now_parts = self.parts_total()
                for sp in self.species:
                    parts_worst[sp] = max(parts_worst[sp], abs(now_parts[sp] - parts0[sp]))
            for name in trace:
                traces[name].append(self.head(name))
            if (n + 1) % sample_every == 0 or n + 1 == steps:
                row = {"t": round((n + 1) * self.dt, 9)}
                row.update({f"h:{k}": self.head(k) for k in tank_names})
                for prefix, ch in self.channels.items():
                    row[f"h:{prefix}"] = [self.head(c) for c in ch["cells"]]
                row.update({f"q:{i}": from_q(qdiv(q, self.dt_q)) for i, q in enumerate(last)})
                if self.species is not None:
                    for k in tank_names:
                        row[f"parts:{k}"] = {sp: from_q(q) for sp, q in self.nodes[k]["parts"].items()}
                bodies = {tr["name"]: self.body(tr) for tr in self.tracers
                          if tr["length"] > 0 and tr["state"] == "node"}
                if bodies:
                    row["bodies"] = bodies
                series.append(row)
        return {
            "fluid": self.fluid.key, "temp_c": self.fluid.temp_c, "dt": self.dt, "steps": steps, "series": series, "trace": traces,
            "final": series[-1] if series else {}, "total_q": total0, "conservation_error_q": worst,
            "volumes_q": {k: nd["v"] for k, nd in self.nodes.items()},
            "species": list(self.species) if self.species is not None else None,
            "parts_q": {k: dict(nd["parts"]) for k, nd in self.nodes.items() if "parts" in nd},
            "parts_error_q": parts_worst,
            "profile_samples": {k: nd["samples"] for k, nd in self.nodes.items()
                                if nd["kind"] == "tank" and nd["profile"] is not None},
            "links": [self._link_summary(i, ln) for i, ln in enumerate(self.links)],
            "channels": {p: self._channel_summary(ch) for p, ch in self.channels.items()},
            "tracers": [{"name": tr["name"], "state": tr["state"], "arrival": tr["arrival"], "path": tr["path"],
                         "size": tr["size"], "length": tr["length"], "strands": tr["strands"],
                         "stuck_time": tr["stuck_steps"] * self.dt, "stuck_at": tr["stuck_at"], "until": tr["until"],
                         "end": tr["end"], "spilled_at": tr.get("spilled_at"), "held_since": tr["held_since"]}
                        for tr in self.tracers],
            "catchers": [{k: c[k] for k in ("name", "node", "t_grab", "grab", "cooldown", "seed", "caught", "missed")}
                         for c in self.catchers],
            "lips": [{"index": i, "src": ln["src"], "dst": ln["dst"], "lip": ln["lip"], "width": ln["width"],
                      "releases": ln["releases"]} for i, ln in enumerate(self.links) if ln["kind"] == "lip"],
        }

    def _link_summary(self, i: int, ln: dict) -> dict:
        fa, ff = ln["first_arrival_step"], ln["first_flow_step"]
        out = {"index": i, "kind": ln["kind"], "src": ln["src"], "dst": ln["dst"], "delay": ln["delay"],
               "first_flow": None if ff is None else ff * self.dt,
               "first_arrival": None if fa is None else fa * self.dt, "volume_moved": from_q(ln["flow_sum"]),
               "volume_moved_q": ln["flow_sum"]}
        if ln["kind"] == "spill":
            out["max_outer"] = ln["max_outer"]
        if ln["kind"] == "boil":
            out.update(last_x=ln["last_x"], last_y=ln["last_y"], last_abv=ln["last_abv"], cuts=ln["cuts"])
        if "parts_sum" in ln:
            out["parts_moved_q"] = dict(ln["parts_sum"])
        if ln["kind"] == "tip":
            out.update(tips=[dict(t) for t in ln["tips"]], t_dump=ln["t_dump"], dump_rate=ln["dump_rate"])
        if ln["kind"] == "siphon":
            out.update(primes=list(ln["primes"]), breaks=[dict(b) for b in ln["breaks"]], primed=ln["primed"])
        if ln["kind"] == "mariotte":
            out.update(constant_flow=from_q(qdiv(ln["const_out"], self.dt_q)), constant_time=ln["const_steps"] * self.dt,
                       falling_from=ln["falling_from"])
        if ln["kind"] == "pipe":
            h = from_q(ln["max_head"])
            v = ln["radius"] ** 2 * self.g * h / (8.0 * self.fluid.nu * ln["length"])
            re = reynolds(self.fluid.nu, v, 2.0 * ln["radius"])
            out.update(tau=ln["tau"], reynolds_max=re, laminar=re < RE_LAMINAR)
        if ln["kind"] == "orifice":
            h = from_q(ln["max_head"])
            if ln["viscous"]:
                cd, re = viscous_orifice(h, ln["diameter"], self.fluid.nu, ln["cd"], ln["delta"], self.g)
            else:
                cd = ln["cd"]
                re = cd * math.sqrt(2.0 * self.g * h) * ln["diameter"] / self.fluid.nu
            out.update(reynolds_max=re, viscous=ln["viscous"], cd_at_max_head=cd, min_cd=ln["min_cd"],
                       viscous_regime=re < ORIFICE_RE_STEADY)
        return out

    def _channel_summary(self, ch: dict) -> dict:
        """樋の安定性（Courant 数）とあふれの有無。Courant = c·dt/Δx、c = dQ/dA（表の傾き）。"""
        mine = [ln for ln in self.links if ln.get("channel") == ch["prefix"]]
        v_max = max((from_q(ln["max_v"]) for ln in mine if ln["kind"] == "rating"), default=0.0)
        spills = [ln for ln in mine if ln["kind"] == "spill"]
        celerity = 0.0
        for (v0, q0), (v1, q1) in zip(ch["table"], ch["table"][1:]):
            if v0 > v_max:
                break
            celerity = max(celerity, (q1 - q0) / (v1 - v0) * ch["cell_length"])
        lips = [ln for ln in mine if ln["kind"] == "lip"]
        return {"max_cell_volume": v_max, "overflow": v_max >= ch["volume_cap"], "celerity": celerity,
                "pond_max_depth": max((ch_depth(ln, self) for ln in lips), default=None),
                "courant": celerity * self.dt / ch["cell_length"], "depth_cap": ch["depth_cap"],
                "spilled": sum(from_q(ln["flow_sum"]) for ln in spills),
                "leaked": sum(from_q(ln["flow_sum"]) for ln in mine if ln.get("crack")),
                "rained_q": sum(ln["flow_sum"] for ln in mine if ln.get("rain")),
                "max_outer_depth": max((ln["max_outer"] for ln in spills), default=None)}


def ch_depth(ln: dict, net: "Network") -> float:
    """lip の手前のセルにたまった最大の水深。"""
    nd = net.nodes[ln["src"]]
    return nd["section"].depth(from_q(ln["max_v"]) / nd["length"])


# --- JSON の仕様から組む --------------------------------------------------

LINK_BUILDERS = ("pipe", "orifice", "weir", "cam", "pump", "valve", "float_valve", "tip", "siphon", "mariotte", "boil",
                 "evaporate")
SPEC_KEYS = {"dt", "fluid", "temp_c", "mu", "mu_scale", "species", "tanks", "channels", "flumes", "rains", "links",
             "tracers", "t_end", "trace", "sample_every", "outlets", "inflows", "rounding"}


def build(spec: dict) -> Network:
    """仕様（dict / JSON）から回路を組む。

    {"dt": 0.02, "fluid": "water",
     "tanks": [{"name": "sea", "area": 1, "volume": 1}],
     "channels": [{"prefix": "flume", "dst": "sea", "length": 2, "cells": 20, "slope": 0.05,
                   "section": {"shape": "semicircle", "size": 0.04}, "law": "manning"}],
     "links": [{"law": "pump", "src": "sea", "dst": "flume:0", "rate": 5e-5}],
     "tracers": [{"name": "somen", "node": "flume:0", "release": 30, "stop": "sea"}]}

    channels・flumes は tanks の後、links の前に作る（樋のセルを管の行き先に使えるように）。
    flumes は {"prefix", "dst", "segments": [...], "spill"}（Network.flume を参照）。
    rains は {"prefix", "intensity", "src", "start", "until"}（Network.rain。樋の全セルに一様な雨）。
    links の law が LINK_BUILDERS 以外なら、register_law で登録した法則として net.link に渡す。
    """
    extra = set(spec) - SPEC_KEYS
    if extra:
        raise ValueError(f"unknown spec keys: {sorted(extra)}")
    species = spec.get("species")
    volume_of = None
    if species and set(species) == {"ethanol", "water"}:  # 🍷 エタノール・水なら密度表で体積を出す
        from fluids import mix_volume
        volume_of = lambda p: mix_volume(p["ethanol"], p["water"])
    net = Network(spec.get("dt", 0.1), spec.get("fluid", "water"), temp_c=spec.get("temp_c"), mu=spec.get("mu"),
                  mu_scale=spec.get("mu_scale"), species=species, volume_of=volume_of,
                  outlets=spec.get("outlets", "exact"), inflows=spec.get("inflows", "uniform"),
                  rounding=spec.get("rounding", "q48"))
    for t in spec.get("tanks", []):
        net.tank(**t)
    for c in spec.get("channels", []):
        net.channel(**c)
    for fl in spec.get("flumes", []):
        net.flume(**fl)
    for rn in spec.get("rains", []):  # 🌧️ {"prefix", "intensity", "src", "start", "until"}
        net.rain(**rn)
    for ln in spec.get("links", []):
        ln = dict(ln)
        law = ln.pop("law")
        if law in LINK_BUILDERS:
            getattr(net, law)(**ln)
        else:
            net.link(law, ln.pop("src"), ln.pop("dst"), **ln)
    for tr in spec.get("tracers", []):
        net.tracer(**tr)
    return net


def simulate(spec: dict) -> dict:
    """build(spec) して spec["t_end"] まで流す。Node から JSON のまま呼べる:
    ``node py_launch.cjs --call simulate -- network.py '<spec の JSON>'``

    2^53 を超える整数（Q の体積）は JSON で桁落ちしないよう 10 進文字列にして返す。
    """
    if "t_end" not in spec:
        raise ValueError("spec needs t_end")
    net = build(spec)
    res = net.run(spec["t_end"], sample_every=spec.get("sample_every"), trace=spec.get("trace", ()))
    return jsonable(res)


def jsonable(obj):
    """JSON 越しでも値が変わらない形へ（大きな整数 → 文字列、inf/nan → 文字列）。"""
    if isinstance(obj, bool) or obj is None or isinstance(obj, str):
        return obj
    if isinstance(obj, int):
        return obj if abs(obj) <= 2 ** 53 else str(obj)
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else ("inf" if obj > 0 else "-inf" if obj < 0 else "nan")
    if isinstance(obj, dict):
        return {k: jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [jsonable(v) for v in obj]
    return obj


# --- 解析の物差し ---------------------------------------------------------

def critical_gain(capacity: float, leak: float, tau_p: float, stages: int) -> float:
    """浮き弁ループ K / ((capacity·s + leak)(1 + tau_p·s)^n) の限界ゲイン（位相 −π で |L| = 1）。

    n ≤ 1 は位相が −π に届かず inf（どれだけ強くしても振動が育たない）。
    """
    if stages <= 1:
        return math.inf

    def phase(w: float) -> float:
        return math.atan2(capacity * w, leak) + stages * math.atan(tau_p * w)

    lo, hi = 0.0, 1.0
    while phase(hi) < math.pi:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if phase(mid) < math.pi else (lo, mid)
    w = 0.5 * (lo + hi)
    return math.hypot(capacity * w, leak) * (1 + (tau_p * w) ** 2) ** (stages / 2)


def inverse_table(points, y: float) -> float:
    """単調な区分線形表の逆引き（流量 → 水位）。"""
    pts = sorted((float(a), float(b)) for a, b in points)
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        if min(y0, y1) <= y <= max(y0, y1) and y0 != y1:
            return x0 + (x1 - x0) * (y - y0) / (y1 - y0)
    raise ValueError(f"{y} is outside the table range")


def window(values: list[float], frac_start: float, frac_end: float = 1.0) -> list[float]:
    n = len(values)
    a = int(n * frac_start)
    return values[a: max(a + 1, int(n * frac_end))]


def peak_to_peak(values: list[float]) -> float:
    return max(values) - min(values) if values else 0.0


def settling_time(values: list[float], target: float, dt: float, band: float = 0.02) -> float:
    """|x − target| ≤ band·|target| に入ったまま出なくなる時刻。最後まで入らなければ inf。"""
    tol = band * abs(target)
    for i in range(len(values) - 1, -1, -1):
        if abs(values[i] - target) > tol:
            return math.inf if i == len(values) - 1 else (i + 1) * dt
    return 0.0
