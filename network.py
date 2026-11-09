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
- どの管にも ``delay``（輸送の遅れ）を付けられる。管の中の水も総量に数える

運ばれる物 ``tracer``（そうめん・葉っぱ・浮き）: 入った節点に先にいた水が全部出たら、いちばん多く
流れた管で次へ進む（押し出し流れ・先入れ先出し）。水槽では平均滞留時間 V/Q、樋では L/v に当たる。
太さ ``size`` を持たせると、継ぎ目の縁の越流水深がそれより浅いあいだ引っかかる（``PASSABLE``）。

1 ステップでは全部の管の流量を同じ状態から求め（同時更新）、出す側の水が足りなければ整数のまま按分する。
移動は引いた整数をそのまま足すだけなので、総量は 1 ulp も変わらない。``build(spec)`` で JSON から組める。
"""

from __future__ import annotations

import math
from typing import Callable

from fixedq import Cam, DelayLine, clamp, from_q, lag_gain, lag_step, qdiv, qmul, qsqrt, to_q
from fluids import G, RE_LAMINAR, Section, channel_flow, fall_time, get_fluid, lip_flow, reynolds, side_spill
from fluids import superelevation

LAWS: dict[str, Callable] = {}
PASSABLE: dict[str, Callable] = {}  # 運ばれる物がその管を通れるか（net, link, tracer）→ bool。無ければ常に通れる


def register_law(kind: str):
    """管の法則を登録する。関数は (net, link, n) を受け、このステップに送りたい量（Q）を返す。"""

    def deco(fn):
        LAWS[kind] = fn
        return fn

    return deco


@register_law("pipe")
def _law_pipe(net, ln, n):
    ln["max_head"] = max(ln["max_head"], net.head_q(ln["src"]))
    return qmul(net.nodes[ln["src"]]["v"], ln["k"])


@register_law("orifice")
def _law_orifice(net, ln, n):
    h = net.head_q(ln["src"])
    return qmul(ln["coef"], qsqrt(h)) if h > 0 else 0


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
    """継ぎ目の段（次の竹の縁が lip だけ高い）。セルの水深が縁を越えた分だけ流れる。"""
    v = net.nodes[ln["src"]]["v"]
    ln["max_v"] = max(ln["max_v"], v)
    h = from_q(net.head_q(ln["src"]))
    return qmul(to_q(lip_flow(ln["section"], h, ln["lip"], ln["cd"], net.g)), net.dt_q)


def _lip_passable(net, ln, tr):
    """越流水深（水深 − 縁）が運ばれる物の太さ以上なら越えられる。"""
    return from_q(net.head_q(ln["src"])) - ln["lip"] >= tr["size"]


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
    if a > 0 and ln["radius"] is not None:
        v = from_q(qdiv(net.last_flow[ln["main"]], net.dt_q)) / a  # そのセルの本流（前のステップ）の流速
        h_out = h + superelevation(v, sec.top_width(h), ln["radius"], net.g) / 2.0
        ln["max_outer"] = max(ln["max_outer"], h_out)
        amount = qmul(to_q(side_spill(nd["length"], h_out, ln["rim"], ln["cd"], net.g)), net.dt_q)
    return max(amount, nd["v"] - ln["v_cap"])


@register_law("float_valve")
def _law_float_valve(net, ln, n):
    u = net.head_q(ln["sensor"])
    for i in range(ln["stages"]):
        ln["z"][i] = lag_step(ln["z"][i], u, ln["k"])
        u = ln["z"][i]
    q = clamp(qmul(ln["gain"], ln["setpoint"] - u), 0, ln["qmax"])
    return qmul(q, net.dt_q)


class Network:
    def __init__(self, dt: float, fluid: str = "water", g: float = G):
        if dt <= 0:
            raise ValueError("dt must be positive")
        self.dt = dt
        self.dt_q = to_q(dt)
        self.fluid = get_fluid(fluid)
        self.g = g
        self.nodes: dict[str, dict] = {}
        self.links: list[dict] = []
        self.channels: dict[str, dict] = {}
        self.tracers: list[dict] = []
        self._groups: tuple[int, dict] | None = None
        self.last_flow: list[int] = []  # 前のステップに各管で動いた量（Q）

    # --- 節点 -----------------------------------------------------------
    def _node(self, name: str, **data) -> dict:
        if name in self.nodes:
            raise ValueError(f"duplicate node {name!r}")
        self.nodes[name] = data
        return data

    def tank(self, name: str, area: float, head: float = 0.0, volume: float | None = None) -> "Network":
        if area <= 0:
            raise ValueError("area must be positive")
        v = head * area if volume is None else volume
        self._node(name, kind="tank", area=area, area_q=to_q(area), v=to_q(v))
        return self

    def channel(self, prefix: str, dst: str, *, length: float, cells: int, slope: float, section,
                law: str = "manning", n_manning: float = 0.010, beta: float | None = None,
                max_depth: float | None = None, table_points: int = 256, delay: float = 0.0,
                lip: float = 0.0, bend_radius: float | None = None, spill: str | None = None,
                cd: float = 0.62) -> dict:
        """傾いた樋。区間セル <prefix>:0 … <prefix>:cells-1 を作り、末端を dst へつなぐ。

        各セルの流出は「セル体積 → 流量」の表（rating カム）で引く。表は等流の式から一度だけ作り、
        水深の浅い側ほど点を密にする（(i/N)² 刻み）。表の上端（半円なら縁、長方形なら max_depth）を
        超えた分は流量が頭打ちになるので、結果の overflow で知らせる。

        - lip > 0: 末端は「次の竹の縁が lip だけ高い」継ぎ目。末端セルに水がたまり、越えた分だけ流れる
        - bend_radius: 曲がり（半径 R）。外側の水面が縁を越えた分は各セルから spill の節点へあふれる
        - delay: 末端の管の遅れ（段差の自由落下などに使う）
        """
        if cells < 1 or length <= 0 or slope <= 0:
            raise ValueError("channel needs cells >= 1, length > 0, slope > 0")
        if dst not in self.nodes:
            raise ValueError(f"unknown node {dst!r}")
        if bend_radius is not None and spill is None:
            raise ValueError("a bend needs a spill node")
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
                self._link("lip", name, nxt, delay, section=sec, lip=lip, cd=cd, max_v=0, channel=prefix)
            else:
                self._link("rating", name, nxt, delay if nxt == dst else 0.0, cam=cam, max_v=0, channel=prefix)
            if bend_radius is not None or (spill is not None and nxt == dst and lip > 0):
                # 曲がりの全セルと、縁の手前でたまるセルには縁からあふれる口を付ける
                self._link("spill", name, spill, 0.0, main=len(self.links) - 1, radius=bend_radius, rim=h_cap,
                           v_cap=to_q(a_max * lc), cd=cd, max_outer=0.0, channel=prefix)
        info = {"prefix": prefix, "cells": names, "length": length, "cell_length": lc, "slope": slope,
                "section": sec.to_dict(), "law": law, "n_manning": n_manning, "beta": beta, "depth_cap": h_cap,
                "volume_cap": a_max * lc, "table": table, "lip": lip, "bend_radius": bend_radius}
        self.channels[prefix] = info
        return info

    def flume(self, prefix: str, dst: str, segments: list[dict], spill: str | None = None) -> list[dict]:
        """竹をつないだ流しそうめんの台。segments を上流から順に並べ、<prefix>/<i>:<cell> を作る。

        各 segment は channel() の引数（length・cells・slope・section・law …）に加えて
        ``joint``（次の竹との継ぎ目）を持てる: ``{"drop": 段差 m}``（自由落下の遅れ）か
        ``{"lip": 縁の高さ m}``（たまってから越える）。``bend_radius`` の竹は spill へあふれうる。
        """
        if not segments:
            raise ValueError("flume needs at least one segment")
        infos = [None] * len(segments)
        nxt = dst
        for i in range(len(segments) - 1, -1, -1):  # 下流から作る（行き先の節点を先に用意する）
            seg = dict(segments[i])
            joint = seg.pop("joint", None) or {}
            unknown = set(joint) - {"drop", "lip"}
            if unknown:
                raise ValueError(f"unknown joint keys: {sorted(unknown)}")
            infos[i] = self.channel(f"{prefix}/{i}", nxt, delay=fall_time(joint.get("drop", 0.0), self.g),
                                    lip=joint.get("lip", 0.0), spill=spill, **seg)
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
        link = {
            "kind": kind, "src": src, "dst": dst, "delay": steps * self.dt, "delay_steps": steps,
            "line": DelayLine(steps - 1) if steps else None,
            "pending": 0, "flow_sum": 0, "max_head": 0, "first_arrival_step": None, **extra,
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
        tau = self.nodes[src]["area"] / c
        return self._link("pipe", src, dst, delay, radius=radius, length=length,
                          conductance=c, tau=tau, k=lag_gain(self.dt, tau))

    def orifice(self, src: str, dst: str, *, area: float, cd: float = 0.61, delay: float = 0.0) -> dict:
        coef = to_q(cd * area * math.sqrt(2.0 * self.g) * self.dt)
        return self._link("orifice", src, dst, delay, coef=coef, orifice_area=area, cd=cd)

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

    def pump(self, src: str, dst: str, *, rate: float, start: float = 0.0, until: float | None = None,
             delay: float = 0.0) -> dict:
        """一定流量。start から until まで動く。"""
        return self._link("pump", src, dst, delay, amount=qmul(to_q(rate), self.dt_q), rate=rate,
                          start_step=int(round(start / self.dt)),
                          until_step=None if until is None else int(round(until / self.dt)))

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

    # --- 運ばれる物 -------------------------------------------------------
    def tracer(self, name: str, node: str, release: float = 0.0, stop=(), size: float = 0.0) -> dict:
        """運ばれる物を node に release 時刻で放す。stop（節点名か、その並び）に入るか、出口のない節点で止まる。

        size は太さ [m]。継ぎ目の段（lip）では越流水深がこれより浅いと引っかかって待つ。
        """
        stops = (stop,) if isinstance(stop, str) else tuple(stop)
        for nd in (node, *stops):
            if nd not in self.nodes:
                raise ValueError(f"unknown node {nd!r}")
        tr = {"name": name, "node": node, "release_step": int(round(release / self.dt)), "state": "waiting",
              "credit": 0, "pipe_left": 0, "pipe_dst": None, "pipe_t": None, "pipe_f": 0.0, "path": [],
              "arrival": None, "stop": stops, "size": size, "stuck_steps": 0, "stuck_at": None,
              "held": False}
        self.tracers.append(tr)
        return tr

    # --- 状態 -----------------------------------------------------------
    def head_q(self, name: str) -> int:
        nd = self.nodes[name]
        if nd["kind"] == "tank":
            return qdiv(nd["v"], nd["area_q"])
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

    # --- 1 ステップ -----------------------------------------------------
    def step(self, n: int) -> list[int]:
        if len(self.last_flow) != len(self.links):
            self.last_flow = [0] * len(self.links)
        req = [max(LAWS[ln["kind"]](self, ln, n), 0) for ln in self.links]
        # 出す側の水が足りなければ、その節点から出る管で整数のまま按分する
        for src, idx in self._by_src().items():
            want = sum(req[i] for i in idx)
            have = self.nodes[src]["v"]
            if want > have:
                for i in idx:
                    req[i] = req[i] * have // want
                # 床で余った端数は最大の管へ（足りない水槽は空になるまで出し切る）
                req[max(idx, key=lambda i: req[i])] += have - sum(req[i] for i in idx)
        for i, ln in enumerate(self.links):
            self.nodes[ln["src"]]["v"] -= req[i]
            ln["flow_sum"] += req[i]
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

    def _move_tracers(self, n: int, req: list[int]) -> None:
        """先入れ先出し。先にいた水（credit）を流出で減らし、0 を切ったステップ内の位置 f まで時刻を補間する。"""
        by_src = self._by_src()
        for tr in self.tracers:
            if tr["state"] == "waiting":
                if n < tr["release_step"]:
                    continue
                # 放した時点でその節点にいた水（このステップの流入を含む）が「先にいた水」
                tr["state"] = "node"
                tr["credit"] = self.nodes[tr["node"]]["v"] + sum(req[i] for i in by_src.get(tr["node"], []))
                tr["path"].append({"node": tr["node"], "t": n * self.dt})
            if tr["state"] == "pipe":
                tr["pipe_left"] -= 1
                if tr["pipe_left"] <= 0:
                    self._enter(tr, tr["pipe_dst"], tr["pipe_t"], tr["pipe_f"], req)
                continue
            if tr["state"] != "node":
                continue
            outs = by_src.get(tr["node"], [])
            if not outs:
                tr["state"], tr["arrival"] = "done", tr["path"][-1]["t"]
                continue
            out_now = sum(req[i] for i in outs)
            tr["credit"] -= out_now
            if tr["credit"] <= 0 and out_now > 0:
                ln = self.links[max(outs, key=lambda i: req[i])]
                ok = PASSABLE.get(ln["kind"])
                if ok is not None and not ok(self, ln, tr):
                    tr["stuck_steps"] += 1  # 引っかかり: 先にいた水は出たが、ここは越えられない
                    tr["stuck_at"], tr["held"] = tr["node"], True
                    continue
                # このステップのうち、越えた後に残る割合（引っかかりから外れた直後は 0）
                f = 0.0 if tr["held"] else -tr["credit"] / out_now
                tr["held"] = False
                t_cross = (n + 1 - f) * self.dt
                if ln["delay_steps"]:
                    tr.update(state="pipe", pipe_left=ln["delay_steps"], pipe_dst=ln["dst"],
                              pipe_t=t_cross + ln["delay"], pipe_f=f)
                else:
                    self._enter(tr, ln["dst"], t_cross, f, req)

    def _enter(self, tr: dict, node: str, t: float, f: float, req: list[int]) -> None:
        tr["node"], tr["state"] = node, "node"
        outs = self._by_src().get(node, [])
        # 入った後のステップの残り f の間に、すでに先へ出ていった分を差し引く
        tr["credit"] = self.nodes[node]["v"] - int(f * sum(req[i] for i in outs))
        tr["path"].append({"node": node, "t": t})
        if node in tr["stop"] or not outs:
            tr["state"], tr["arrival"] = "done", t

    def run(self, t_end: float, *, sample_every: int | None = None, trace=()) -> dict:
        """t_end まで流す。series はプロッター記録（抽出）、trace は指定した節点の毎ステップ水位。"""
        steps = int(round(t_end / self.dt))
        if sample_every is None:
            sample_every = max(1, steps // 200)
        total0 = self.total_q()
        worst = 0
        series = []
        traces = {name: [] for name in trace}
        tank_names = [k for k, nd in self.nodes.items() if nd["kind"] == "tank"]
        for n in range(steps):
            last = self.step(n)
            worst = max(worst, abs(self.total_q() - total0))
            for name in trace:
                traces[name].append(self.head(name))
            if (n + 1) % sample_every == 0 or n + 1 == steps:
                row = {"t": round((n + 1) * self.dt, 9)}
                row.update({f"h:{k}": self.head(k) for k in tank_names})
                for prefix, ch in self.channels.items():
                    row[f"h:{prefix}"] = [self.head(c) for c in ch["cells"]]
                row.update({f"q:{i}": from_q(qdiv(q, self.dt_q)) for i, q in enumerate(last)})
                series.append(row)
        return {
            "fluid": self.fluid.key, "dt": self.dt, "steps": steps, "series": series, "trace": traces,
            "final": series[-1] if series else {}, "total_q": total0, "conservation_error_q": worst,
            "volumes_q": {k: nd["v"] for k, nd in self.nodes.items()},
            "links": [self._link_summary(i, ln) for i, ln in enumerate(self.links)],
            "channels": {p: self._channel_summary(ch) for p, ch in self.channels.items()},
            "tracers": [{"name": tr["name"], "state": tr["state"], "arrival": tr["arrival"], "path": tr["path"],
                         "size": tr["size"], "stuck_time": tr["stuck_steps"] * self.dt, "stuck_at": tr["stuck_at"]}
                        for tr in self.tracers],
        }

    def _link_summary(self, i: int, ln: dict) -> dict:
        fa = ln["first_arrival_step"]
        out = {"index": i, "kind": ln["kind"], "src": ln["src"], "dst": ln["dst"], "delay": ln["delay"],
               "first_arrival": None if fa is None else fa * self.dt, "volume_moved": from_q(ln["flow_sum"])}
        if ln["kind"] == "pipe":
            h = from_q(ln["max_head"])
            v = ln["radius"] ** 2 * self.g * h / (8.0 * self.fluid.nu * ln["length"])
            re = reynolds(self.fluid.nu, v, 2.0 * ln["radius"])
            out.update(tau=ln["tau"], reynolds_max=re, laminar=re < RE_LAMINAR)
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
                "max_outer_depth": max((ln["max_outer"] for ln in spills), default=None)}


def ch_depth(ln: dict, net: "Network") -> float:
    """lip の手前のセルにたまった最大の水深。"""
    nd = net.nodes[ln["src"]]
    return nd["section"].depth(from_q(ln["max_v"]) / nd["length"])


# --- JSON の仕様から組む --------------------------------------------------

LINK_BUILDERS = ("pipe", "orifice", "weir", "cam", "pump", "float_valve")
SPEC_KEYS = {"dt", "fluid", "tanks", "channels", "flumes", "links", "tracers", "t_end", "trace", "sample_every"}


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
    links の law が LINK_BUILDERS 以外なら、register_law で登録した法則として net.link に渡す。
    """
    extra = set(spec) - SPEC_KEYS
    if extra:
        raise ValueError(f"unknown spec keys: {sorted(extra)}")
    net = Network(spec.get("dt", 0.1), spec.get("fluid", "water"))
    for t in spec.get("tanks", []):
        net.tank(**t)
    for c in spec.get("channels", []):
        net.channel(**c)
    for fl in spec.get("flumes", []):
        net.flume(**fl)
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
