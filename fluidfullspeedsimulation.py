"""🌊 FluidFullSpeedSimulation — 水・エタノールなどの作動流体で「流れ・遅れ・検証」を整数演算で回す。

MONIAC（フィリップス・マシン）の物理部品（水槽＝積分器、管・小孔・堰・カム、浮き弁、遅延管、2 台連結）を
流体ごとに再現し、閉じた式と突き合わせる。経済の解釈は扱わない。

名前で呼べる（``__main__`` を通らない）:

- Python: ``from fluidfullspeedsimulation import run_experiment; run_experiment("E02")``
- Node:   ``node py_launch.cjs --call run_experiment -- fluidfullspeedsimulation.py E02``
          （親リポジトリの py_launch.cjs。kwargs は ``call(script, "run_experiment", {args, kwargs})``）
- CLI:    ``python fluidfullspeedsimulation.py run E02``

終了コード: 0 = ✅ 全検証が通った / 1 = ❌ どれかが外れた / 2 = 引数の誤り（argparse）。
"""

from __future__ import annotations

import argparse
import json
import math
import random
import sys
from pathlib import Path

from fixedq import ONE, from_q, lag_step, qmul, shift_lag_step, to_q
from fluids import FLUIDS, REAL_FLUIDS, RE_LAMINAR, Section, channel_reynolds, fall_time, fluid_table, get_fluid, lip_head
from fluids import normal_depth
from fluids import poiseuille_head, poiseuille_tau, torricelli_empty_time, torricelli_head, transport_delay
from network import Network, critical_gain, inverse_table, jsonable, peak_to_peak, settling_time, window
import pixelwater as pw

FLUID_KEYS = REAL_FLUIDS  # 実験で回す流体（仮置きの流体は入れない）
ALL_FLUID_KEYS = tuple(FLUIDS)  # 描画で選べる流体（仮置きも含む）
SEA = 1.0  # 外の水（海）。総量を閉じるための大きな水槽 [m³]


# --- 検証の記録 -------------------------------------------------------------

def rel_err(value: float, expected: float) -> float:
    return abs(value - expected) / abs(expected) if expected else abs(value - expected)


def check(cid: str, icon: str, ok: bool, *, value=None, expected=None, tol=None, note: str = "") -> dict:
    delta = None
    if isinstance(value, (int, float)) and isinstance(expected, (int, float)):
        delta = value - expected
    return {"id": cid, "icon": icon, "ok": bool(ok), "value": jsonable(value), "expected": jsonable(expected),
            "delta": jsonable(delta), "tol": jsonable(tol), "note": note}


def close(cid: str, icon: str, value: float, expected: float, tol: float, note: str = "") -> dict:
    return check(cid, icon, rel_err(value, expected) <= tol, value=value, expected=expected, tol=tol, note=note)


def conserved(cid: str, res: dict) -> dict:
    return check(cid, "⚖️", res["conservation_error_q"] == 0, value=res["conservation_error_q"], expected=0,
                 tol=0, note="Σ水槽 + 管の中（Q 整数、ulp）")


def laminar(cid: str, res: dict) -> dict:
    pipes = [ln for ln in res["links"] if ln["kind"] == "pipe"]
    worst = max((ln["reynolds_max"] for ln in pipes), default=0.0)
    return check(cid, "〰️", worst < RE_LAMINAR, value=worst, expected=RE_LAMINAR, note="Re < 2300（層流の前提）")


def _sea(net: Network) -> Network:
    return net.tank("sea", 1.0, volume=SEA)


# --- E01 水槽は積分器 -------------------------------------------------------

def e01_tank_integrator(dt: float = 0.5, t_end: float = 100.0, rate: float = 1e-4, area: float = 0.01,
                        seed: int = 0) -> dict:
    fill = _sea(Network(dt)).tank("tank", area)
    fill.pump("sea", "tank", rate=rate)
    r1 = fill.run(t_end)
    per_step = qmul(to_q(rate), to_q(dt))

    tau_dt = 1.0
    drain = _sea(Network(tau_dt)).tank("tank", 0.01, head=0.30)
    pipe = drain.pipe("tank", "sea", radius=0.75e-3, length=0.50)
    r2 = drain.run(2.0 * pipe["tau"], trace=["tank"])
    hs = r2["trace"]["tank"]
    worst = max(abs(h - poiseuille_head(0.30, (n + 1) * tau_dt, pipe["tau"])) / 0.30 for n, h in enumerate(hs))

    rng = random.Random(seed)
    trials = [(rng.randrange(-ONE * 100, ONE * 100), rng.randrange(-ONE * 100, ONE * 100), rng.randint(0, 31))
              for _ in range(2000)]
    same = all(lag_step(y, x, ONE >> s) == shift_lag_step(y, x, s) for y, x, s in trials)

    checks = [
        check("integrate-exact", "🧮", r1["volumes_q"]["tank"] == r1["steps"] * per_step,
              value=r1["volumes_q"]["tank"], expected=r1["steps"] * per_step, tol=0, note="一定流入 → 体積 = n·(Q·dt) を整数で厳密に"),
        check("integrate-float", "📏", abs(from_q(r1["volumes_q"]["tank"]) - rate * t_end) <= (r1["steps"] + 1) / ONE,
              value=from_q(r1["volumes_q"]["tank"]), expected=rate * t_end, tol=(r1["steps"] + 1) / ONE,
              note="V = Q·t（Q·dt の丸め 1 ulp × n ステップ以内）"),
        check("drain-exp", "📉", worst < 1e-6, value=worst, expected=0.0, tol=1e-6, note="h = h0·e^(−t/τ) との最大ずれ / h0"),
        check("shift-lag", "🔀", same, value=len(trials), note="k = 2^−s の遅れ = ビットシフトだけ（2000 例で一致）"),
        conserved("conserve-fill", r1),
        conserved("conserve-drain", r2),
    ]
    return {"params": {"dt": dt, "t_end": t_end, "rate": rate, "area": area, "seed": seed},
            "summary": {"tau_drain": pipe["tau"], "drain_max_rel": worst},
            "runs": {"fill": r1["series"], "drain": r2["series"]}, "checks": checks}


# --- E02 流体ごとの排水 -----------------------------------------------------

def e02_fluid_drain(h0: float = 0.30, area: float = 0.01, radius: float = 0.75e-3, length: float = 0.50,
                    orifice_d: float = 4e-3, cd: float = 0.61, dt_pipe: float = 2.0, dt_orifice: float = 0.1) -> dict:
    a_orf = math.pi * (orifice_d / 2) ** 2
    t_empty = torricelli_empty_time(h0, area, a_orf, cd)
    checks, runs, summary, tor_traces, halves = [], {}, {}, {}, {}
    for key in FLUID_KEYS:
        fl = get_fluid(key)
        net = _sea(Network(dt_pipe, key)).tank("tank", area, head=h0)
        pipe = net.pipe("tank", "sea", radius=radius, length=length)
        tau = pipe["tau"]
        rp = net.run(1.5 * tau, trace=["tank"])
        hs = rp["trace"]["tank"]
        worst = max(abs(h - poiseuille_head(h0, (n + 1) * dt_pipe, tau)) / h0 for n, h in enumerate(hs))
        half = next((n + 1) * dt_pipe for n, h in enumerate(hs) if h <= h0 / 2)
        halves[key] = half

        net_o = _sea(Network(dt_orifice, key)).tank("tank", area, head=h0)
        net_o.orifice("tank", "sea", area=a_orf, cd=cd)
        ro = net_o.run(0.5 * t_empty, trace=["tank"])
        tor_traces[key] = ro["trace"]["tank"]
        h_mid = ro["trace"]["tank"][-1]
        h_mid_exact = torricelli_head(h0, ro["steps"] * dt_orifice, area, a_orf, cd)

        checks += [
            check(f"{key}:poiseuille", fl.icon + "📉", worst < 1e-6, value=worst, expected=0.0, tol=1e-6,
                  note="細管: h0·e^(−t/τ) との最大ずれ / h0"),
            close(f"{key}:torricelli", fl.icon + "🕳️", h_mid, h_mid_exact, 2e-3, note="小孔: t = T/2 の水位（前進 Euler）"),
            laminar(f"{key}:laminar", rp),
            conserved(f"{key}:conserve", rp),
            conserved(f"{key}:conserve-orifice", ro),
        ]
        runs[f"{key}:pipe"] = rp["series"]
        runs[f"{key}:orifice"] = ro["series"]
        summary[key] = {"nu": fl.nu, "tau": tau, "t_half": half, "reynolds": rp["links"][0]["reynolds_max"],
                        "torricelli_mid": h_mid, "torricelli_mid_exact": h_mid_exact}
    w = get_fluid("water")
    for key in FLUID_KEYS:
        if key == "water":
            continue
        fl = get_fluid(key)
        checks.append(close(f"{key}:speed", fl.icon + "⏱️", halves[key] / halves["water"], fl.nu / w.nu, 2e-3,
                            note="半減時間の比 = 動粘度の比 ν/ν_water"))
        checks.append(check(f"{key}:torricelli-same", fl.icon + "🟰", tor_traces[key] == tor_traces["water"],
                            note="小孔は ρ・μ によらない → 水とビット単位で同じ水位列"))
    summary["t_empty_torricelli"] = t_empty
    return {"params": {"h0": h0, "area": area, "radius": radius, "length": length, "orifice_d": orifice_d, "cd": cd,
                       "dt_pipe": dt_pipe, "dt_orifice": dt_orifice},
            "summary": summary, "runs": runs, "checks": checks}


# --- E03 輸送の遅れ ---------------------------------------------------------

def e03_transport_delay(h0: float = 0.30, area: float = 0.01, radius: float = 0.75e-3, length: float = 0.50,
                        dt: float = 0.01, t_end: float = 5.0) -> dict:
    checks, runs, summary, measured = [], {}, {}, {}
    for key in FLUID_KEYS:
        fl = get_fluid(key)
        delay = transport_delay(fl.nu, h0, radius, length)
        net = Network(dt, key).tank("A", area, head=h0).tank("B", area)
        net.pipe("A", "B", radius=radius, length=length, delay=delay)
        res = net.run(t_end)
        link = res["links"][0]
        measured[key] = link["first_arrival"]
        moved = link["volume_moved"]
        arrived = from_q(res["volumes_q"]["B"])
        checks += [
            check(f"{key}:arrival", fl.icon + "⏳", abs(link["first_arrival"] - delay) <= dt / 2 + 1e-12,
                  value=link["first_arrival"], expected=delay, tol=dt / 2, note="最初に届く時刻 = L / v"),
            check(f"{key}:in-pipe", fl.icon + "🚰", moved > arrived > 0, value=moved - arrived,
                  note="出た量 − 届いた量 = 管の中の水（> 0）"),
            conserved(f"{key}:conserve", res),
        ]
        runs[key] = res["series"]
        summary[key] = {"delay": delay, "first_arrival": link["first_arrival"], "in_pipe": moved - arrived}
    w = get_fluid("water")
    for key in FLUID_KEYS:
        if key != "water":
            fl = get_fluid(key)
            checks.append(close(f"{key}:delay-ratio", fl.icon + "⏱️", measured[key] / measured["water"],
                                fl.nu / w.nu, 2 * dt / measured["water"], note="遅れの比 = ν/ν_water"))
    return {"params": {"h0": h0, "area": area, "radius": radius, "length": length, "dt": dt, "t_end": t_end},
            "summary": summary, "runs": runs, "checks": checks}


# --- E04 堰とカム（非線形の出口） -------------------------------------------

CAM_TABLE = [[0.0, 0.0], [0.02, 1e-4], [0.05, 4e-4]]


def e04_weir(inflow: float = 2e-4, area: float = 0.05, width: float = 0.02, dt: float = 0.01,
             t_end: float = 150.0) -> dict:
    checks, runs, summary = [], {}, {}
    finals = {}
    for shape in ("rect", "v", "cam"):
        for key in ("water", "ethanol"):
            net = _sea(Network(dt, key)).tank("tank", area)
            net.pump("sea", "tank", rate=inflow)
            ln = net.cam("tank", "sea", table=CAM_TABLE) if shape == "cam" else net.weir(
                "tank", "sea", shape=shape, width=width)
            res = net.run(t_end)
            finals[(shape, key)] = res["volumes_q"]["tank"]
            if key != "water":
                continue
            if shape == "cam":
                expected = inverse_table(CAM_TABLE, inflow)
            else:
                expected = (inflow / ln["k"]) ** (2.0 / 3.0 if shape == "rect" else 2.0 / 5.0)
            h = res["final"]["h:tank"]
            out = res["final"]["q:1"]
            checks += [
                close(f"{shape}:head", "🪜", h, expected, 1e-3, note="定常水位（流入 = 堰の流出）"),
                close(f"{shape}:balance", "🔁", out, inflow, 1e-3, note="定常で流出 = 流入"),
                conserved(f"{shape}:conserve", res),
            ]
            runs[shape] = res["series"]
            summary[shape] = {"head": h, "head_exact": expected}
    for shape in ("rect", "v", "cam"):
        checks.append(check(f"{shape}:fluid-free", "🟰", finals[(shape, "ethanol")] == finals[(shape, "water")],
                            note="堰・カムは粘度によらない → 🍶 と 💧 でビット単位で同じ"))
    return {"params": {"inflow": inflow, "area": area, "width": width, "dt": dt, "t_end": t_end, "cam": CAM_TABLE},
            "summary": summary, "runs": runs, "checks": checks}


# --- E05 浮き弁と遅れ -------------------------------------------------------

def e05_float_valve(area: float = 0.01, tau_out: float = 20.0, setpoint: float = 0.20, qmax: float = 5e-3,
                    tau_link: float = 2.0, gain_low: float = 2e-3, gain_high: float = 2e-2, dt: float = 0.02,
                    t_end: float = 400.0) -> dict:
    """浮き → ワイヤー → 弁のつながりに遅れ（stages 段、各 tau_link）があると、強い弁ほど揺れが育つ。"""
    leak = area / tau_out  # 出口の線形コンダクタンス [m²/s]
    kc = critical_gain(area, leak, tau_link, 3)
    checks, runs, summary = [], {}, {}
    p2p, err = {}, {}
    for stages in (1, 3):
        for name, gain in (("low", gain_low), ("high", gain_high)):
            net = _sea(Network(dt)).tank("tank", area)
            net.float_valve("sea", "tank", sensor="tank", setpoint=setpoint, gain=gain, qmax=qmax,
                            stages=stages, tau=tau_link)
            # 線形の出口（時定数 tau_out）を細管で作る: τ = A·8νL/(π r⁴ g) を r について解く
            nu = net.fluid.nu
            length = 0.5
            radius = (area * 8 * nu * length / (math.pi * net.g * tau_out)) ** 0.25
            net.pipe("tank", "sea", radius=radius, length=length)
            res = net.run(t_end, trace=["tank"])
            late = window(res["trace"]["tank"], 0.75)
            tag = f"n{stages}:{name}"
            p2p[tag] = peak_to_peak(late)
            err[tag] = abs(sum(late) / len(late) - setpoint)
            runs[tag] = res["series"]
            summary[tag] = {"gain": gain, "stages": stages, "late_p2p": p2p[tag], "late_mean": sum(late) / len(late),
                            "steady_exact": gain * setpoint / (gain + leak)}
            checks.append(conserved(f"{tag}:conserve", res))
            if not (stages == 3 and gain > kc):
                checks.append(close(f"{tag}:steady", "🎯", sum(late) / len(late), gain * setpoint / (gain + leak), 1e-3,
                                    note="定常水位 h* = K·h_set/(K + A/τ_out)"))
    still = 1e-4 * setpoint
    checks += [
        check("gain-order", "📐", gain_low < kc < gain_high, value=kc, note="限界ゲイン K_c（3 段の遅れ、位相 −π）"),
        check("n1:stronger-closer", "🎯", err["n1:high"] < err["n1:low"], value=err["n1:high"], expected=err["n1:low"],
              note="遅れ 1 段: 強い弁ほど設定値に近い"),
        check("n1:high-still", "🫧", p2p["n1:high"] < still, value=p2p["n1:high"], tol=still, note="遅れ 1 段: 揺れは消える"),
        check("n3:low-still", "🫧", p2p["n3:low"] < still, value=p2p["n3:low"], tol=still, note="遅れ 3 段・K < K_c: 揺れは消える"),
        check("n3:high-swing", "🌊", p2p["n3:high"] > 0.1 * setpoint, value=p2p["n3:high"], tol=0.1 * setpoint,
              note="遅れ 3 段・K > K_c: 揺れが残る（弁の 0..qmax で頭打ち）"),
    ]
    summary["critical_gain"] = kc
    return {"params": {"area": area, "tau_out": tau_out, "setpoint": setpoint, "qmax": qmax, "tau_link": tau_link,
                       "gain_low": gain_low, "gain_high": gain_high, "dt": dt, "t_end": t_end},
            "summary": summary, "runs": runs, "checks": checks}


# --- E06 2 台連結（遅延つきの環） --------------------------------------------

def e06_coupled_loop(area: float = 0.01, h0: float = 0.05, radius: float = 2e-3, length_ab: float = 0.5,
                     length_ba: float = 1.0, delay_ab: float = 5.0, delay_ba: float = 10.0, dt: float = 0.1,
                     t_end: float = 2000.0) -> dict:
    net = Network(dt).tank("A", area, head=h0).tank("B", area)
    ab = net.pipe("A", "B", radius=radius, length=length_ab, delay=delay_ab)
    ba = net.pipe("B", "A", radius=radius, length=length_ba, delay=delay_ba)
    res = net.run(t_end)
    total = h0 * area
    q = total / (ab["tau"] + ba["tau"] + ab["delay"] + ba["delay"])
    va, vb = from_q(res["volumes_q"]["A"]), from_q(res["volumes_q"]["B"])
    la, lb = res["links"]
    checks = [
        close("A:steady", "🅰️", va, q * ab["tau"], 1e-3, note="V_A* = q·τ_AB、q = V/(τ_AB + τ_BA + d_AB + d_BA)"),
        close("B:steady", "🅱️", vb, q * ba["tau"], 1e-3, note="V_B* = q·τ_BA"),
        close("in-pipes", "🚰", total - va - vb, q * (ab["delay"] + ba["delay"]), 1e-3, note="管の中 = q·(d_AB + d_BA)"),
        check("A→B:arrival", "⏳", abs(la["first_arrival"] - delay_ab) <= dt / 2 + 1e-12, value=la["first_arrival"],
              expected=delay_ab, tol=dt / 2, note="B に最初に届く時刻"),
        check("B→A:echo", "↩️", abs(lb["first_arrival"] - (delay_ab + delay_ba)) <= dt + 1e-12,
              value=lb["first_arrival"], expected=delay_ab + delay_ba, tol=dt, note="跳ね返りが A に戻る時刻 = 往復の遅れ（同時更新で +1 ステップまで）"),
        laminar("laminar", res),
        conserved("conserve", res),
    ]
    return {"params": {"area": area, "h0": h0, "radius": radius, "length_ab": length_ab, "length_ba": length_ba,
                       "delay_ab": delay_ab, "delay_ba": delay_ba, "dt": dt, "t_end": t_end},
            "summary": {"q": q, "tau_ab": ab["tau"], "tau_ba": ba["tau"], "VA": va, "VB": vb},
            "runs": {"loop": res["series"]}, "checks": checks}


# --- E07 作動流体で変わる回路の速さ -----------------------------------------

def e07_fluid_speed(inflow: float = 2e-6, area: float = 0.01, radius: float = 1.2e-3, length: float = 0.5,
                    dt: float = 0.5, horizon: float = 12.0) -> dict:
    """ポンプ → 水槽1 → 細管 → 水槽2 → 細管 → 海。流量は流体によらず、水位と速さは ν に比例する。"""
    tau_max = max(poiseuille_tau(get_fluid(k).nu, area, radius, length) for k in FLUID_KEYS)
    t_end = horizon * tau_max
    checks, runs, summary = [], {}, {}
    for key in FLUID_KEYS:
        fl = get_fluid(key)
        net = _sea(Network(dt, key)).tank("t1", area).tank("t2", area)
        net.pump("sea", "t1", rate=inflow)
        p1 = net.pipe("t1", "t2", radius=radius, length=length)
        net.pipe("t2", "sea", radius=radius, length=length)
        res = net.run(t_end, trace=["t2"])
        h2_exact = inflow * p1["tau"] / area
        h2 = res["final"]["h:t2"]
        summary[key] = {"nu": fl.nu, "tau": p1["tau"], "h2": h2, "h2_exact": h2_exact,
                        "settle": settling_time(res["trace"]["t2"], h2_exact, dt)}
        checks += [
            close(f"{key}:outflow", fl.icon + "🔁", res["final"]["q:2"], inflow, 1e-3, note="定常の流出 = ポンプ流量"),
            close(f"{key}:head", fl.icon + "📏", h2, h2_exact, 1e-3,
                  note="定常水位 h* = Q·τ/A（δ ≈ dt/2τ は流入をステップ頭に入れる離散化の分）"),
            laminar(f"{key}:laminar", res),
            conserved(f"{key}:conserve", res),
        ]
        runs[key] = res["series"]
    w = get_fluid("water")
    for key in FLUID_KEYS:
        if key != "water":
            fl = get_fluid(key)
            ratio = summary[key]["settle"] / summary["water"]["settle"]
            checks.append(close(f"{key}:settle-ratio", fl.icon + "🏎️", ratio, fl.nu / w.nu, 1e-2,
                                note="2% 整定時間の比 = ν/ν_water"))
    return {"params": {"inflow": inflow, "area": area, "radius": radius, "length": length, "dt": dt,
                       "horizon": horizon, "t_end": t_end},
            "summary": summary, "runs": runs, "checks": checks}


# --- E08 傾いた樋（流しそうめんの台の下ごしらえ） ----------------------------

FLUME_CASES = {
    # 名前: (断面, 法則, 流量 [m³/s], 表の上限水深, dt, 流体)
    "rect-laminar": ({"shape": "rect", "size": 0.05}, "laminar", 2e-6, 0.005, 0.05, FLUID_KEYS),
    "bamboo-manning": ({"shape": "semicircle", "size": 0.04}, "manning", 5e-5, None, 0.01, ("water", "ethanol")),
}


def e08_flume(length: float = 2.0, cells: int = 20, slope: float = 0.05, n_manning: float = 0.010) -> dict:
    """ポンプ → 樋（cells 区間）→ 海。等流水深に落ち着いたところでトレーサー（そうめん役）を放し、
    樋を抜ける時間 L/v と比べる。層流の薄い流れは ν に効き、Manning（乱流）は流体によらない。"""
    checks, runs, summary = [], {}, {}
    travel, volumes = {}, {}
    for case, (section, law, inflow, max_depth, dt, fluids) in FLUME_CASES.items():
        sec = Section(section["shape"], section["size"])
        h_cap = sec.size if sec.shape == "semicircle" else max_depth
        for key in fluids:
            fl = get_fluid(key)
            hn = normal_depth(sec, inflow, slope, law, nu=fl.nu, h_max=h_cap, n_manning=n_manning)
            vn = inflow / sec.area(hn)
            t_pass = length / vn
            release = 3.0 * t_pass
            net = _sea(Network(dt, key))
            ch = net.channel("flume", "sea", length=length, cells=cells, slope=slope, section=section, law=law,
                             n_manning=n_manning, max_depth=max_depth)
            net.pump("sea", "flume:0", rate=inflow)
            net.tracer("somen", "flume:0", release=release, stop="sea")
            res = net.run(release + 2.0 * t_pass)
            info = res["channels"]["flume"]
            last = ch["cells"][-1]
            h_last = net.head(last)
            tr = res["tracers"][0]
            took = None if tr["arrival"] is None else tr["arrival"] - tr["path"][0]["t"]
            tag = f"{case}:{key}"
            travel[tag] = took
            volumes[tag] = [res["volumes_q"][c] for c in ch["cells"]]
            re = channel_reynolds(sec, hn, inflow, fl.nu)
            tol_t = 2 * dt / t_pass
            checks += [
                close(f"{tag}:depth", fl.icon + "📏", h_last, hn, 2e-3, note="末端セルの水深 = 等流水深"),
                close(f"{tag}:outflow", fl.icon + "🔁", res["final"][f"q:{cells - 1}"], inflow, 1e-3,
                      note="末端の流出 = 流入"),
                check(f"{tag}:somen", fl.icon + "🍜", took is not None and rel_err(took, t_pass) <= tol_t,
                      value=took, expected=t_pass, tol=tol_t, note="トレーサーが樋を抜ける時間 = L/v（ステップ内を補間、2·dt まで）"),
                check(f"{tag}:courant", "🧭", info["courant"] <= 1.0, value=info["courant"], expected=1.0,
                      note="Courant = c·dt/Δx ≤ 1（陽解法の安定）"),
                check(f"{tag}:no-overflow", "🫗", not info["overflow"], note="表の上限（縁）を超えていない"),
                check(f"{tag}:regime", "〰️", (re < 500) if law == "laminar" else (re > 500), value=re,
                      expected=500.0, note="Re = v·R/ν：層流なら < 500、Manning なら > 500"),
                conserved(f"{tag}:conserve", res),
            ]
            runs[tag] = res["series"]
            summary[tag] = {"normal_depth": hn, "velocity": vn, "t_pass": t_pass, "travel": took,
                            "reynolds": re, "courant": info["courant"]}
    w = get_fluid("water")
    for key in FLUID_KEYS:
        if key == "water":
            continue
        fl = get_fluid(key)
        a, b = f"rect-laminar:{key}", "rect-laminar:water"
        expected = summary[a]["t_pass"] / summary[b]["t_pass"]
        checks.append(close(f"{a}:speed", fl.icon + "⏱️", travel[a] / travel[b], expected, 1e-2,
                            note=f"層流は ν に効く（広い流れなら (ν/ν_water)^(1/3) = {(fl.nu / w.nu) ** (1 / 3):.4f}）"))
    checks.append(check("bamboo-manning:fluid-free", "🟰",
                        volumes["bamboo-manning:ethanol"] == volumes["bamboo-manning:water"],
                        note="Manning は粘度によらない → 🍶 と 💧 で各セルの体積がビット単位で同じ"))
    return {"params": {"length": length, "cells": cells, "slope": slope, "n_manning": n_manning,
                       "cases": {k: {"section": v[0], "law": v[1], "inflow": v[2], "max_depth": v[3], "dt": v[4],
                                     "fluids": list(v[5])} for k, v in FLUME_CASES.items()}},
            "summary": summary, "runs": runs, "checks": checks}


# --- E09 ドット絵の水（ビット盤） ---------------------------------------------

def e09_pixel_water(seed: int = 0, base_rate: float = 0.5, scale: int = 8, rows: int = pw.SECTION_ROWS) -> dict:
    """🟦 を下から積む断面図（行 = 高さ、列 = 区間・水槽）で、溜まる・流れる・遅れるが読めるかを確かめる。

    断面図: 段数 = round(h / 縦いっぱい × rows) で、段数から水深を ±½ 段まで読み戻せる。先端（最初の 1 段）が
    届く時刻がセルごとに増え、その間隔 ≈ Δx/v。じょうごには水が溜まって段が積み上がり、🍜 は管の遅れの間 ⏳。
    ドット絵: ドット数の保存・水面の平ら・流体ごとの広がりの順・再現性・PNG/GIF の読み戻し。
    """
    checks, runs, summary, settle = [], {}, {}, {}
    spacing = {}
    for key in ("water", "ethanol"):
        fl = get_fluid(key)
        res, cols = pw.section_case("flume", key)
        view = pw.network_frames(res, cols, rows=rows)
        arr = view["arrivals"]
        sec = Section("rect", 0.05)
        hn = normal_depth(sec, 2e-6, 0.05, "laminar", nu=fl.nu, h_max=0.005)
        dx_v = 0.1 / (2e-6 / sec.area(hn))
        gaps = [b - a for a, b in zip(arr[2:], arr[3:])] if None not in arr else []
        spacing[key] = sum(gaps) / len(gaps) if gaps else math.nan
        moves = [c for c in (w[0] for w in view["tracer_columns"]) if c is not None]
        checks += [
            check(f"section:{key}:readback", fl.icon + "📏", view["readback_error"] <= view["half_block"] + 1e-15,
                  value=view["readback_error"], expected=0.0, tol=view["half_block"],
                  note="段数 × 1 段の高さで水深を読み戻せる（±½ 段）"),
            check(f"section:{key}:front-order", fl.icon + "➡️",
                  None not in arr and all(a < b for a, b in zip(arr, arr[1:])),
                  value=arr[-1], note="最初の 1 段が積まれる時刻がセルごとに増える（流れ）"),
            close(f"section:{key}:front-gap", fl.icon + "⏱️", spacing[key], dx_v, 2e-2,
                  note="届く時刻の間隔 ≈ Δx/v（遅れ。v は等流の流速）"),
            check(f"section:{key}:somen-down", fl.icon + "🍜",
                  bool(moves) and moves == sorted(moves) and moves[-1] == 19,
                  value=moves[-1] if moves else None, expected=19, note="🍜 の列は下流へ進むだけ（戻らない）"),
            conserved(f"section:{key}:conserve", res),
        ]
        runs[f"section:{key}"] = [{"t": fr["t"], "lines": fr["lines"]} for fr in view["frames"]]
        summary[f"section:{key}"] = {"dx_over_v": dx_v, "front_gap": spacing[key], "arrivals": arr,
                                     "caps": view["caps"][:1], "rows": rows}
    checks.append(close("section:ethanol:slower", "🍶⏱️", spacing["ethanol"] / spacing["water"],
                        (get_fluid("ethanol").nu / get_fluid("water").nu) ** (1 / 3), 3e-2,
                        note="層流の樋は 🍶 の先端が遅い（広い流れなら (ν/ν_water)^(1/3)）"))
    res, cols = pw.section_case("flume-hopper", "water")
    view = pw.network_frames(res, cols, rows=rows, frames=None)
    hop = [cnt[0] for cnt in view["counts"]]
    pipe = next((a, b) for a, b, place in view["tracer_spans"][0] if place == "hopper>flume:0")
    pipe_col = view["columns"].index("hopper>flume:0")
    shown = [fr["t"] for fr in view["frames"] if any(line[pipe_col] == pw.IN_PIPE for line in fr["lines"])]
    first = view["arrivals"][view["columns"].index("flume:0")]
    checks += [
        check("section:hopper:fills", "🫙", hop == sorted(hop) and hop[-1] == rows, value=hop[-1], expected=rows,
              note="じょうごに水が溜まり、段が下から積み上がる（減らない）"),
        close("section:hopper:pipe-delay", "⏳", pipe[1] - pipe[0], 3.0, 1e-9, note="🍜 が管の中にいる時間 = 管の遅れ"),
        check("section:hopper:pipe-shown", "⏳", bool(shown) and abs(len(shown) * 0.05 - 3.0) <= 0.1,
              value=len(shown) * 0.05, expected=3.0, tol=0.1, note="⏳ が出る時刻の長さ ≈ 遅れ"),
        check("section:hopper:first-water", "🚿", first is not None and first >= 3.0, value=first, expected=3.0,
              note="樋の最初のセルに水が届くのは管の遅れの後"),
        check("section:hopper:readback", "📏", view["readback_error"] <= view["half_block"] + 1e-15,
              value=view["readback_error"], tol=view["half_block"], note="段数から水深を ±½ 段で読み戻せる"),
        conserved("section:hopper:conserve", res),
    ]
    runs["section:hopper"] = [{"t": fr["t"], "lines": fr["lines"]} for fr in view["frames"][:: max(1, len(view["frames"]) // 6)]]
    summary["section:hopper"] = {"pipe": pipe, "first_water": first, "hopper_blocks_final": hop[-1]}
    for scene in pw.SCENES:
        for key in FLUID_KEYS:
            fl = get_fluid(key)
            res = pw.run_pixels(scene, key, seed=seed, base_rate=base_rate)
            tag = f"{scene}:{key}"
            settle[tag] = res["settle_tick"]
            checks += [
                check(f"{tag}:count", fl.icon + "🔢", res["count_error"] == 0, value=res["count_error"], expected=0,
                      tol=0, note="毎ティックのドット数（popcount）が不変・壁と重ならない"),
                check(f"{tag}:flat", fl.icon + "🟰", res["settle_tick"] is not None and res["flatness"] <= 1,
                      value=res["flatness"], expected=1, note="静まった水面: 列の高さの最大 − 最小 ≤ 1"),
            ]
            if scene == "dam-break":
                # 毎ティック取り直して先端を読む。1 ティックに動ける列は 斜め 1 + 横 ceil(rate) まで
                fine = pw.run_pixels(scene, key, seed=seed, base_rate=base_rate, frame_every=1)
                arr = pw.pixel_arrivals(fine)[8:]
                slack = 1 + math.ceil(fine["lateral_rate"])
                checks.append(check(f"{tag}:front", fl.icon + "➡️", None not in arr and
                                    all(a <= b + slack for a, b in zip(arr, arr[1:])),
                                    value=arr[-1], tol=slack,
                                    note="水の先端が右の壁へ順に届く（1 ティックで飛び越せる列の分だけ前後を許す）"))
                settle[f"{tag}:wall"] = arr[-1]
            if scene == "moniac-drain":
                upper = sum(row.bit_count() for row in res["frames"][-1][0][:14])
                checks.append(check(f"{tag}:drained", fl.icon + "🕳️", upper == 0, value=upper, expected=0,
                                    note="上の槽は穴から全部落ちきる"))
            runs[tag] = [{"t": res["ticks"], "final": pw.render_text(res, -1, ascii_only=True)}]
            summary[tag] = {"pixels": res["pixels"], "settle_tick": res["settle_tick"],
                            "lateral_rate": res["lateral_rate"], "lateral_passes": res["lateral_passes"],
                            "first_floor_tick": res["first_floor_tick"], "hash": res["hash"]}
    nu_w = get_fluid("water").nu
    for scene in pw.SCENES:
        e, w_, m = (settle[f"{scene}:{k}"] for k in ("ethanol", "water", "methanol"))
        checks.append(check(f"{scene}:order", "⏱️", None not in (e, w_, m) and e > w_ > m,
                            value=[m, w_, e], note="静まるまで: 🧴 < 💧 < 🍶（ν が大きいほど横に広がる回数が減る）"))
        if scene == "dam-break":
            walls = [settle[f"dam-break:{k}:wall"] for k in ("methanol", "water", "ethanol")]
            checks.append(check("dam-break:wall-order", "🧱", None not in walls and walls[0] <= walls[1] <= walls[2],
                                value=walls, note="先端が右の壁に届く: 🧴 ≤ 💧 ≤ 🍶"))
        if None not in (e, w_, m):
            summary[f"{scene}:ratio"] = {"ethanol/water": e / w_, "methanol/water": m / w_,
                                         "nu ethanol/water": get_fluid("ethanol").nu / nu_w,
                                         "nu methanol/water": get_fluid("methanol").nu / nu_w}
    a = pw.run_pixels("dam-break", "water", seed=seed, base_rate=base_rate)
    b = pw.run_pixels("dam-break", "water", seed=seed, base_rate=base_rate)
    c = pw.run_pixels("dam-break", "water", seed=seed + 1, base_rate=base_rate)
    checks += [
        check("same-seed", "🎲", a["hash"] == b["hash"], value=a["hash"][:12], expected=b["hash"][:12],
              note="同じ seed → 全フレームが同じ"),
        check("other-seed", "🎲", a["hash"] != c["hash"], value=c["hash"][:12], note="seed を変えると動きが変わる"),
    ]
    short = pw.run_pixels("moniac-drain", "water", ticks=40, seed=seed, base_rate=base_rate)
    gif = pw.read_gif(pw.gif_bytes(short, scale=scale))
    gif_ok = len(gif["frames"]) == len(short["frames"]) and all(
        frame[(y * scale) * gif["width"] + x * scale] == row[x]
        for frame, i in zip(gif["frames"], range(len(short["frames"])))
        for y, row in enumerate(pw.indices(short, i)) for x in range(len(row)))
    png = pw.read_png(pw.png_bytes(short, -1, scale=scale))
    png_ok = (png["width"], png["height"]) == (short["width"] * scale, short["height"] * scale) and all(
        png["rows"][y * scale][x * scale] == row[x] for y, row in enumerate(pw.indices(short, -1)) for x in range(len(row)))
    checks += [
        check("gif-roundtrip", "🎞️", gif_ok, value=len(gif["frames"]), expected=len(short["frames"]),
              note="GIF89a を読み戻すと全フレームの色番号が一致（LZW・差分矩形）"),
        check("png-roundtrip", "🖼️", png_ok, value=png["chunks"], note="PNG の CRC・大きさ・色番号が一致"),
    ]
    return {"params": {"seed": seed, "base_rate": base_rate, "scale": scale, "scenes": list(pw.SCENES)},
            "summary": summary, "runs": runs, "checks": checks}


# --- E10 流しそうめん（竹の継ぎ目・段差・曲がり・引っかかり） ----------------------

BAMBOO = {"shape": "semicircle", "size": 0.04}  # 割った竹（半径 4 cm）
NOODLES = {"somen": 1.0e-3, "hiyamugi": 1.5e-3}  # 太さの目安。機械製麺の区分はそうめん < 1.3 mm ≤ ひやむぎ < 1.7 mm


def somen_layout(q: float, *, law: str = "manning", beta: float | None = None, dt: float = 0.01,
                 drop: float = 0.05, lip: float = 3e-3, bend_radius: float = 0.12, bucket: float = 0.5,
                 fluid: str = "water") -> Network:
    """桶 → ポンプ → 竹 3 本（直線 → 段差 → 曲がり → 縁の継ぎ目 → 直線）→ 受け皿 → ホースで桶へ戻す。
    曲がりや継ぎ目からあふれた水は床（floor）へ。"""
    net = Network(dt, fluid).tank("bucket", 0.05, volume=bucket).tank("basin", 0.05).tank("floor", 1.0)
    kw = {"law": law, "beta": beta, "section": BAMBOO}
    net.flume("take", "basin", [
        {"length": 2.0, "cells": 20, "slope": 0.05, "joint": {"drop": drop}, **kw},
        {"length": 1.0, "cells": 10, "slope": 0.03, "bend_radius": bend_radius, "joint": {"lip": lip}, **kw},
        {"length": 1.0, "cells": 10, "slope": 0.05, **kw},
    ], spill="floor")
    net.pump("bucket", "take/0:0", rate=q)
    net.pump("basin", "bucket", rate=q, delay=2.0)
    return net


def _into(res: dict, dst: str) -> float:
    """最後のステップに dst へ流れ込んだ流量の合計。"""
    return sum(res["final"][f"q:{ln['index']}"] for ln in res["links"] if ln["dst"] == dst and ln["kind"] != "pump")


def e10_somen_flume(q_flow: float = 5e-5, q_flood: float = 1.2e-3, q_trickle: float = 3e-6,
                    q_flush: float = 3e-5) -> dict:
    sec = Section(BAMBOO["shape"], BAMBOO["size"])
    rim = sec.size
    w = get_fluid("water")
    checks, runs, summary = [], {}, {}

    # 🎋 ふつうに流す: 等流水深・継ぎ目のたまり・段差の落下・あふれなし・追い越しなし
    net = somen_layout(q_flow)
    for i in range(3):
        net.tracer(f"somen-{i + 1}", "take/0:0", release=20.0 + i, stop="basin", size=NOODLES["somen"])
    res = net.run(40.0)
    h0 = normal_depth(sec, q_flow, 0.05, "manning", nu=w.nu, h_max=rim)
    pond = lip_head(sec, q_flow, 3e-3)
    drop_link = next(ln for ln in res["links"] if ln["src"] == "take/0:19")
    order = [tr["name"] for tr in sorted(res["tracers"], key=lambda t: t["arrival"] or math.inf)]
    chs = res["channels"]
    checks += [
        close("flow:depth", "📏", net.head("take/0:19"), h0, 2e-3, note="1 本目の竹の等流水深（Manning）"),
        close("flow:pond", "🫙", net.head("take/1:9"), pond, 2e-3, note="縁の継ぎ目の手前にたまる水深 = 縁 + 越流水深"),
        check("flow:drop", "🪂", abs(drop_link["delay"] - fall_time(0.05)) <= net.dt / 2 + 1e-12,
              value=drop_link["delay"], expected=fall_time(0.05), tol=net.dt / 2, note="段差の自由落下 √(2h/g)"),
        check("flow:no-spill", "🫗", chs["take/1"]["spilled"] == 0 and chs["take/1"]["max_outer_depth"] < rim,
              value=chs["take/1"]["max_outer_depth"], expected=rim, note="曲がりの外側の水面 < 竹の縁"),
        check("flow:order", "🍜", order == ["somen-1", "somen-2", "somen-3"] and all(t["state"] == "done"
                                                                                  for t in res["tracers"]),
              value=order, note="同じ太さのそうめんは放した順に届く（追い越しなし）"),
        check("flow:courant", "🧭", max(c["courant"] for c in chs.values()) <= 1.0,
              value=max(c["courant"] for c in chs.values()), expected=1.0, note="Courant ≤ 1"),
        conserved("flow:conserve", res),
    ]
    runs["flow"] = res["series"]
    summary["flow"] = {"normal_depth": h0, "pond_depth": pond, "arrivals": {t["name"]: t["arrival"]
                                                                           for t in res["tracers"]}}

    # 🌊 流しすぎ: 曲がりの外側と継ぎ目のたまりが縁を越えて床へあふれる
    net = somen_layout(q_flood)
    res = net.run(40.0)
    spilled = _into(res, "floor")
    checks += [
        check("flood:spill", "💦", spilled > 0, value=spilled, note="曲がり・継ぎ目から床へあふれる"),
        close("flood:balance", "🔁", _into(res, "basin") + spilled, q_flood, 1e-3, note="受け皿 + 床 = 流した量"),
        conserved("flood:conserve", res),
    ]
    runs["flood"] = res["series"]
    summary["flood"] = {"to_basin": _into(res, "basin"), "to_floor": spilled,
                        "max_outer_depth": res["channels"]["take/1"]["max_outer_depth"]}

    # 🥢 ちょろちょろ流す: 継ぎ目の越流水深より太いひやむぎだけ引っかかり、後から来たそうめんに抜かれる。
    #     途中で水を足すと流れ出す。細い流れなので層流の式（β = 3）で回す
    dt = 0.05
    net = somen_layout(q_trickle, law="laminar", beta=3.0, dt=dt)
    net.pump("bucket", "take/0:0", rate=q_flush, start=200.0)
    net.tracer("hiyamugi", "take/0:0", release=80.0, stop="basin", size=NOODLES["hiyamugi"])
    net.tracer("somen", "take/0:0", release=100.0, stop="basin", size=NOODLES["somen"])
    res = net.run(260.0)
    head = lip_head(sec, q_trickle, 3e-3) - 3e-3
    h_t = normal_depth(sec, q_trickle, 0.05, "laminar", nu=w.nu, h_max=rim, beta=3.0)
    re = channel_reynolds(sec, h_t, q_trickle, w.nu)
    tr = {t["name"]: t for t in res["tracers"]}
    checks += [
        check("trickle:premise", "📐", NOODLES["somen"] <= head < NOODLES["hiyamugi"], value=head,
              note="越流水深がそうめん以上・ひやむぎ未満になる流し方"),
        check("trickle:regime", "〰️", re < 500, value=re, expected=500.0, note="細い流れは層流（Re < 500）"),
        check("trickle:snag", "🪝", tr["hiyamugi"]["stuck_time"] > 0 and tr["hiyamugi"]["stuck_at"] == "take/1:9",
              value=tr["hiyamugi"]["stuck_time"], note="ひやむぎは縁の継ぎ目で引っかかる"),
        check("trickle:somen-free", "🍜", tr["somen"]["stuck_time"] == 0, value=tr["somen"]["stuck_time"],
              expected=0.0, note="そうめんは引っかからない"),
        check("trickle:overtake", "🏁", (tr["somen"]["arrival"] or math.inf) < (tr["hiyamugi"]["arrival"] or math.inf),
              value=[tr["somen"]["arrival"], tr["hiyamugi"]["arrival"]], note="後から放したそうめんが先に届く"),
        check("trickle:flush", "🚿", tr["hiyamugi"]["arrival"] is not None and tr["hiyamugi"]["arrival"] > 200.0,
              value=tr["hiyamugi"]["arrival"], expected=200.0, note="水を足すとひやむぎも流れ出す"),
        conserved("trickle:conserve", res),
    ]
    runs["trickle"] = res["series"]
    summary["trickle"] = {"head_over_lip": head, "reynolds": re,
                          "tracers": {k: {"arrival": v["arrival"], "stuck_time": v["stuck_time"]} for k, v in tr.items()}}
    return {"params": {"q_flow": q_flow, "q_flood": q_flood, "q_trickle": q_trickle, "q_flush": q_flush,
                       "bamboo": BAMBOO, "noodles": NOODLES},
            "summary": summary, "runs": runs, "checks": checks}


# --- E11 MONIAC 型の水回り（据え置きの試験台） ---------------------------------

RIG_OPEN = {"store": 0.5, "side": 0.3, "through": 0.2}  # 頭の水槽から分ける 3 つの弁の開き具合
RIG_STORE_CAM = [[0.0, 0.0], [0.05, 0.0], [0.25, 1.0]]  # 貯める槽の浮き → 自分の出口の弁（5 cm までは閉じたまま）
RIG_SIDE_CAM = [[0.0, 0.0], [0.05, 0.0], [0.15, 1.0]]  # 貯める槽の浮き → 横の槽の出口の弁（別の槽を動かす）


def moniac_rig(dt: float = 0.05, q1: float = 2e-4, q2: float = 3e-4, t_step: float = 600.0, stages: int = 1,
               fluid: str = "water") -> Network:
    """MONIAC の水回りを物理のまま組んだ試験台（経済の読み替えはしない）。

    🛢️ 下の槽 → ⛽ ポンプ → 頭の槽 → 3 つの弁で分流（貯める槽・横の槽・下へ素通り）。
    貯める槽の浮きがカムを通して、自分の出口の弁（遅れ stages 段）と横の槽の出口の弁を開く。
    貯める槽からの戻りは 3 秒の管。t_step でポンプの流量を q1 → q2 に上げる。
    """
    net = Network(dt, fluid).tank("sump", 1.0, volume=0.2).tank("header", 0.01).tank("store", 0.02).tank("side", 0.01)
    net.pump("sump", "header", rate=q1, until=t_step)
    net.pump("sump", "header", rate=q2, start=t_step)
    for name, o in RIG_OPEN.items():
        net.valve("header", "sump" if name == "through" else name, tau=10.0, opening=o)
    net.valve("store", "sump", tau=20.0, sensor="store", cam=RIG_STORE_CAM, stages=stages, link_tau=2.0, delay=3.0)
    net.valve("side", "sump", tau=15.0, sensor="store", cam=RIG_SIDE_CAM)
    return net


def _solve(f, lo: float, hi: float) -> float:
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if f(mid) < 0 else (lo, mid)
    return 0.5 * (lo + hi)


def rig_steady(q: float) -> dict:
    """試験台の定常（閉じた式）。頭の槽 h = Q·τ/(A·Σo)、分流は開き具合の比、
    貯める槽は cam(h)·A·h/τ = Q·o_store/Σo を二分法、横の槽は V = q_side·τ / cam_side(h_store)。"""
    total = sum(RIG_OPEN.values())
    flows = {k: q * o / total for k, o in RIG_OPEN.items()}
    store = _solve(lambda h: interp(RIG_STORE_CAM, h) * 0.02 * h / 20.0 - flows["store"], 0.0, 10.0)
    side = flows["side"] * 15.0 / min(1.0, interp(RIG_SIDE_CAM, store)) / 0.01
    return {"header": q * 10.0 / (0.01 * total), "store": store, "side": side, "flows": flows}


def interp(points, x: float) -> float:
    pts = sorted((float(a), float(b)) for a, b in points)
    if x <= pts[0][0]:
        return pts[0][1]
    if x >= pts[-1][0]:
        return pts[-1][1]
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        if x <= x1:
            return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
    return pts[-1][1]


def e11_moniac_rig(dt: float = 0.05, q1: float = 2e-4, q2: float = 3e-4, t_step: float = 600.0,
                   t_end: float = 1200.0) -> dict:
    checks, runs, summary = [], {}, {}
    net = moniac_rig(dt, q1, q2, t_step)
    res = net.run(t_end, trace=["store", "header"])
    links = res["links"]
    mid = int(round(t_step / dt)) - 1
    tr = res["trace"]
    tol = dt / 10.0  # 流入をステップ頭に入れる離散化の分（≈ dt/2τ、τ ≥ 10 s）を見込む
    for phase, q, idx in (("q1", q1, mid), ("q2", q2, len(tr["store"]) - 1)):
        st = rig_steady(q)
        checks += [
            close(f"{phase}:header", "🛢️", tr["header"][idx], st["header"], tol, note="頭の槽 h = Q·τ/(A·Σo)"),
            close(f"{phase}:store", "🫙", tr["store"][idx], st["store"], tol, note="貯める槽: 浮き → 弁のカムで決まる水位"),
        ]
        summary[phase] = {"steady": st, "header": tr["header"][idx], "store": tr["store"][idx]}
    st2 = rig_steady(q2)
    checks.append(close("q2:side", "🪣", res["final"]["h:side"], st2["side"], tol,
                        note="横の槽: 別の槽の浮きが出口の弁を開く"))
    f = res["final"]
    split = [f["q:2"], f["q:3"], f["q:4"]]
    want = [RIG_OPEN[k] / sum(RIG_OPEN.values()) for k in ("store", "side", "through")]
    checks.append(check("split", "🔱", all(rel_err(a / sum(split), b) <= 1e-8 for a, b in zip(split, want)),
                        value=[a / sum(split) for a in split], expected=want, tol=1e-8,
                        note="3 つの弁の分流 = 開き具合の比（床の丸め 1 ulp 分のずれまで）"))

    # 溜まる: 貯める槽の弁が開く前（浮きが 5 cm に届く前）は、入った量がそのまま貯まる（整数で一致）
    probe = moniac_rig(dt, q1, q2, t_step)
    early = probe.run(10.0)
    moved_in = early["links"][2]["volume_moved_q"]
    checks.append(check("fill:exact", "🫙", early["volumes_q"]["store"] == moved_in and
                        early["links"][5]["volume_moved_q"] == 0,
                        value=early["volumes_q"]["store"], expected=moved_in, tol=0,
                        note="弁が閉じている間: 貯める槽の体積 = 入った量（ulp まで）"))
    # 浮きが別の槽の弁を動かす: 横の槽の出口が開くのは、貯める槽の水位が 5 cm を越えたステップ
    cross = next(i for i, h in enumerate(tr["store"]) if interp(RIG_SIDE_CAM, h) > 0)
    opened = links[6]["first_arrival"]
    checks.append(check("float:cross", "🎈", opened is not None and abs(opened - (cross + 1) * dt) <= dt + 1e-12,
                        value=opened, expected=(cross + 1) * dt, tol=dt,
                        note="横の槽の弁は、貯める槽の浮きが 5 cm を越えた次のステップで開く"))
    # 遅れ: 貯める槽からの戻りは 3 秒の管を通る
    gap = links[5]["first_arrival"] - links[5]["first_flow"]
    checks.append(check("delay", "⏳", abs(gap - 3.0) <= dt / 2 + 1e-12, value=gap, expected=3.0, tol=dt / 2,
                        note="貯める槽の弁が開いてから下の槽に届くまで = 戻りの管の遅れ 3 秒"))
    late = window(tr["store"], 0.9)
    checks += [
        check("still", "🫧", peak_to_peak(late) < 1e-6, value=peak_to_peak(late), tol=1e-6,
              note="遅れ 1 段の浮き弁は揺れずに落ち着く"),
        check("step-up", "📈", summary["q2"]["store"] > summary["q1"]["store"], value=summary["q2"]["store"],
              expected=summary["q1"]["store"], note="ポンプを強めると貯める槽の水位が上がる"),
        conserved("conserve", res),
    ]
    runs["rig"] = res["series"]
    summary["golden"] = rig_fingerprint(res)
    return {"params": {"dt": dt, "q1": q1, "q2": q2, "t_step": t_step, "t_end": t_end, "open": RIG_OPEN,
                       "store_cam": RIG_STORE_CAM, "side_cam": RIG_SIDE_CAM},
            "summary": summary, "runs": runs, "checks": checks}


def rig_fingerprint(res: dict) -> dict:
    """据え置きの比較用: 各槽の最後の水位を有効 9 桁で（OS の libm の差で末尾が揺れても変わらない桁）。"""
    return {k: float(f"{v:.9g}") for k, v in sorted(res["final"].items()) if k.startswith("h:")}


# --- 登録簿 ---------------------------------------------------------------

EXPERIMENTS = {
    "E01": ("tank-integrator", "🧮", "水槽は積分器（整数で厳密）", e01_tank_integrator),
    "E02": ("fluid-drain", "💧", "流体ごとの排水（細管 / 小孔）", e02_fluid_drain),
    "E03": ("transport-delay", "⏳", "管を通る遅れ", e03_transport_delay),
    "E04": ("weir", "🪜", "堰とカム（非線形の出口）", e04_weir),
    "E05": ("float-valve", "🎈", "浮き弁と遅れ", e05_float_valve),
    "E06": ("coupled-loop", "🌐", "2 台連結（遅延つきの環）", e06_coupled_loop),
    "E07": ("fluid-speed", "🏎️", "作動流体で変わる回路の速さ", e07_fluid_speed),
    "E08": ("flume", "🎋", "傾いた樋とトレーサー（流しそうめんの下ごしらえ）", e08_flume),
    "E09": ("pixel-water", "🟦", "🟦 を積む断面図で溜まる・流れる・遅れる（ドット絵の水も）", e09_pixel_water),
    "E10": ("somen-flume", "🥢", "流しそうめん（竹の継ぎ目・段差・曲がり・引っかかり）", e10_somen_flume),
    "E11": ("moniac-rig", "🏛️", "MONIAC 型の水回り（分流・貯める・浮きで弁・遅れ）を据え置きの試験に", e11_moniac_rig),
}
SLUGS = {slug: eid for eid, (slug, *_rest) in EXPERIMENTS.items()}


def resolve(key: str) -> str:
    """'E02'・'e02'・'fluid-drain'・'E02:fluid-drain' のどれでも ID に直す。"""
    k = key.strip()
    if ":" in k:
        k = k.split(":", 1)[0]
    if k.upper() in EXPERIMENTS:
        return k.upper()
    if k in SLUGS:
        return SLUGS[k]
    raise KeyError(f"unknown experiment {key!r}")


def list_experiments() -> list[dict]:
    return [{"id": eid, "slug": slug, "icon": icon, "title": title} for eid, (slug, icon, title, _f) in EXPERIMENTS.items()]


def run_experiment(key: str, **overrides) -> dict:
    """実験を回す。戻り値は JSON 越しでも変わらない形（Node から ``--call run_experiment`` で呼べる）。"""
    eid = resolve(key)
    slug, icon, title, func = EXPERIMENTS[eid]
    out = func(**overrides)
    checks = out["checks"]
    return jsonable({"id": eid, "slug": slug, "icon": icon, "title": title, "ok": all(c["ok"] for c in checks),
                     "passed": sum(c["ok"] for c in checks), "total": len(checks), **out})


def verify_all(keys=None) -> dict:
    results = [run_experiment(k) for k in (keys or EXPERIMENTS)]
    return {"ok": all(r["ok"] for r in results), "results": results}


# --- 表示 -----------------------------------------------------------------

def _fmt(x) -> str:
    if isinstance(x, float):
        return f"{x:.6g}"
    return str(x)


def render(result: dict, verbose: bool = False) -> str:
    mark = "✅" if result["ok"] else "❌"
    lines = [f"{mark} {result['id']} {result['icon']} {result['slug']}  {result['passed']}/{result['total']}"]
    for c in result["checks"]:
        if verbose or not c["ok"]:
            m = "✅" if c["ok"] else "❌"
            vals = ""
            if c["value"] is not None:
                vals = f"  {_fmt(c['value'])}"
                if c["expected"] is not None:
                    vals += f" ≟ {_fmt(c['expected'])}"
                if c["delta"] is not None:
                    vals += f"  δ={_fmt(c['delta'])}"
            lines.append(f"   {m} {c['icon']} {c['id']}{vals}")
    return "\n".join(lines)


def render_fluids() -> str:
    rows = [f"{f['icon']} {pw.EMOJI.get(f['key'], '🟦')} {f['key']:<15} ρ={f['rho']:g} kg/m³  μ={f['mu'] * 1e3:g} mPa·s"
            f"  ν={f['nu'] * 1e6:.4g} mm²/s" + ("  🏷️ 仮置き（物性は水）" if f["placeholder"] else "")
            for f in fluid_table()]
    return "\n".join(rows)


def write_jsonl(path: str, results: list[dict]) -> None:
    """プロッター記録（runs の各行）を JSONL で書く。1 行 = 1 サンプル。"""
    with Path(path).open("w", encoding="utf-8") as fh:
        for r in results:
            for run, rows in r["runs"].items():
                for row in rows:
                    fh.write(json.dumps({"id": r["id"], "slug": r["slug"], "run": run, **row}, ensure_ascii=False))
                    fh.write("\n")


def rig_section(scene: str, fluid: str = "water") -> tuple[dict, list[str]]:
    """断面図で見る試験台（run の結果, 列）。moniac-rig は頭・貯める・戻りの管・横の槽（下の槽は水源なので省く）、
    somen-flume は竹 3 本と受け皿・床。"""
    if scene == "moniac-rig":
        return moniac_rig(fluid=fluid).run(1200.0), ["header", "store", "store>sump", "side"]
    if scene == "somen-flume":
        net = somen_layout(5e-5, fluid=fluid)
        for i in range(3):
            net.tracer(f"somen-{i + 1}", "take/0:0", release=20.0 + 2 * i, stop="basin", size=NOODLES["somen"])
        return net.run(40.0, sample_every=20), ["take/0", "take/1", "take/2", "basin", "floor"]
    raise KeyError(f"unknown rig scene {scene!r}")


RIG_SCENES = ("moniac-rig", "somen-flume")


def _draw(args) -> int:
    if args.scene in RIG_SCENES:
        if args.view == "files":
            build_parser().error(f"{args.scene} は --view frames だけ")
        res, cols = rig_section(args.scene, args.fluid)
        print(pw.render_frames(pw.network_frames(res, cols, frames=args.frames, rows=args.rows,
                                                 cap_mode=args.cap_mode, fill=args.fill)))
        return 0 if res["conservation_error_q"] == 0 else 1
    if args.scene in pw.SECTION_CASES:
        if args.view == "files":
            build_parser().error(f"{args.scene} は --view frames だけ")
        res, cols = pw.section_case(args.scene, args.fluid)
        print(pw.render_frames(pw.network_frames(res, cols, frames=args.frames, rows=args.rows,
                                                 cap_mode=args.cap_mode, fill=args.fill)))
        return 0 if res["conservation_error_q"] == 0 else 1
    if args.view == "files":
        formats = ("gif", "png") if args.format == "both" else (args.format,)
        out = pw.draw(args.scene, args.fluid, args.out_dir, ticks=args.ticks, seed=args.seed, scale=args.scale,
                      frame_every=args.frame_every, formats=formats, file_prefix=args.file_prefix,
                      file_suffix=args.file_suffix)
        icons = {"gif": "🎞️", "png": "🖼️", "json": "🧾"}
        for kind, path in out["paths"].items():
            print(f"{icons[kind]} {path}")
        print(f"🔢 {out['pixels']} px  🎬 {out['frames']}  📐 {out['canvas']}  ⏱️ {out['settle_tick']}  ⚖️ {out['count_error']}")
        return 0 if out["count_error"] == 0 else 1
    res = pw.run_pixels(args.scene, args.fluid, ticks=args.ticks, seed=args.seed, frame_every=args.frame_every)
    frames = 1 if args.view == "frame" else args.frames
    print(pw.render_frames(pw.pixel_frames(res, frames=frames, fill=args.fill)))
    return 0 if res["count_error"] == 0 else 1


def build_parser() -> argparse.ArgumentParser:
    choices = sorted(EXPERIMENTS) + sorted(SLUGS)
    p = argparse.ArgumentParser(prog="fluidfullspeedsimulation", description="🌊 流れ・遅れ・検証")
    sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("list", help="🗂️ 実験の一覧")
    sub.add_parser("fluids", help="🧪 流体の物性")
    run = sub.add_parser("run", help="▶️ 実験を回す")
    run.add_argument("experiment", choices=choices, help="ID か識別子")
    ver = sub.add_parser("verify", help="✅ 全実験を回して検証")
    drw = sub.add_parser("draw", help="🎨 🟦 を積む断面図で溜まる・流れる・遅れるを見る（ドット絵・GIF・PNG も）")
    drw.add_argument("scene", choices=RIG_SCENES + pw.SECTION_CASES + pw.SCENES,
                     help="場面（試験台・樋の回路・ドット絵）")
    drw.add_argument("--view", choices=("frames", "frame", "files"), default="frames",
                     help="frames = 時刻ごとの断面図 / frame = 最後の 1 枚 / files = GIF・PNG（ドット絵だけ）")
    drw.add_argument("--fluid", choices=ALL_FLUID_KEYS, default="water", help="作動流体（🏷️ 仮置きも選べる）")
    drw.add_argument("--fill", help="水のブロック（色は自由。既定は流体ごと: 🟦 水 🟨 エタノール 🟧 メタノール "
                                    "🟩 Namekian Water 🟪 Poison Water）")
    drw.add_argument("--frames", type=int, default=6, help="並べる時刻の数")
    drw.add_argument("--rows", type=int, default=pw.SECTION_ROWS, help="断面図の縦の段数（既定 8）")
    drw.add_argument("--cap-mode", choices=pw.CAP_MODES, default="run-max",
                     help="縦いっぱいの水深: run-max = 走らせた中の最大 / rim = 樋の縁")
    drw.add_argument("--ticks", type=int, help="ドット絵のティック数（既定は場面ごと）")
    drw.add_argument("--seed", type=int, default=0, help="ドット絵の向きの初期値の乱数 seed")
    drw.add_argument("--scale", type=int, default=8, help="1 マスのピクセル数（既定 8。Ironmate の scale と同じ意味）")
    drw.add_argument("--frame-every", type=int, default=2, help="何ティックごとに 1 フレーム")
    drw.add_argument("--format", choices=("gif", "png", "both"), default="both", help="files で書き出す形式")
    drw.add_argument("--out-dir", default="build/pixel", help="書き出し先")
    drw.add_argument("--file-prefix", default="", help="ファイル名の頭に付ける文字列")
    drw.add_argument("--file-suffix", default="", help="ファイル名の末尾（拡張子の前）に付ける文字列")
    for sp in (run, ver):
        sp.add_argument("--json", action="store_true", help="結果を JSON で出す")
        sp.add_argument("--jsonl", metavar="PATH", help="プロッター記録を JSONL へ")
        sp.add_argument("--verbose", "-v", action="store_true", help="通った検証も並べる")
    return p


def _utf8_stdout() -> None:
    """絵文字を出すので、UTF-8 でない端末（Windows のパイプ等）では stdout を UTF-8 に切り替える。"""
    enc = (getattr(sys.stdout, "encoding", None) or "").lower().replace("-", "")
    if enc != "utf8" and hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    _utf8_stdout()
    if args.command == "list":
        for e in list_experiments():
            print(f"{e['icon']} {e['id']} {e['slug']}  {e['title']}")
        return 0
    if args.command == "fluids":
        print(render_fluids())
        return 0
    if args.command == "draw":
        return _draw(args)
    results = [run_experiment(args.experiment)] if args.command == "run" else verify_all()["results"]
    if args.jsonl:
        write_jsonl(args.jsonl, results)
    if args.json:
        slim = [{k: v for k, v in r.items() if k != "runs"} for r in results]
        print(json.dumps(slim[0] if args.command == "run" else slim, ensure_ascii=False, indent=2))
    else:
        for r in results:
            print(render(r, args.verbose))
    return 0 if all(r["ok"] for r in results) else 1


if __name__ == "__main__":
    sys.exit(main())
