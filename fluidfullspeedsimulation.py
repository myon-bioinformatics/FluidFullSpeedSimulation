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
import importlib.util
import json
import math
import random
import sys
from pathlib import Path

from fixedq import ONE, from_q, lag_step, qmul, shift_lag_step, to_q
from fluids import FLUIDS, REAL_FLUIDS, RE_LAMINAR, Section, channel_reynolds, fall_time, fluid_table, get_fluid, lip_head
from fluids import TEMP_TABLES, lip_flow, normal_depth, superelevation, table_value
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
# 太さの目安。機械製麺の区分はそうめん < 1.3 mm ≤ ひやむぎ < 1.7 mm ≤ うどん
NOODLES = {"somen": 1.0e-3, "hiyamugi": 1.5e-3, "udon": 2.0e-3}


def somen_layout(q: float, *, law: str = "manning", beta: float | None = None, dt: float = 0.01,
                 drop: float = 0.05, lip: float = 3e-3, bend_radius: float = 0.12, bucket: float = 0.5,
                 fluid: str = "water", temp_c: float | None = None, tail_cells: int = 10) -> Network:
    """桶 → ポンプ → 竹 3 本（直線 → 段差 → 曲がり → 縁の継ぎ目 → 直線）→ 受け皿 → ホースで桶へ戻す。
    曲がりや継ぎ目からあふれた水は床（floor）へ。temp_c で水温を変えられる（🧊 E15）。"""
    net = Network(dt, fluid, temp_c=temp_c).tank("bucket", 0.05, volume=bucket).tank("basin", 0.05)
    net.tank("floor", 1.0)
    kw = {"law": law, "beta": beta, "section": BAMBOO}
    net.flume("take", "basin", [
        {"length": 2.0, "cells": 20, "slope": 0.05, "joint": {"drop": drop}, **kw},
        {"length": 1.0, "cells": 10, "slope": 0.03, "bend_radius": bend_radius, "joint": {"lip": lip}, **kw},
        {"length": 1.0, "cells": tail_cells, "slope": 0.05, **kw},
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


# --- E12 長さのある麺（体が区間をまたぐ・縁に掛かる・曲がりで飛び出す） ------------------

STRAND = {"size": 1.5e-3, "density": 1100.0, "friction": 0.3}  # ひやむぎ程度の太さ。密度・摩擦は仮定の値


def bend_outer(q: float, slope: float = 0.03, radius: float = 0.12, n_manning: float = 0.010) -> float:
    """曲がりの等流で、外側の壁の水面 h + Δh/2（Δh = v²B/(gR)）。"""
    sec = Section(BAMBOO["shape"], BAMBOO["size"])
    h = normal_depth(sec, q, slope, "manning", nu=get_fluid("water").nu, h_max=sec.size, n_manning=n_manning)
    v = q / sec.area(h)
    return h + superelevation(v, sec.top_width(h), radius) / 2.0


def eject_flow(size: float, lo: float = 5e-5, hi: float = 2e-3) -> float:
    """外側の水面が「縁 − 太さ/2」になる流量（水はこぼれず、麺だけが縁に届く窓のまん中）。"""
    target = BAMBOO["size"] - size / 2.0
    return _solve(lambda q: bend_outer(q) - target, lo, hi)


def drape_flume(q: float = 3e-6, splash: tuple | None = None, dt: float = 0.01) -> Network:
    """細い流れ（層流 β = 3）の竹 2 本を縁 3 mm の継ぎ目でつなぐ（区間 5 cm）。
    splash = (始め, 終わり, 流量): 継ぎ目の手前のセルへ手で水をかけ、頭だけ縁を越えさせる。"""
    net = Network(dt).tank("bucket", 0.05, volume=0.5).tank("basin", 0.05)
    kw = {"law": "laminar", "beta": 3.0, "section": BAMBOO}
    net.flume("take", "basin", [
        {"length": 2.0, "cells": 40, "slope": 0.05, "joint": {"lip": 3e-3}, **kw},
        {"length": 1.0, "cells": 20, "slope": 0.05, **kw},
    ])
    net.pump("bucket", "take/0:0", rate=q)
    net.pump("basin", "bucket", rate=q, delay=2.0)
    if splash:
        t0, t1, rate = splash
        net.pump("bucket", "take/0:39", rate=rate, start=t0, until=t1)
    return net


def e12_noodle_length(length: float = 0.3, long_drape: float = 0.5) -> dict:
    """🍜 長さのある麺。体は頭の後ろに通ってきた区間を合計 length まで引く（区間の分解能で）。

    - 体の区間数 = ceil(length / Δx)、長さは速さを変えない（水と一緒に先入れ先出し）、尾は L/v 遅れて通る
    - 体のどこかが「越流水深 < 太さ」の縁に掛かっていれば全体が止まる（頭が越えた後でも尾で引っかかる）
    - 曲がりでは外側の水面 + 太さが縁に届くと、麺だけ縁の外へはじき出される（水はまだこぼれない）
    """
    checks, runs, summary = [], {}, {}
    sec = Section(BAMBOO["shape"], BAMBOO["size"])
    w = get_fluid("water")

    # 🍜 ふつうに流す: 体の区間数・長さで速さが変わらない・尾の遅れ
    net = somen_layout(5e-5)
    net.tracer("long", "take/0:0", release=20.0, stop="basin", length=length, **STRAND)
    net.tracer("point", "take/0:0", release=20.0, stop="basin", size=STRAND["size"])
    res = net.run(40.0, sample_every=5)
    tr = {t["name"]: t for t in res["tracers"]}
    dx = 0.1
    want_cells = math.ceil(length / dx - 1e-9)
    spans = [len(r["bodies"]["long"]) for r in res["series"] if "long" in r.get("bodies", {})
             and all(c.startswith("take/0:") for c in r["bodies"]["long"])]
    enter = {p["node"]: p["t"] for p in tr["long"]["path"]}
    h0 = normal_depth(sec, 5e-5, 0.05, "manning", nu=w.nu, h_max=sec.size)
    v0 = 5e-5 / sec.area(h0)
    lag = enter["take/0:15"] - enter["take/0:12"]  # 頭が 3 区間（= 0.3 m）進む時間 = 尾がそこを通る遅れ
    view = pw.network_frames(res, ["take/0"], frames=None)
    wide = max((sum(line.count(pw.CARRIED) for line in fr["lines"]) for fr in view["frames"]), default=0)
    checks += [
        check("span", "📏", bool(spans) and max(spans) == want_cells, value=max(spans) if spans else None,
              expected=want_cells, note="体の区間数 = ceil(長さ / Δx)"),
        check("same-speed", "🟰", tr["long"]["arrival"] is not None and
              abs(tr["long"]["arrival"] - tr["point"]["arrival"]) <= net.dt + 1e-12,
              value=tr["long"]["arrival"], expected=tr["point"]["arrival"], tol=net.dt,
              note="長さは速さを変えない（水と一緒に流れる）"),
        close("tail-lag", "⏱️", lag, length / v0, 2e-2, note="尾は頭から L/v 遅れて同じ所を通る"),
        check("section-span", "🍜", wide >= want_cells, value=wide, expected=want_cells,
              note="断面図で 🍜 が体の区間すべてに出る"),
        conserved("flow:conserve", res),
    ]
    runs["flow"] = res["series"]
    summary["flow"] = {"span_cells": want_cells, "lag": lag, "L_over_v": length / v0}

    # 🪝 尾で引っかかる: 細い流れで頭も尾も縁で待つ → 継ぎ目に 1 秒だけ水をかけて頭を越えさせる → 止めると尾が残る
    splash = (150.0, 151.0, 1e-6)
    net = drape_flume(splash=splash, dt=0.02)
    net.tracer("long", "take/0:0", release=80.0, stop="basin", length=long_drape, **STRAND)
    net.tracer("point", "take/0:0", release=80.0, stop="basin", size=STRAND["size"])
    res = net.run(175.0, sample_every=50)
    tr = {t["name"]: t for t in res["tracers"]}
    waited = [next((p["t"] for p in t["path"] if p["node"] == "take/0:39"), math.inf) for t in tr.values()]
    head = tr["long"]["path"][-1]["node"]
    pond_before = lip_head(sec, 3e-6, 3e-3) - 3e-3
    pond_splash = lip_head(sec, 3e-6 + splash[2], 3e-3) - 3e-3
    checks += [
        check("drape:premise", "📐", pond_before < STRAND["size"] <= pond_splash, value=[pond_before, pond_splash],
              expected=STRAND["size"], note="ふだんの越流水深 < 太さ ≤ 水をかけている間の越流水深"),
        check("drape:wait", "⏳", max(waited) < splash[0], value=max(waited), expected=splash[0],
              note="水をかける前に、どちらも縁の手前まで来ている"),
        check("drape:point-free", "🍜", tr["point"]["state"] == "done", value=tr["point"]["arrival"],
              note="点の麺（長さ 0）は水をかけた間に越えて届く"),
        check("drape:tail-held", "🪝", tr["long"]["state"] == "node" and tr["long"]["stuck_at"] == "take/0:39"
              and head.startswith("take/1:"), value=[head, tr["long"]["stuck_at"]],
              note="長い麺は頭が越えた後、縁に掛かった尾で止まる"),
        conserved("drape:conserve", res),
    ]
    runs["drape"] = res["series"]
    summary["drape"] = {"head_over_rim": pond_before, "head_over_rim_splash": pond_splash, "long_head": head,
                        "long_stuck_time": tr["long"]["stuck_time"]}

    # 💨 曲がりで飛び出す: 外側の水面が「縁 − 太さ/2」になる流量で、麺だけ外へ（水は曲がりからこぼれない）
    q_fast = eject_flow(STRAND["size"])
    outer = bend_outer(q_fast)
    net = somen_layout(q_fast, bucket=1.0)
    net.tracer("long", "take/0:0", release=10.0, stop="basin", length=length, **STRAND)
    res = net.run(20.0, sample_every=20)
    t_long = res["tracers"][0]
    bend_spill = sum(ln["volume_moved"] for ln in res["links"] if ln["kind"] == "spill"
                     and ln["src"].startswith("take/1:") and ln["src"] != "take/1:9")
    measured = max(ln["max_outer"] for ln in res["links"] if ln["kind"] == "spill"
                   and ln["src"].startswith("take/1:") and ln["src"] != "take/1:9")
    checks += [
        check("eject:premise", "📐", outer < BAMBOO["size"] <= outer + STRAND["size"], value=outer,
              expected=BAMBOO["size"], note="外側の水面 < 縁 ≤ 外側の水面 + 太さ"),
        check("eject:out", "💨", t_long["state"] == "spilled" and (t_long["spilled_at"] or "").startswith("take/1:"),
              value=t_long["spilled_at"], note="麺は曲がりで縁の外（床）へはじき出される"),
        check("eject:water-stays", "🫗", bend_spill == 0.0, value=bend_spill, expected=0.0,
              note="曲がりの区間からは水はこぼれていない（麺だけ）"),
        close("eject:outer", "🌀", measured, outer, 2e-3, note="曲がりの外側の水面（シミュレーションの最大）= 閉じた式"),
        conserved("eject:conserve", res),
    ]
    summary["eject"] = {"q": q_fast, "outer_depth": outer, "measured_outer_max": measured,
                        "spilled_at": t_long["spilled_at"]}
    runs["eject"] = res["series"]
    return {"params": {"length": length, "long_drape": long_drape, "strand": STRAND, "splash": list(splash)},
            "summary": summary, "runs": runs, "checks": checks}


# --- E13 麺のたまり（縁に積もってせき止め、一気に押し流される） ----------------------

def clog_prediction(q: float, strands: int, size: float, length: float, density: float, friction: float,
                    lip: float = 3e-3, cd: float = 0.62, fluid: str = "water") -> dict:
    """閉じた式。縁の幅 w = 縁の高さでの水面幅（たまりの間も一定とみなす）。

    - 越流水深 head = (Q/K)^(2/3)、K = Cd·(2/3)·w·√(2g)（頂が上がっても幅は同じなので head は変わらない）
    - 1 束が積もる高さ Δ₁ = 本数·πd²/4 / (PACKING·w)、k 束で Δ_k = k·Δ₁、たまりの水深 H_k = lip + Δ_k + head
    - 押し流し: ρ·g·w·Δ·(H − lip − Δ/2) ≥ Σ μ·(ρ_麺 − ρ)·g·a·L ⇔ Δ/2 + head ≥ h* = μ·PACKING·(ρ_麺/ρ − 1)·L
      → 流れ出す束の数 k* = ceil(2·(h* − head) / Δ₁)
    """
    from network import PACKING
    sec = Section(BAMBOO["shape"], BAMBOO["size"])
    w = sec.top_width(lip)
    g = 9.80665
    k = cd * (2.0 / 3.0) * w * math.sqrt(2.0 * g)
    head = (q / k) ** (2.0 / 3.0)
    d1 = strands * math.pi * size ** 2 / 4.0 / (PACKING * w)
    rho = get_fluid(fluid).rho
    h_star = friction * PACKING * (density / rho - 1.0) * length
    k_star = max(1, math.ceil(2.0 * (h_star - head) / d1 - 1e-12))
    return {"width": w, "head": head, "delta_1": d1, "h_star": h_star, "k_star": k_star,
            "pond": [lip + i * d1 + head for i in range(k_star + 1)]}


def e13_noodle_clog(q: float = 3e-6, every: float = 5.0, bunches: int = 34, strands: int = 10,
                    length: float = 0.2, first: float = 80.0, dt: float = 0.05,
                    friction: float = STRAND["friction"], density: float = STRAND["density"]) -> dict:
    """🧱 一口分の束（ひやむぎ 10 本）を every 秒ごとに流す。越流水深が太さより浅いので縁で引っかかり、
    積もった分だけ縁が高くなって水が溜まる。たまりが押す力が摩擦を超えると、全部まとめて流れ出す。"""
    net = somen_layout(q, law="laminar", beta=3.0, dt=dt)
    for i in range(bunches):
        net.tracer(f"bunch-{i + 1:02d}", "take/0:0", release=first + i * every, stop="basin", length=length,
                   strands=strands, size=STRAND["size"], density=density, friction=friction)
    t_end = first + bunches * every + 60.0
    res = net.run(t_end, sample_every=20, trace=["take/1:9"])
    pred = clog_prediction(q, strands, STRAND["size"], length, density, friction)
    lip = res["lips"][0]
    rel = lip["releases"]
    depth = res["trace"]["take/1:9"]
    tr = {t["name"]: t for t in res["tracers"]}
    joined = sorted((t["held_since"], name) for name, t in tr.items() if t["held_since"] is not None)
    checks, summary = [], {"prediction": pred, "releases": rel}
    # たまりの水深: k 束が積もって落ち着いたところ（次の束が来る直前）= lip + Δ_k + head
    ladder = []
    first_release = rel[0]["t"] if rel else math.inf
    for k in range(1, pred["k_star"]):
        t_next = joined[k][0]
        if t_next >= first_release:
            break
        ladder.append((k, depth[int(round(t_next / dt)) - 2], pred["pond"][k]))  # 次の束が積もる直前
    worst = max((rel_err(got, want) for _k, got, want in ladder), default=math.inf)
    checks += [
        check("held", "🪝", pred["head"] < STRAND["size"], value=pred["head"], expected=STRAND["size"],
              note="越流水深 < 太さ（束は縁で引っかかる）"),
        check("pond-ladder", "🪜", len(ladder) == pred["k_star"] - 1 and worst <= 5e-3, value=worst, tol=5e-3,
              note="k 束積もったたまりの水深 = 縁 + Δ_k + 越流水深（k = 1 … k*−1。次の束まで 5 秒で落ち着く途中の分を見込む）"),
        check("release-count", "🧮", bool(rel) and all(len(r["members"]) == pred["k_star"] for r in rel),
              value=[len(r["members"]) for r in rel], expected=pred["k_star"],
              note="流れ出すのは k* = ceil(2(h* − head)/Δ₁) 束たまったとき"),
        check("release-force", "⚖️", bool(rel) and all(r["push"] >= r["hold"] for r in rel),
              value=[round(r["push"] / r["hold"], 4) for r in rel], expected=1.0,
              note="流れ出す瞬間: 押す力 ≥ 支える力（記録の値）"),
        check("sawtooth", "📈", len(rel) >= min(2, bunches // (2 * pred["k_star"] + 4) + 1), value=len(rel),
              expected=2, note="溜まる → 一気に流れる、を繰り返す（束が足りるだけ）"),
    ]
    spreads, drops = [], []
    for r in rel:
        arr = [tr[m]["arrival"] for m in r["members"]]
        spreads.append(max(arr) - min(arr) if None not in arr else math.inf)
        i0 = int(round(r["t"] / dt)) - 1
        after = depth[i0: i0 + int(round(every / dt))]
        drops.append(depth[i0 - 1] - min(after))
    checks += [
        check("bunched", "🍜", bool(spreads) and max(spreads) < every / 2, value=max(spreads, default=None),
              expected=every, note=f"押し流された束は固まって届く（間隔 {every} 秒より詰まる）"),
        check("drop", "📉", bool(drops) and min(drops) >= 0.5 * pred["k_star"] * pred["delta_1"],
              value=min(drops, default=None), expected=pred["k_star"] * pred["delta_1"],
              note="流れ出すとたまりの水深がたまりの高さの半分以上下がる"),
        check("no-overflow", "🫗", max(depth) < BAMBOO["size"], value=max(depth), expected=BAMBOO["size"],
              note="たまりは竹の縁を越えない"),
        conserved("conserve", res),
    ]
    swept = sorted(name for name, t in tr.items() if t["state"] == "done" and t["stuck_time"] == 0)
    summary.update({"ladder": ladder, "spreads": spreads, "drops": drops, "swept_along": swept})
    return {"params": {"q": q, "every": every, "bunches": bunches, "strands": strands, "length": length,
                       "first": first, "dt": dt, "size": STRAND["size"], "density": density, "friction": friction},
            "summary": summary, "runs": {"clog": res["series"]}, "checks": checks}


# --- E14 箸ですくう（遅い麺ほど取れる） --------------------------------------------

CATCH_ZONE = "take/2:1"  # 最後の竹（4 区間 × 25 cm）の 2 つ目の区間に手を構える
GRAB = (0.3, 0.8)  # 手を伸ばして挟むまでの時間 [s] の幅（一様。仮定の値）


def catch_run(q: float, *, noodles: int = 30, every: float = 1.0, cooldown: float = 0.0, seed: int = 7,
              grab=GRAB, dt: float = 0.02) -> dict:
    net = somen_layout(q, tail_cells=4, dt=dt)
    for i in range(noodles):
        net.tracer(f"somen-{i + 1:02d}", "take/0:0", release=20.0 + i * every, stop="basin",
                   size=NOODLES["somen"])
    net.catcher("hand", CATCH_ZONE, grab=grab, cooldown=cooldown, seed=seed)
    return net.run(20.0 + noodles * every + 20.0, sample_every=50)


def e14_chopsticks(flows=(2e-5, 5e-5, 1.2e-4, 2e-4), noodles: int = 30, seed: int = 7) -> dict:
    """🥢 取り手は CATCH_ZONE の区間で、入ってきた麺に t_grab（麺ごとに GRAB の一様乱数、seed で再現）
    たってから箸を入れる。麺がその前に区間を出れば取り逃し。区間にいる時間 = Δx / v なので、流れが速いほど
    取れない。取れる割合の目安は P(t_grab ≤ Δx/v) = (Δx/v − 0.3) / 0.5（0..1 に収める）。"""
    sec = Section(BAMBOO["shape"], BAMBOO["size"])
    w = get_fluid("water")
    dx = 0.25
    checks, runs, summary = [], {}, {}
    caught_counts = []
    for q in flows:
        res = catch_run(q, noodles=noodles, seed=seed)
        c = res["catchers"][0]
        h = normal_depth(sec, q, 0.05, "manning", nu=w.nu, h_max=sec.size)
        stay = dx / (q / sec.area(h))
        draws = random.Random(seed)
        expect = sum(draws.uniform(*GRAB) <= stay for _ in range(noodles))
        exact = all(abs(x["t_catch"] - x["t_in"] - x["t_grab"]) <= 1e-9 for x in c["caught"]) and \
            all(m["why"] == "fast" and m["residence"] < m["t_grab"] for m in c["missed"])
        tag = f"q{q:.1e}"
        re = channel_reynolds(sec, h, q, w.nu)
        checks += [
            check(f"{tag}:rule", "📏", exact, value=[len(c["caught"]), len(c["missed"])],
                  note="取れた麺は入ってちょうど t_grab 後、取り逃した麺は t_grab より先に区間を出ている"),
            check(f"{tag}:count", "🥢", abs(len(c["caught"]) - expect) <= 1, value=len(c["caught"]), expected=expect,
                  tol=1, note="取れた数 ≈ #(t_grab ≤ Δx/v)（等流の v で。境目の 1 本は時刻の刻みで入れ替わりうる）"),
            check(f"{tag}:all-accounted", "🧾", len(c["caught"]) + len(c["missed"]) == noodles,
                  value=len(c["caught"]) + len(c["missed"]), expected=noodles, note="全部の麺が取れたか取り逃しか"),
            check(f"{tag}:regime", "〰️", re > 500, value=re, expected=500.0, note="Manning の前提（Re > 500）"),
            conserved(f"{tag}:conserve", res),
        ]
        caught_counts.append(len(c["caught"]))
        summary[tag] = {"q": q, "stay": stay, "caught": len(c["caught"]), "missed": len(c["missed"]),
                        "expected": expect, "fraction_model": min(1.0, max(0.0, (stay - GRAB[0]) / (GRAB[1] - GRAB[0])))}
        runs[tag] = res["series"]
    checks.append(check("trend", "📉", all(a >= b for a, b in zip(caught_counts, caught_counts[1:]))
                        and caught_counts[0] > caught_counts[-1], value=caught_counts,
                        note="流れが速いほど取れる数が減る"))
    # ✋ 手が空くまで（口へ運ぶ間）は取れない: 1 秒おきに来る麺を、取ったら 1.5 秒休む
    res = catch_run(flows[0], noodles=noodles, cooldown=1.5, seed=seed)
    c = res["catchers"][0]
    busy = [m for m in c["missed"] if m["why"] == "busy"]
    gaps = [b["t_catch"] - a["t_catch"] for a, b in zip(c["caught"], c["caught"][1:])]
    checks += [
        check("busy:missed", "✋", len(busy) > 0 and all(m["free_at"] > m["t_in"] + m["t_grab"] for m in busy),
              value=len(busy), note="手がふさがっている間に来た麺は取り逃す"),
        check("busy:spacing", "⏱️", all(g >= 1.5 - 1e-9 for g in gaps), value=min(gaps, default=None), expected=1.5,
              note="取る間隔は休む時間以上"),
        check("busy:fewer", "🥢", len(c["caught"]) < caught_counts[0], value=len(c["caught"]),
              expected=caught_counts[0], note="休む時間があると取れる数が減る"),
        conserved("busy:conserve", res),
    ]
    summary["busy"] = {"caught": len(c["caught"]), "busy_missed": len(busy),
                       "fast_missed": len(c["missed"]) - len(busy)}
    return {"params": {"flows": list(flows), "noodles": noodles, "seed": seed, "grab": list(GRAB),
                       "zone": CATCH_ZONE, "zone_length": dx},
            "summary": summary, "runs": runs, "checks": checks}


# --- E15 冷たい水（🧊 水温で粘度・密度が変わる） ---------------------------------

def _path_time(tr: dict, node: str) -> float | None:
    return next((p["t"] for p in tr["path"] if p["node"] == node), None)


def e15_cold_water(cold: float = 10.0, warm: float = 20.0, q_trickle: float = 3e-6, q_flow: float = 5e-5) -> dict:
    """流しそうめんの水は冷たい（10 °C 前後）。温度表（TEMP_TABLES）で μ・ρ を動かし、
    細い流れ（層流）が深く・遅くなること、縁の越流水深と引っかかりは変わらないこと（越流の式は粘度によらない）、
    Manning の流れはビット単位で同じことを確かめる。"""
    sec = Section(BAMBOO["shape"], BAMBOO["size"])
    rim = sec.size
    checks, runs, summary = [], {}, {}

    # 📚 温度表: 別の出典との幅、ν の比
    for t, ref in TEMP_TABLES["water"]["cross"].items():
        checks.append(close(f"table:water:{t:g}C", "📚", get_fluid("water", t).mu * 1e3, ref, 3e-3,
                            note="Anton Paar（IAPWS 2008）と chem-casts の幅"))
    w_c, w_w = get_fluid("water", cold), get_fluid("water", warm)
    (m_c, r_c), (m_w, r_w) = table_value("water", cold), table_value("water", warm)
    checks.append(close("table:nu-ratio", "🧊", w_c.nu / w_w.nu, (m_c / r_c) / (m_w / r_w), 1e-4,
                        note=f"ν({cold:g} °C)/ν({warm:g} °C) = 表の比（20 °C を FLUIDS の値にそろえた分のずれまで）"))
    checks.append(check("placeholder", "🏷️", get_fluid("namekian-water", cold).nu == w_c.nu,
                        value=get_fluid("namekian-water", cold).nu, expected=w_c.nu,
                        note="仮置きの流体は水と同じ温度表で動く"))
    summary["nu"] = {"cold": w_c.nu, "warm": w_w.nu, "ratio": w_c.nu / w_w.nu}
    summary["ethanol_spread"] = {f"{t:g}C": {"ours": get_fluid("ethanol", t).mu * 1e3, "kaye_laby": v}
                                 for t, v in TEMP_TABLES["ethanol"]["cross"].items()}

    # 🥢 細い流れ（層流）を 20 °C と 10 °C で
    got = {}
    for tag, temp in (("warm", warm), ("cold", cold)):
        fl = get_fluid("water", temp)
        net = somen_layout(q_trickle, law="laminar", beta=3.0, dt=0.05, temp_c=temp)
        net.tracer("hiyamugi", "take/0:0", release=40.0, stop="basin", size=NOODLES["hiyamugi"])
        net.tracer("somen", "take/0:0", release=45.0, stop="basin", size=NOODLES["somen"])
        res = net.run(110.0)
        hn = normal_depth(sec, q_trickle, 0.05, "laminar", nu=fl.nu, h_max=rim, beta=3.0)
        vn = q_trickle / sec.area(hn)
        tr = {t["name"]: t for t in res["tracers"]}
        t0, t19 = _path_time(tr["somen"], "take/0:0"), _path_time(tr["somen"], "take/0:19")
        took = None if t19 is None else t19 - t0
        want = 19 * 0.1 / vn  # 1 本目の竹の 19 区間を抜ける時間
        re = channel_reynolds(sec, hn, q_trickle, fl.nu)
        got[tag] = {"depth": net.head("take/0:19"), "hn": hn, "took": took, "want": want,
                    "pond": net.head("take/1:9"), "tracers": tr, "re": re}
        checks += [
            close(f"{tag}:depth", "📏", got[tag]["depth"], hn, 2e-3, note=f"{temp:g} °C の等流水深（層流）"),
            check(f"{tag}:pass", "⏱️", took is not None and rel_err(took, want) <= 2 * 0.05 / want + 1e-3,
                  value=took, expected=want, tol=2 * 0.05 / want + 1e-3, note=f"{temp:g} °C でそうめんが 19 区間を抜ける時間"),
            check(f"{tag}:regime", "〰️", re < 500, value=re, expected=500.0, note="層流（Re < 500）"),
            conserved(f"{tag}:conserve", res),
        ]
        runs[f"trickle:{tag}"] = res["series"]
    c, w = got["cold"], got["warm"]
    checks += [
        close("cold:deeper", "🧊📏", c["depth"] / w["depth"], c["hn"] / w["hn"], 1e-3,
              note="冷たいと深い（等流水深の比どおり。広い流れなら (ν 比)^(1/3)）"),
        close("cold:slower", "🧊⏱️", c["took"] / w["took"], c["want"] / w["want"], 2e-2,
              note="冷たいと遅い（流速の逆比どおり）"),
        check("cold:same-head", "🫙", abs(c["pond"] - w["pond"]) <= 1e-9, value=c["pond"] - w["pond"], expected=0.0,
              tol=1e-9, note="縁の越流水深は水温で変わらない（越流の式に粘度が入らない）"),
        check("cold:same-snag", "🪝", all(g["tracers"]["hiyamugi"]["stuck_at"] == "take/1:9" and
                                         g["tracers"]["somen"]["state"] == "done" for g in (c, w)),
              note="ひやむぎは両方で縁に引っかかり、そうめんは両方で通る"),
    ]
    summary["trickle"] = {k: {"depth": v["depth"], "normal_depth": v["hn"], "pass": v["took"], "pond": v["pond"],
                              "head_over_lip": v["pond"] - 3e-3, "reynolds": v["re"]} for k, v in got.items()}

    # 🎋 Manning（乱流・粗さ支配）は水温によらない → 各セルの体積がビット単位で同じ
    a = somen_layout(q_flow).run(30.0)["volumes_q"]
    b = somen_layout(q_flow, temp_c=cold).run(30.0)["volumes_q"]
    checks.append(check("manning:same", "🟰", a == b, note="Manning の竹は 10 °C でも 20 °C と同じ（ビット単位）"))
    return {"params": {"cold": cold, "warm": warm, "q_trickle": q_trickle, "q_flow": q_flow,
                       "sources": list(TEMP_TABLES["water"]["sources"])},
            "summary": summary, "runs": runs, "checks": checks}


# --- E16 竹の傷み（🎋 割れ目からの漏れ・節の削り残し） ----------------------------

WEAR_NODE = {"at": 0.5, "height": 2e-3}  # 竹の途中に残った節（縁 2 mm）
WEAR_CRACK = {"cell": 15, "area": 1e-6, "cd": 0.61}  # 割れ目（1 mm²、薄い縁の穴の Cd）


def bamboo_layout(q: float, *, law: str = "manning", beta: float | None = None, dt: float = 0.01,
                  node: dict | None = None, crack: dict | None = None, temp_c: float | None = None,
                  fluid: str = "water") -> Network:
    """桶 → ポンプ → 竹 1 本（途中に節の削り残し、下流に割れ目）→ 受け皿。漏れは床（floor）へ。"""
    node = WEAR_NODE if node is None else node
    crack = WEAR_CRACK if crack is None else crack
    net = Network(dt, fluid, temp_c=temp_c).tank("bucket", 0.05, volume=0.5).tank("basin", 0.05).tank("floor", 1.0)
    net.flume("take", "basin", [{"length": 2.0, "cells": 20, "slope": 0.05, "section": BAMBOO, "law": law,
                                 "beta": beta, "node_remnant": node, "cracks": [crack]}], spill="floor")
    net.pump("bucket", "take/0:0", rate=q)
    return net


def e16_bamboo_wear(q_flow: float = 5e-5, q_trickle: float = 3e-6) -> dict:
    """割れ目は Torricelli（Cd·a·√(2gh)、h はそのセルの水深）で床へ漏れ、節の削り残しは竹の途中の低い縁として
    手前に水をためる。漏れた分は床へちょうど届き、細い流れでは太い麺（うどん）が節で引っかかる。"""
    sec = Section(BAMBOO["shape"], BAMBOO["size"])
    cut = round(20 * WEAR_NODE["at"])
    node_cell = f"take/0:{cut - 1}"
    crack_cell = f"take/0b:{WEAR_CRACK['cell'] - cut}"
    checks, runs, summary = [], {}, {}

    net = bamboo_layout(q_flow)
    res = net.run(30.0)
    crack_ln = next(ln for ln in res["links"] if ln["kind"] == "orifice")
    h_c = net.head(crack_cell)
    leak = res["final"][f"q:{crack_ln['index']}"]
    leak_want = WEAR_CRACK["cd"] * WEAR_CRACK["area"] * math.sqrt(2 * 9.80665 * h_c)
    to_floor = [ln for ln in res["links"] if ln["dst"] == "floor"]
    checks += [
        close("flow:node-pond", "🎋", net.head(node_cell), lip_head(sec, q_flow, WEAR_NODE["height"]), 2e-3,
              note="節の手前にたまる水深 = 節の高さ + 越流水深"),
        close("flow:leak", "💧", leak, leak_want, 2e-3, note="割れ目の漏れ = Cd·a·√(2gh)（h は割れ目のセルの水深）"),
        close("flow:balance", "🔁", _into(res, "basin") + leak, q_flow, 1e-3, note="受け皿 + 漏れ = 流した量"),
        check("flow:floor-exact", "⚖️", res["volumes_q"]["floor"] == crack_ln["volume_moved_q"] and
              all(ln["volume_moved_q"] == 0 for ln in to_floor if ln is not crack_ln and ln["kind"] != "orifice"),
              value=res["volumes_q"]["floor"], expected=crack_ln["volume_moved_q"], tol=0,
              note="床に落ちたのは割れ目から漏れた分だけ（ulp まで）"),
        conserved("flow:conserve", res),
    ]
    runs["flow"] = res["series"]
    summary["flow"] = {"node_pond": net.head(node_cell), "crack_depth": h_c, "leak": leak, "leak_closed": leak_want,
                       "leak_cd_0.82": 0.82 / WEAR_CRACK["cd"] * leak_want}

    net = bamboo_layout(q_trickle, law="laminar", beta=3.0, dt=0.05)
    for i, name in enumerate(("udon", "hiyamugi", "somen")):
        net.tracer(name, "take/0:0", release=40.0 + i, stop="basin", size=NOODLES[name])
    res = net.run(100.0)
    head = lip_head(sec, q_trickle, WEAR_NODE["height"]) - WEAR_NODE["height"]
    tr = {t["name"]: t for t in res["tracers"]}
    checks += [
        check("trickle:premise", "📐", NOODLES["hiyamugi"] <= head < NOODLES["udon"], value=head,
              note="節の越流水深がひやむぎ以上・うどん未満になる流し方"),
        check("trickle:udon-snag", "🪝", tr["udon"]["stuck_at"] == node_cell and tr["udon"]["state"] != "done",
              value=tr["udon"]["stuck_time"], note="うどんは節で引っかかる"),
        check("trickle:pass", "🍜", tr["hiyamugi"]["state"] == "done" and tr["somen"]["state"] == "done",
              value=[tr["hiyamugi"]["arrival"], tr["somen"]["arrival"]], note="ひやむぎ・そうめんは節を越える"),
        conserved("trickle:conserve", res),
    ]
    runs["trickle"] = res["series"]
    summary["trickle"] = {"head_over_node": head, "tracers": {k: {"arrival": v["arrival"], "stuck_time": v["stuck_time"],
                                                                  "stuck_at": v["stuck_at"]} for k, v in tr.items()}}
    return {"params": {"q_flow": q_flow, "q_trickle": q_trickle, "node": WEAR_NODE, "crack": WEAR_CRACK},
            "summary": summary, "runs": runs, "checks": checks}


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
    "E12": ("noodle-length", "🍜", "長さのある麺（体が区間をまたぐ・尾で引っかかる・曲がりで飛び出す）", e12_noodle_length),
    "E13": ("noodle-clog", "🧱", "麺のたまり（縁に積もってせき止め、一気に押し流される）", e13_noodle_clog),
    "E14": ("chopsticks", "🍽️", "箸ですくう（遅い麺ほど取れる、手が空くまでは取れない）", e14_chopsticks),
    "E15": ("cold-water", "🧊", "冷たい水（水温で粘度・密度が変わる）", e15_cold_water),
    "E16": ("bamboo-wear", "🪵", "竹の傷み（割れ目からの漏れ・節の削り残し）", e16_bamboo_wear),
}
SLUGS = {slug: eid for eid, (slug, *_rest) in EXPERIMENTS.items()}


def _plain(token: str) -> str:
    """絵文字の比較用: 前後の空白と異体字セレクタ（U+FE0F）を落とす（🏛️ と 🏛 を同じに扱う）。"""
    return token.strip().replace("\ufe0f", "")


def resolve(key: str) -> str:
    """'E02'・'e02'・'fluid-drain'・'E02:fluid-drain'・'💧'（実験のアイコン）のどれでも ID に直す。"""
    k = _plain(key)
    if ":" in k:
        k = k.split(":", 1)[0]
    if k.upper() in EXPERIMENTS:
        return k.upper()
    if k in SLUGS:
        return SLUGS[k]
    for eid, (_slug, icon, _t, _f) in EXPERIMENTS.items():
        if _plain(icon) == k:
            return eid
    raise KeyError(f"unknown experiment {key!r}")


# --- 絵文字をそのままコマンドに -------------------------------------------------

VENDOR = Path(__file__).resolve().parent / "vendor"


def vendored(name: str):
    """vendor/<name>.py（vendor.lock.json で GitHub の commit に固定したファイル）を読み込む。"""
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, VENDOR / f"{name}.py")
    if spec is None or spec.loader is None:
        raise ImportError(f"vendor/{name}.py not found")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def fluid_aliases() -> dict[str, str]:
    """流体のキー・アイコン（💧 🍶 …）・ブロック（🟦 🟨 …）→ キー。"""
    out = {}
    for key, f in FLUIDS.items():
        for token in (key, f.icon, pw.EMOJI.get(key, "")):
            if token:
                out[_plain(token)] = key
    return out


SCENE_ICONS = {"🏛": "moniac-rig", "🥢": "somen-flume", "🎋": "flume", "🫙": "flume-hopper", "🌊": "dam-break",
               "🕳": "moniac-drain", "🧊": "cold-water", "🪵": "bamboo-wear", "🍜": "noodle-length",
               "🧱": "noodle-clog", "🍽": "chopsticks"}


def fluid_arg(text: str) -> str:
    """--fluid: キーでも 💧 でも 🟦 でもよい。"""
    return fluid_aliases().get(_plain(text), text)


def scene_arg(text: str) -> str:
    """draw の場面: 名前でも 🏛️ 🥢 🎋 🫙 🌊 🕳️ 🧊 🪵 🍜 🧱 🍽️ でもよい。"""
    return SCENE_ICONS.get(_plain(text), text)


def experiment_arg(text: str) -> str:
    try:
        return resolve(text)
    except KeyError:
        return text


def block_arg(text: str) -> tuple[str, str]:
    """draw の 2 つ目: 流体のブロック・アイコンなら流体、それ以外の絵文字ならそのまま塗る色。"""
    key = fluid_aliases().get(_plain(text))
    if key:
        return ("fluid", key)
    if text.strip() and not any(c.isascii() and c.isalnum() for c in text):
        return ("fill", text.strip())
    raise argparse.ArgumentTypeError(f"{text!r}: 流体（{', '.join(sorted(FLUIDS))}）か絵文字のブロック")


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


def rig_section(scene: str, fluid: str = "water", temp_c: float | None = None) -> tuple[dict, list[str]]:
    """断面図で見る試験台（run の結果, 列）。moniac-rig は頭・貯める・戻りの管・横の槽（下の槽は水源なので省く）、
    somen-flume は竹 3 本と受け皿・床、cold-water は 10 °C の細い流れ、bamboo-wear は節と割れ目のある竹。
    temp_c は水温（somen-flume・cold-water・bamboo-wear。cold-water は省くと 10 °C）。"""
    if scene == "moniac-rig":
        return moniac_rig(fluid=fluid).run(1200.0), ["header", "store", "store>sump", "side"]
    if scene == "somen-flume":
        net = somen_layout(5e-5, fluid=fluid, temp_c=temp_c)
        for i in range(3):
            net.tracer(f"somen-{i + 1}", "take/0:0", release=20.0 + 2 * i, stop="basin", size=NOODLES["somen"])
        return net.run(40.0, sample_every=20), ["take/0", "take/1", "take/2", "basin", "floor"]
    if scene == "noodle-length":
        net = somen_layout(5e-5, fluid=fluid, temp_c=temp_c)
        net.tracer("long", "take/0:0", release=20.0, stop="basin", length=0.3, **STRAND)
        return _from(net.run(32.0, sample_every=10), 20.0), ["take/0", "take/1", "take/2"]
    if scene == "noodle-clog":
        net = somen_layout(3e-6, law="laminar", beta=3.0, dt=0.05, fluid=fluid, temp_c=temp_c)
        for i in range(16):
            net.tracer(f"bunch-{i + 1:02d}", "take/0:0", release=80.0 + 5.0 * i, stop="basin", length=0.2,
                       strands=10, **STRAND)
        return _from(net.run(200.0, sample_every=20), 100.0), ["take/1", "take/2", "basin"]
    if scene == "chopsticks":
        return _from(catch_run(5e-5, noodles=10), 24.0), ["take/1", "take/2", "basin"]
    if scene == "cold-water":
        net = somen_layout(3e-6, law="laminar", beta=3.0, dt=0.05, fluid=fluid, temp_c=10.0 if temp_c is None else temp_c)
        net.tracer("hiyamugi", "take/0:0", release=40.0, stop="basin", size=NOODLES["hiyamugi"])
        net.tracer("somen", "take/0:0", release=45.0, stop="basin", size=NOODLES["somen"])
        return net.run(110.0, sample_every=10), ["take/0", "take/1", "take/2", "basin"]
    if scene == "bamboo-wear":
        net = bamboo_layout(3e-6, law="laminar", beta=3.0, dt=0.05, fluid=fluid, temp_c=temp_c)
        for i, name in enumerate(("udon", "somen")):
            net.tracer(name, "take/0:0", release=40.0 + i, stop="basin", size=NOODLES[name])
        return net.run(100.0, sample_every=10), ["take/0", "take/0b", "basin", "floor"]
    raise KeyError(f"unknown rig scene {scene!r}")


def _from(res: dict, t0: float) -> dict:
    """断面図に出す時刻を t0 から後に絞る（水が行き渡って麺を流し始めてから）。"""
    return {**res, "series": [row for row in res["series"] if row["t"] >= t0]}


RIG_SCENES = ("moniac-rig", "somen-flume", "noodle-length", "noodle-clog", "chopsticks", "cold-water",
              "bamboo-wear")


def _draw(args) -> int:
    if args.scene in RIG_SCENES:
        if args.view == "files":
            build_parser().error(f"{args.scene} は --view frames だけ")
        res, cols = rig_section(args.scene, args.fluid, args.temp_c)
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
    """引数は argparse の choices / type で宣言する（vendor/cli_args.py の make_parser。省略形は受け付けない）。
    実験・場面・流体は名前のほか絵文字をそのまま書ける: ``run 🥢``、``draw 🏛️ 🟩``、``draw 🥢 🟥``。"""
    cli = vendored("cli_args")
    p = cli.make_parser(suggest_unknown=True, prog="fluidfullspeedsimulation", description="🌊 流れ・遅れ・検証")
    sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("list", help="🗂️ 実験の一覧", allow_abbrev=False)
    sub.add_parser("fluids", help="🧪 流体の物性", allow_abbrev=False)
    run = sub.add_parser("run", help="▶️ 実験を回す", allow_abbrev=False)
    run.add_argument("experiment", type=experiment_arg, choices=sorted(EXPERIMENTS),
                     help="ID・識別子・アイコン（例 E10 / somen-flume / 🥢）")
    ver = sub.add_parser("verify", help="✅ 全実験を回して検証", allow_abbrev=False)
    drw = sub.add_parser("draw", help="🎨 🟦 を積む断面図で溜まる・流れる・遅れるを見る（ドット絵・GIF・PNG も）",
                         allow_abbrev=False)
    scenes = RIG_SCENES + pw.SECTION_CASES + pw.SCENES
    drw.add_argument("scene", type=scene_arg, choices=scenes,
                     help="場面（名前か 🏛️ 🥢 🧊 🪵 🎋 🫙 🌊 🕳️）")
    drw.add_argument("block", nargs="?", type=block_arg,
                     help="流体のブロック・アイコン（🟦 🟨 🟧 🟩 🟪 / 💧 🍶 🧴 🟢 ☠️）か、塗りたい好きな絵文字")
    drw.add_argument("--view", choices=("frames", "frame", "files"), default="frames",
                     help="frames = 時刻ごとの断面図 / frame = 最後の 1 枚 / files = GIF・PNG（ドット絵だけ）")
    drw.add_argument("--fluid", type=fluid_arg, choices=ALL_FLUID_KEYS, default=None,
                     help="作動流体（キーか絵文字。🏷️ 仮置きも選べる。既定 water）")
    drw.add_argument("--fill", help="水のブロック（色は自由。既定は流体ごと: 🟦 水 🟨 エタノール 🟧 メタノール "
                                    "🟩 Namekian Water 🟪 Poison Water）")
    drw.add_argument("--temp-c", type=float, help="🧊 水温 °C（🥢 somen-flume・🧊 cold-water・🪵 bamboo-wear。表の範囲内）")
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


def _draw_choices(parser: argparse.ArgumentParser, args) -> None:
    """2 つ目の絵文字（block）と --fluid / --fill をまとめる。"""
    kind, value = args.block if args.block else (None, None)
    if kind == "fluid":
        if args.fluid not in (None, value):
            parser.error(f"流体が 2 つ: {value} と --fluid {args.fluid}")
        args.fluid = value
    elif kind == "fill":
        if args.fill not in (None, value):
            parser.error(f"色が 2 つ: {value} と --fill {args.fill}")
        args.fill = value
    args.fluid = args.fluid or "water"


def _utf8_stdout() -> None:
    """絵文字を出すので、UTF-8 でない端末（Windows のパイプ等）では stdout を UTF-8 に切り替える。"""
    enc = (getattr(sys.stdout, "encoding", None) or "").lower().replace("-", "")
    if enc != "utf8" and hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "draw":
        _draw_choices(parser, args)
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
