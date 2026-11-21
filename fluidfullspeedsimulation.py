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
import hashlib
import importlib.util
import json
import math
import random
import sys
from copy import deepcopy
from pathlib import Path

from fixedq import ONE, from_q, lag_gain, lag_step, qdiv, qmul, random_scales, shift_lag_step, to_q
from fluids import FLUIDS, PLACEHOLDER_MU_SCALE, REAL_FLUIDS, RE_LAMINAR, Section, channel_reynolds, fall_time
from fluids import bisection, fluid_table, get_fluid, lip_head
from fluids import ORIFICE_DELTA, TEMP_TABLES, normal_depth, superelevation, table_value, viscous_orifice
from fluids import overland_equilibrium, overland_power_law, overland_recession, overland_rise, overland_t_e_power_law
from fluids import G, poiseuille_head, poiseuille_tau, torricelli_empty_time, torricelli_head, transport_delay
from fluids import mariotte_constant_time, mariotte_flow, siphon_drain_time, siphon_drain_time_with_inflow, siphon_flow
from fluids import tip_period
from fluids import AZEOTROPE_X, abv, bubble_t, mix_volume, mole_fraction, moles_for_abv, rayleigh_log, vle_y
from network import OUTLET_MODES, Network, critical_gain, inverse_table, jsonable, peak_to_peak, register_law, settling_time, window
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


def sweep_points(seed: int, span: tuple) -> list[float]:
    """🎚️ 仮定を振る点: 出典の幅 (下限, 使う値, 上限) の両端と、内側の 1 点（A×10^-1 の A を乱数で）。"""
    lo, _use, hi = span
    return [lo, lo + random_scales(seed, 1, n_min=1, n_max=1)[0]["value"] * (hi - lo), hi]


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
                 fluid: str = "water", temp_c: float | None = None, tail_cells: int = 10, mu: float | None = None,
                 mu_scale: float | None = None) -> Network:
    """桶 → ポンプ → 竹 3 本（直線 → 段差 → 曲がり → 縁の継ぎ目 → 直線）→ 受け皿 → ホースで桶へ戻す。
    曲がりや継ぎ目からあふれた水は床（floor）へ。temp_c で水温を変えられる（🧊 E15）。"""
    net = Network(dt, fluid, temp_c=temp_c, mu=mu, mu_scale=mu_scale).tank("bucket", 0.05, volume=bucket).tank("basin", 0.05)
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


def rig_steady(q: float) -> dict:
    """試験台の定常（閉じた式）。頭の槽 h = Q·τ/(A·Σo)、分流は開き具合の比、
    貯める槽は cam(h)·A·h/τ = Q·o_store/Σo を二分法、横の槽は V = q_side·τ / cam_side(h_store)。"""
    total = sum(RIG_OPEN.values())
    flows = {k: q * o / total for k, o in RIG_OPEN.items()}
    store = bisection(lambda h: interp(RIG_STORE_CAM, h) * 0.02 * h / 20.0 - flows["store"] < 0, 0.0, 10.0)
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
    return bisection(lambda q: bend_outer(q) - target < 0, lo, hi)


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
    nk = PLACEHOLDER_MU_SCALE["namekian-water"]
    checks.append(close("placeholder", "🏷️", get_fluid("namekian-water", cold).nu, w_c.nu * nk, 1e-12,
                        note=f"仮置きの流体は水と同じ温度表で動く（粘度は水の {nk:g} 倍のまま）"))
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
                  fluid: str = "water", mu: float | None = None, mu_scale: float | None = None) -> Network:
    """桶 → ポンプ → 竹 1 本（途中に節の削り残し、下流に割れ目）→ 受け皿。漏れは床（floor）へ。"""
    node = WEAR_NODE if node is None else node
    crack = WEAR_CRACK if crack is None else crack
    net = Network(dt, fluid, temp_c=temp_c, mu=mu, mu_scale=mu_scale).tank("bucket", 0.05, volume=0.5).tank(
        "basin", 0.05).tank("floor", 1.0)
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


# --- E17 粘度のつまみ（🍯） ------------------------------------------------------

VISC_SCALES = (0.5, 2.0, 8.0)  # 水の何倍か（さらさら〜どろっと）


def e17_viscosity(q_trickle: float = 3e-6, scales=VISC_SCALES, drain_scales=(1.0, 30.0, 1000.0)) -> dict:
    """粘度を直接変える（mu_scale）。効くところ（層流の樋・細管）だけが動き、効かないところ（Manning・小孔・堰・縁）
    はビット単位で同じ。粘度が上がって小孔の Re が下がったら、viscous=True の Cd(Re) で抜け方が遅くなる。"""
    checks, runs, summary = [], {}, {}
    sec = Section(BAMBOO["shape"], BAMBOO["size"])
    w = get_fluid("water")
    base = None
    for k in (1.0, *scales):
        fl = get_fluid("water", mu_scale=k)
        net = somen_layout(q_trickle, law="laminar", beta=3.0, dt=0.05, mu_scale=k)
        net.tracer("somen", "take/0:0", release=60.0, stop="basin", size=NOODLES["somen"])
        res = net.run(160.0)
        hn = normal_depth(sec, q_trickle, 0.05, "laminar", nu=fl.nu, h_max=sec.size, beta=3.0)
        tr = res["tracers"][0]
        took = None if tr["arrival"] is None else tr["arrival"] - tr["path"][0]["t"]
        tag = f"x{k:g}"
        summary[tag] = {"nu": fl.nu, "depth": net.head("take/0:19"), "depth_exact": hn, "travel": took,
                        "pond": net.head("take/1:9")}
        checks += [
            close(f"{tag}:depth", "📏", net.head("take/0:19"), hn, 2e-3, note="層流の等流水深（ν に効く）"),
            conserved(f"{tag}:conserve", res),
        ]
        if base is None:
            base = summary[tag]
        else:
            checks.append(check(f"{tag}:pond-same", "🫙", abs(summary[tag]["pond"] - base["pond"]) <= 1e-12,
                                value=summary[tag]["pond"], expected=base["pond"],
                                note="縁を越える堰の式に粘度は入らない → 継ぎ目のたまりは同じ"))
            if took is not None and base["travel"]:
                checks.append(check(f"{tag}:slower", "🐌", (took > base["travel"]) == (k > 1.0), value=took,
                                    expected=base["travel"], note="粘いほどそうめんが遅く着く（さらさらなら速い）"))
    # Manning（乱流）の竹は粘度によらない
    vols = []
    for k in (1.0, 8.0):
        res = somen_layout(5e-5, mu_scale=k).run(20.0)
        vols.append([res["volumes_q"][f"take/0:{i}"] for i in range(20)])
    checks.append(check("manning-free", "🟰", vols[0] == vols[1], note="Manning の竹は ×8 でも各区間の体積がビット単位で同じ"))
    # 小孔: Cd 一定なら粘度によらない。viscous=True では Re が下がると Cd が下がって遅くなる
    a_orf = math.pi * 0.002 ** 2
    t_drain = {}
    for k in drain_scales:
        for viscous in (False, True):
            net = _sea(Network(0.1, mu_scale=k)).tank("tank", 0.01, head=0.3)
            net.orifice("tank", "sea", area=a_orf, viscous=viscous)
            res = net.run(200.0)
            ln = res["links"][0]
            tag = f"orifice:x{k:g}:{'visc' if viscous else 'flat'}"
            t_drain[tag] = res["volumes_q"]["tank"]
            summary[tag] = {"left": from_q(res["volumes_q"]["tank"]), "re_max": ln["reynolds_max"],
                            "cd": ln["cd_at_max_head"], "viscous_regime": ln["viscous_regime"]}
            checks.append(conserved(f"{tag}:conserve", res))
    flat = [t_drain[f"orifice:x{k:g}:flat"] for k in drain_scales]
    visc = [t_drain[f"orifice:x{k:g}:visc"] for k in drain_scales]
    checks += [
        check("orifice:flat-free", "🕳️", len(set(flat)) == 1, note="Cd 一定の小孔は粘度によらない（残りの体積が同じ）"),
        check("orifice:visc-slower", "🍯", visc == sorted(visc) and visc[0] < visc[-1], value=[from_q(v) for v in visc],
              note="viscous=True: 粘いほど 200 s 後に多く残る"),
        close("orifice:re-scales", "🚩", summary[f"orifice:x{drain_scales[-1]:g}:flat"]["re_max"] * drain_scales[-1],
              summary[f"orifice:x{drain_scales[0]:g}:flat"]["re_max"] * drain_scales[0], 1e-9,
              note="Cd 一定なら Re ∝ 1/ν。Re < 1e4 で viscous_regime の印（水でもこの 4 mm の穴は 5.9e3 と境目より下）"),
        close("orifice:water-unchanged", "💧", summary[f"orifice:x{drain_scales[0]:g}:visc"]["cd"], 0.61, 1e-3,
              note="水の粘度なら Cd(Re) の補正は効かない（Cd ≈ 0.61 のまま）"),
    ]
    cd, re = viscous_orifice(0.3, math.sqrt(4 * a_orf / math.pi), w.nu * 1e5)
    checks.append(close("orifice:laminar-limit", "📐", cd, ORIFICE_DELTA * math.sqrt(re), 1e-2,
                        note="Re ≪ 1（粘度 ×1e5）では Cd ≈ δ·√Re（δ = 0.2、Merritt の層流の傾き）"))
    # 仮置きの流体: Poison Water は水の 5 倍の粘度
    fl = get_fluid("poison-water")
    net = somen_layout(q_trickle, law="laminar", beta=3.0, dt=0.05, fluid="poison-water")
    res = net.run(100.0)
    hn = normal_depth(sec, q_trickle, 0.05, "laminar", nu=fl.nu, h_max=sec.size, beta=3.0)
    checks.append(close("poison-water:depth", "☠️", net.head("take/0:19"), hn, 2e-3,
                        note=f"☠️ Poison Water（粘度 ×{PLACEHOLDER_MU_SCALE['poison-water']:g}）の層流の水深"))
    return {"params": {"q_trickle": q_trickle, "scales": list(scales), "drain_scales": list(drain_scales)},
            "summary": summary, "runs": runs, "checks": checks}


# --- E19〜E21 溜まる → 閾値 → 一気に（鹿威し・ピタゴラスの杯・マリオットびん） --------------------

TUBE = {"area": 2e-3, "v_tip": 1e-4}  # 🎍 筒（断面 20 cm²、100 mL で傾く）
CUP = {"area": 3e-3, "h_crest": 0.08, "h_drop": 0.05, "tube": 2e-5, "cd": 0.6}  # 🍷 杯（30 cm²、頂 8 cm、出口は底より 5 cm 下）
BOTTLE = {"area": 5e-3, "h0": 0.30, "h_air": 0.10, "orifice": 1e-5, "cd": 0.61}  # 🫗 びん（50 cm²、空気管の下端 10 cm）


def shishi_layout(q: float, *, dt: float = 0.05, t_dump: float | None = None, dump_rate: float | None = None,
                  fluid: str = "water") -> Network:
    """桶 → ポンプ（細い注ぎ込み q）→ 🎍 筒 → 傾いて受け皿へ。"""
    net = Network(dt, fluid).tank("bucket", 0.05, volume=0.5).tank("tube", TUBE["area"]).tank("basin", 0.01)
    net.pump("bucket", "tube", rate=q)
    net.tip("tube", "basin", v_tip=TUBE["v_tip"], t_dump=t_dump, dump_rate=dump_rate)
    return net


def _intervals(times: list[float]) -> list[float]:
    return [b - a for a, b in zip(times, times[1:])]


def e19_shishi_odoshi(q: float = 2e-6, dt: float = 0.05, cycles: float = 4.2, t_dump: float = 1.5,
                      dump_rate: float = 1e-4, seed: int = 19) -> dict:
    """🎍 注ぎ込み q が筒に溜まり、V_tip で傾いて一気に落ちる。周期 T = V_tip/q + t_dump（傾いている間の注ぎ込みも
    一緒に落ちる）。落ちた量は「傾いた瞬間の体積 + 傾いている間に注がれた分」と整数で一致し、音（clack）の間隔は毎回同じ。"""
    checks, runs, summary = [], {}, {}
    v_tip = TUBE["v_tip"]

    def run_case(tag: str, q: float, t_dump_: float | None, rate: float | None, icon: str, dt: float = dt) -> list[dict]:
        net = shishi_layout(q, dt=dt, t_dump=t_dump_, dump_rate=rate)
        amount = qmul(to_q(q), net.dt_q)
        link = net.links[1]
        t_end = cycles * tip_period(v_tip, q, link["t_dump"])
        res = net.run(t_end, trace=["tube"])
        tips = res["links"][1]["tips"]
        closed = [x for x in tips if x["dumped_q"] is not None]
        gaps = _intervals([x["t"] for x in tips])
        expected = tip_period(v_tip, q, link["t_dump"])
        summary[tag] = {"q": q, "t_dump": link["t_dump"], "clacks": [x["t"] for x in tips], "period": expected,
                        "mean_gap": sum(gaps) / len(gaps) if gaps else None, "dump_steps": link["dump_steps"]}
        runs[tag] = res["series"]
        out = [
            check(f"{tag}:count", icon + "🔔", len(tips) == int(cycles) and len(closed) >= int(cycles) - 1,
                  value=len(tips), expected=int(cycles), note="t_end = cycles·T の間に鳴る回数"),
            close(f"{tag}:period", icon + "⏱️", sum(gaps) / len(gaps), expected, 2 * dt / expected,
                  note="clack の間隔 = V_tip/q + t_dump（2 ステップまで）"),
            check(f"{tag}:even", icon + "🎵", max(gaps) - min(gaps) <= 1e-9, value=max(gaps) - min(gaps), tol=1e-9,
                  note="間隔は毎回同じ（戻った直後の残りが毎回 1 ステップ分の注ぎ込みで同じ）"),
            check(f"{tag}:threshold", icon + "📏",
                  all(to_q(v_tip) <= x["volume_q"] < to_q(v_tip) + amount for x in tips),
                  value=from_q(tips[0]["volume_q"]) if tips else None, expected=v_tip,
                  note="傾いた瞬間の体積は V_tip 以上、V_tip + 1 ステップ分未満"),
            check(f"{tag}:dump-exact", icon + "🪣",
                  all(x["dumped_q"] == x["volume_q"] + (link["dump_steps"] - 1) * amount for x in closed),
                  value=from_q(closed[0]["dumped_q"]) if closed else None,
                  note="落ちた量 = 傾いた瞬間の体積 + 傾いている間の注ぎ込み（0 ulp）"),
            check(f"{tag}:returns-empty", icon + "🔁", all(x["left_q"] == amount for x in closed),
                  value=from_q(closed[0]["left_q"]) if closed else None, expected=from_q(amount),
                  note="戻った筒には最後のステップの注ぎ込みだけが残る"),
            check(f"{tag}:basin", icon + "🫙", res["volumes_q"]["basin"] == res["links"][1]["volume_moved_q"]
                  and sum(x["dumped_q"] for x in closed) <= res["links"][1]["volume_moved_q"],
                  value=from_q(res["volumes_q"]["basin"]), note="受け皿に入ったのは筒から落ちた分だけ（0 ulp）"),
            conserved(f"{tag}:conserve", res),
        ]
        return out

    checks += run_case("instant", q, None, None, "⚡")
    checks += run_case("timed", q, t_dump, dump_rate, "🐢")
    need = v_tip + q * t_dump
    checks.append(check("timed:rate-enough", "🐢📐", dump_rate * t_dump >= need, value=dump_rate * t_dump, expected=need,
                        note="傾いている間に V_tip + q·t_dump を落とし切れる流量（前提）"))
    for item in random_scales(seed, 3, n_min=5, n_max=7):  # 1 周期を 1000 ステップに刻んで q を振る
        q_r = item["value"]
        checks += run_case(f"q={item['label']}", q_r, None, None, "🎲", dt=v_tip / q_r / 1000.0)[:5]
    summary["v_tip"] = v_tip
    return {"params": {"q": q, "dt": dt, "cycles": cycles, "t_dump": t_dump, "dump_rate": dump_rate, "tube": TUBE,
                       "seed": seed},
            "summary": summary, "runs": runs, "checks": checks}


def cup_layout(q_in: float, *, dt: float = 0.02, head: float = 0.0, fluid: str = "water") -> Network:
    """桶 → ポンプ → 🍷 杯 → サイフォン → 鉢。q_in = 0 なら注ぎ込みなし（head から排水だけ）。"""
    net = Network(dt, fluid).tank("bucket", 0.05, volume=0.5).tank("cup", CUP["area"], head=head).tank("bowl", 0.05)
    if q_in > 0:
        net.pump("bucket", "cup", rate=q_in)
    net.siphon("cup", "bowl", h_crest=CUP["h_crest"], h_drop=CUP["h_drop"], area=CUP["tube"], cd=CUP["cd"])
    return net


def e20_pythagorean_cup(q_in: float = 5e-6, dt: float = 0.02, dt_drain: float = 0.005, cycles: float = 3.1,
                        seed: int = 20) -> dict:
    """🍷 杯に注ぐと頂（h_crest）まで溜まり、サイフォンが働いて吸い口まで吸い切る。
    注ぎ込みなしの排水時間は閉じた式、注ぎ込みありの周期は「溜まる時間 + RK4 の参照値」と比べる。"""
    checks, runs, summary = [], {}, {}
    A, hc, hd, a, cd = CUP["area"], CUP["h_crest"], CUP["h_drop"], CUP["tube"], CUP["cd"]
    # 注ぎ込みなし: 頂まで入った杯が吸い切られる時間
    net = cup_layout(0.0, dt=dt_drain, head=hc)
    res = net.run(2.0 * siphon_drain_time(A, cd, a, hc, 0.0, hd), trace=["cup"])
    ln = res["links"][0]
    t_drain = siphon_drain_time(A, cd, a, hc, 0.0, hd)
    checks += [
        check("drain:primes-at-once", "🍷", ln["primes"] == [0.0], value=ln["primes"], note="頂まで入っていれば最初のステップで働く"),
        close("drain:time", "⏱️", ln["breaks"][0]["t"], t_drain, 2 * dt_drain / t_drain,
              note="吸い切る時間 = 2A/(Cd·a·√2g)·(√(h_crest + h_drop) − √h_drop)（2 ステップまで）"),
        check("drain:empty", "🫗", res["volumes_q"]["cup"] == 0 and ln["breaks"][0]["left_q"] == 0,
              value=res["volumes_q"]["cup"], expected=0, tol=0, note="吸い口（底）まで吸い切る: 残り 0 ulp"),
        check("drain:breaks-once", "🔂", len(ln["breaks"]) == 1 and not ln["primed"], value=len(ln["breaks"]),
              expected=1, note="空になったら止まり、注ぎ込みがなければ二度と働かない"),
        conserved("drain:conserve", res),
    ]
    runs["drain"] = res["series"]
    summary["drain"] = {"t_break": ln["breaks"][0]["t"], "t_closed_form": t_drain}
    # 注ぎ込みあり: 溜まる → 吸い切る → 溜まる…
    t_fill = A * hc / q_in
    t_ref = siphon_drain_time_with_inflow(A, cd, a, hc, 0.0, hd, q_in)
    period = t_fill + t_ref
    net = cup_layout(q_in, dt=dt)
    res = net.run(cycles * period, trace=["cup"])
    ln = res["links"][1]
    amount = qmul(to_q(q_in), net.dt_q)
    gaps = _intervals(ln["primes"])
    tr = res["trace"]["cup"]
    after_break = [tr[int(round(b["t"] / dt)) - 1] for b in ln["breaks"]]
    checks += [
        check("cycle:siphon-wins", "📐", siphon_flow(cd, a, hd) > q_in, value=siphon_flow(cd, a, hd), expected=q_in,
              note="吸い口での流量 > 注ぎ込み（吸い切れる前提）"),
        close("cycle:fill", "🫙", ln["primes"][0], t_fill, 2 * dt / t_fill, note="最初に頂へ届く時刻 = A·h_crest/q_in"),
        check("cycle:no-flow-before-prime", "🚫", ln["first_flow"] == ln["primes"][0], value=ln["first_flow"],
              expected=ln["primes"][0], note="頂に届くまでサイフォンは流れない"),
        check("cycle:rearms", "🔁", len(ln["primes"]) >= 3 and len(ln["breaks"]) >= 2, value=len(ln["primes"]),
              expected=3, note="吸い切ったあと、また頂まで溜まると働く"),
        close("cycle:period", "⏱️", sum(gaps) / len(gaps), period, 3 * dt / period,
              note="周期 = A·h_crest/q_in + 注ぎ込みありの排水時間（RK4 の参照値。3 ステップまで）"),
        check("cycle:empties", "🫗", all(b["left_q"] == 0 for b in ln["breaks"]) and
              all(h <= from_q(amount) / A + 1e-12 for h in after_break),
              value=max(after_break) if after_break else None, expected=from_q(amount) / A,
              note="止まった直後の杯には、そのステップの注ぎ込みしか残らない"),
        conserved("cycle:conserve", res),
    ]
    runs["cycle"] = res["series"]
    summary["cycle"] = {"primes": ln["primes"], "breaks": [b["t"] for b in ln["breaks"]], "t_fill": t_fill,
                        "t_drain_ref": t_ref, "period": period, "mean_gap": sum(gaps) / len(gaps)}
    # 注ぎ込みを乱数で振る: 溜まる時間
    for item in random_scales(seed, 3, n_min=5, n_max=6):
        qq = item["value"]
        net = cup_layout(qq, dt=dt)
        res = net.run(1.2 * A * hc / qq + 5.0)
        ln = res["links"][1]
        tag = f"q={item['label']}"
        checks += [
            close(f"{tag}:fill", "🎲🫙", ln["primes"][0], A * hc / qq, 2 * dt / (A * hc / qq),
                  note="頂へ届く時刻 = A·h_crest/q_in"),
            conserved(f"{tag}:conserve", res),
        ]
    return {"params": {"q_in": q_in, "dt": dt, "dt_drain": dt_drain, "cycles": cycles, "cup": CUP, "seed": seed},
            "summary": summary, "runs": runs, "checks": checks}


def bottle_layout(h_air: float, *, dt: float = 0.05, h0: float | None = None, fluid: str = "water") -> Network:
    """🫗 びん → マリオットの出口 → 受け皿。"""
    net = Network(dt, fluid).tank("bottle", BOTTLE["area"], head=BOTTLE["h0"] if h0 is None else h0).tank("basin", 0.05)
    net.mariotte("bottle", "basin", area=BOTTLE["orifice"], h_air=h_air, cd=BOTTLE["cd"])
    return net


def e21_mariotte(dt: float = 0.02, seed: int = 21) -> dict:
    """🫗 空気管の下端が水面より下にあるあいだ、流量は Q = Cd·a·√(2g·h_air) で一定（水位が下がっても変わらない）。
    水面が下端まで下がると、ふつうの Torricelli の排水に切り替わる。"""
    checks, runs, summary = [], {}, {}
    A, h0, h_air, a, cd = BOTTLE["area"], BOTTLE["h0"], BOTTLE["h_air"], BOTTLE["orifice"], BOTTLE["cd"]
    q_c = mariotte_flow(cd, a, h_air)
    t_c = mariotte_constant_time(h0, h_air, A, cd, a)
    t_e = torricelli_empty_time(h_air, A, a, cd)
    net = bottle_layout(h_air, dt=dt)
    link = net.links[0]
    res = net.run(t_c + 1.1 * t_e, trace=["bottle"])
    ln = res["links"][0]
    tr = res["trace"]["bottle"]
    n_c = int(round(ln["constant_time"] / dt))
    v0, out_q, area_q = to_q(h0 * A), link["const_out"], net.nodes["bottle"]["area_q"]
    linear = all(tr[i] == from_q(qdiv(v0 - (i + 1) * out_q, area_q)) for i in range(n_c))
    t_fall = ln["falling_from"]
    worst = 0.0
    for frac in (0.25, 0.5, 0.75):
        t = t_fall + frac * t_e
        worst = max(worst, abs(tr[int(round(t / dt)) - 1] - torricelli_head(h_air, t - t_fall, A, a, cd)))
    checks += [
        close("const:flow", "🫗", ln["constant_flow"], q_c, 4.0 / (q_c * dt * ONE),
              note="一定の流量 Q = Cd·a·√(2g·h_air)（1 ステップ分の Q 整数の丸め 4 ulp まで）"),
        check("const:linear", "📉", linear and n_c > 0, value=n_c, note="水位は毎ステップ同じ量だけ下がる（整数で厳密に直線）"),
        close("const:time", "⏱️", ln["constant_time"], t_c, 2 * dt / t_c, note="一定のまま流れる時間 = (h0 − h_air)·A/Q"),
        check("fall:torricelli", "🕳️", worst <= 1e-4, value=worst, tol=1e-4,
              note="下端より下がった後は Torricelli の排水曲線（前進 Euler、1/4・1/2・3/4 で 0.1 mm まで）"),
        close("fall:start", "🔀", t_fall, t_c, 2 * dt / t_c, note="切り替わる時刻 = 一定のまま流れる時間"),
        conserved("conserve", res),
    ]
    runs["bottle"] = res["series"]
    summary["bottle"] = {"q_const": ln["constant_flow"], "q_closed_form": q_c, "t_const": ln["constant_time"],
                         "t_closed_form": t_c, "falling_from": t_fall, "torricelli_worst_m": worst}
    for item in random_scales(seed, 3, n_min=2, n_max=2):  # 空気管の下端 1〜10 cm（h0 = 30 cm より下）
        ha = item["value"]
        qq = mariotte_flow(cd, a, ha)
        tc = mariotte_constant_time(h0, ha, A, cd, a)
        net = bottle_layout(ha, dt=0.1)
        res = net.run(1.05 * tc)
        ln = res["links"][0]
        tag = f"h_air={item['label']}"
        checks += [
            close(f"{tag}:flow", "🎲🫗", ln["constant_flow"], qq, 4.0 / (qq * 0.1 * ONE), note="一定の流量 = Cd·a·√(2g·h_air)"),
            close(f"{tag}:time", "🎲⏱️", ln["constant_time"], tc, 2 * 0.1 / tc, note="一定のまま流れる時間"),
            conserved(f"{tag}:conserve", res),
        ]
    return {"params": {"dt": dt, "bottle": BOTTLE, "seed": seed}, "summary": summary, "runs": runs, "checks": checks}


# --- E22 水時計（🏺 断面が高さで変わる器） --------------------------------------

CLOCK = {"h0": 0.30, "r_top": 0.10, "orifice_d": 2.0e-3, "cd": 0.61, "marks": 12}


def clepsydra_area(h0: float, r_top: float):
    """流出型の水時計の断面 A(h) = π r_top² √(h/h0)（半径 ∝ h^(1/4)、4 次の放物面）。
    Torricelli の流出 Cd·a·√(2gh) と釣り合って dh/dt が一定になる形。"""
    return lambda h: math.pi * r_top * r_top * math.sqrt(max(h, 0.0) / h0)


def cone_area(h0: float, r_top: float, r_bottom: float):
    """切頭円すい（カルナックの水時計のおおよその形）。半径が高さに比例するので A(h) ∝ r(h)²。"""
    return lambda h: math.pi * (r_bottom + (r_top - r_bottom) * max(h, 0.0) / h0) ** 2


def clock_rate(h0: float, r_top: float, orifice_d: float, cd: float, g: float = 9.80665) -> float:
    """水位の下がる速さ k = Cd·a·√(2g h0) / (π r_top²) [m/s]。空になる時間 T = h0 / k。"""
    a = math.pi * (orifice_d / 2) ** 2
    return cd * a * math.sqrt(2.0 * g * h0) / (math.pi * r_top * r_top)


def water_clock(shape: str = "clepsydra", *, dt: float | None = None, samples: int = 512, viscous: bool = False,
                fluid: str = "water", **over) -> Network:
    """水時計 1 つ（clock）と受け皿（sea）。shape: clepsydra（4 次の放物面）/ cylinder（同じ上の半径の筒）/ cone。
    dt を省くと空になる時間の 1/8000。"""
    c = {**CLOCK, **over}
    h0, r_top, a = c["h0"], c["r_top"], math.pi * (c["orifice_d"] / 2) ** 2
    k = clock_rate(h0, r_top, c["orifice_d"], c["cd"])
    net = Network(dt or (h0 / k) / 8000.0, fluid).tank("sea", 1.0)
    if shape == "clepsydra":
        net.tank("clock", h0, head=h0, profile=clepsydra_area(h0, r_top), samples=samples)
    elif shape == "cylinder":
        net.tank("clock", math.pi * r_top * r_top, head=h0)
    elif shape == "cone":
        net.tank("clock", h0, head=h0, profile=cone_area(h0, r_top, c.get("r_bottom", 0.5 * r_top)), samples=samples)
    else:
        raise ValueError("shape must be clepsydra, cylinder or cone")
    net.orifice("clock", "sea", area=a, cd=c["cd"], viscous=viscous)
    return net


def _marks(trace: list[float], dt: float, t_end: float, n: int, h0: float) -> list[float]:
    """時刻 T·j/n（j = 0..n）の水位（trace[i] は (i+1)·dt の値。間は線形補間、0 は h0）。"""
    pts = [h0] + list(trace)
    out = []
    for j in range(n + 1):
        x = (t_end * j / n) / dt
        i = min(int(x), len(pts) - 2)
        out.append(pts[i] + (x - i) * (pts[i + 1] - pts[i]))
    return out


def _linear_dev(trace: list[float], dt: float, h0: float, k: float, until: float) -> float:
    """h0 − k·t からの最大のずれ / h0（t ≤ until）。"""
    n = int(until / dt)
    return max(abs(h - (h0 - k * (i + 1) * dt)) for i, h in enumerate(trace[:n])) / h0


def e22_water_clock(dt_list=(1.0, 0.25, 0.05), samples_list=(64, 512, 4096), seed: int = 22) -> dict:
    """🏺 流出型の水時計。半径 ∝ h^(1/4) の器なら水位が一定の速さで下がり（目盛りが等間隔）、筒だと
    目盛りが詰まっていく。断面が高さで変わる水槽（profile）の精度を、閉じた式 h = h0 − k·t で測る。"""
    checks, runs, summary = [], {}, {}
    h0, r_top, n_marks = CLOCK["h0"], CLOCK["r_top"], CLOCK["marks"]
    k = clock_rate(h0, r_top, CLOCK["orifice_d"], CLOCK["cd"])
    t_empty = h0 / k
    a = math.pi * (CLOCK["orifice_d"] / 2) ** 2
    area_cyl = math.pi * r_top * r_top

    # 🏺 4 次の放物面: 線形に下がる・目盛りが等間隔・空になる時刻
    net = water_clock("clepsydra")
    res = net.run(t_empty * 1.02, trace=["clock"], sample_every=1)
    tr = res["trace"]["clock"]
    dev = _linear_dev(tr, net.dt, h0, k, 0.9 * t_empty)
    marks = _marks(tr, net.dt, t_empty, n_marks, h0)
    gaps = [a_ - b_ for a_, b_ in zip(marks, marks[1:])]
    t_dry = next(((i + 1) * net.dt for i, h in enumerate(tr) if h <= 0.0), None)
    checks += [
        check("clepsydra:linear", "📏", dev <= 2e-4, value=dev, expected=0.0, tol=2e-4,
              note="h = h0 − k·t からの最大のずれ / h0（t ≤ 0.9 T。前進 Euler の分だけ残る）"),
        check("clepsydra:marks-even", "🕰️", max(gaps) - min(gaps) <= 2e-3 * h0 / n_marks,
              value=[round(g_ / (h0 / n_marks), 5) for g_ in gaps], expected=1.0,
              note=f"{n_marks} 等分の時刻の水位の間隔が等しい（h0/{n_marks} の 0.2 % まで）"),
        check("clepsydra:empty-time", "⏱️", t_dry is not None and rel_err(t_dry, t_empty) <= 2e-3, value=t_dry,
              expected=t_empty, tol=2e-3, note="空になる時刻 T = h0/k = π r² h0 / (Cd a √(2g h0))"),
        conserved("clepsydra:conserve", res),
    ]
    runs["clepsydra"] = res["series"]
    summary["clepsydra"] = {"k": k, "t_empty": t_empty, "t_dry": t_dry, "linear_dev": dev, "marks": marks,
                            "dt": net.dt, "samples": res["profile_samples"]["clock"]}

    # 🥫 同じ上の半径の筒: Torricelli で h = (√h0 − c t)²。目盛りの間隔は 23:21:…:1 に詰まる
    c_cyl = CLOCK["cd"] * a * math.sqrt(2.0 * 9.80665) / (2.0 * area_cyl)
    t_cyl = math.sqrt(h0) / c_cyl  # 筒が空になる時間（= 2 T。筒は 🏺 の 1.5 倍の水を持ち、下では遅い）
    net = water_clock("cylinder")
    res = net.run(t_cyl * 1.02, trace=["clock"], sample_every=1)
    tr = res["trace"]["clock"]
    marks_c = _marks(tr, net.dt, t_cyl, n_marks, h0)
    gaps_c = [a_ - b_ for a_, b_ in zip(marks_c, marks_c[1:])]
    ratio_exact = (2 * n_marks - 1)  # 最初の間隔 / 最後の間隔 = (2n−1)/1
    dev_c = _linear_dev(tr, net.dt, h0, h0 / t_cyl, 0.9 * t_cyl)
    checks += [
        close("cylinder:gap-ratio", "🥫", gaps_c[0] / gaps_c[-1], ratio_exact, 2e-2,
              note="筒の目盛り: 最初の間隔 / 最後の間隔 = (2n−1)/1（h ∝ (1 − t/T)²）"),
        check("cylinder:gaps-shrink", "📉", all(x > y for x, y in zip(gaps_c, gaps_c[1:])), value=round(gaps_c[-1] / gaps_c[0], 4),
              note="筒は時間が進むほど目盛りが詰まる（上ほど速く下がる）"),
        check("cylinder:not-linear", "〽️", dev_c > 50 * dev, value=dev_c, expected=dev, note="筒は直線から大きく外れる（🏺 の 50 倍以上）"),
        conserved("cylinder:conserve", res),
    ]
    runs["cylinder"] = res["series"]
    summary["cylinder"] = {"t_empty": t_cyl, "gap_ratio": gaps_c[0] / gaps_c[-1], "gap_ratio_exact": ratio_exact,
                           "linear_dev": dev_c, "marks": marks_c}

    # 🔺 切頭円すい（古代の形）: ほぼ線形だが、4 次の放物面ほどではない
    net = water_clock("cone")
    res = net.run(t_empty * 1.5, trace=["clock"], sample_every=1)
    tr = res["trace"]["clock"]
    t_cone = next(((i + 1) * net.dt for i, h in enumerate(tr) if h <= 0.0), None)
    dev_cone = _linear_dev(tr, net.dt, h0, h0 / t_cone, 0.9 * t_cone) if t_cone else math.inf
    checks += [
        check("cone:between", "🔺", dev < dev_cone < dev_c, value=dev_cone, note="円すいの直線からのずれは 🏺 と 🥫 の間"),
        conserved("cone:conserve", res),
    ]
    runs["cone"] = res["series"]
    summary["cone"] = {"t_empty": t_cone, "linear_dev": dev_cone}

    # 📐 収束: dt を小さくするとずれが 1 次で減る。表の刻み（samples）はそれより小さい
    conv = {}
    for dt in dt_list:
        net = water_clock("clepsydra", dt=dt)
        conv[dt] = _linear_dev(net.run(0.9 * t_empty, trace=["clock"], sample_every=1)["trace"]["clock"], dt, h0, k, 0.9 * t_empty)
    devs = [conv[dt] for dt in dt_list]
    orders = [math.log(devs[i] / devs[i + 1]) / math.log(dt_list[i] / dt_list[i + 1]) for i in range(len(dt_list) - 1)]
    tbl = {}
    for n_s in samples_list:
        net = water_clock("clepsydra", dt=0.25, samples=n_s)
        tbl[n_s] = _linear_dev(net.run(0.9 * t_empty, trace=["clock"], sample_every=1)["trace"]["clock"], 0.25, h0, k, 0.9 * t_empty)
    checks += [
        check("converge:dt", "📐", all(x > y for x, y in zip(devs, devs[1:])) and all(0.7 <= o <= 1.3 for o in orders),
              value=[round(o, 3) for o in orders], expected=1.0, note="dt を 1/4 にするとずれも約 1/4（前進 Euler の 1 次）"),
        check("converge:samples", "🧮", tbl[samples_list[-1]] <= tbl[samples_list[0]] and tbl[samples_list[-1]] < 2 * conv[0.25],
              value={str(n_s): v for n_s, v in tbl.items()}, note="表の刻みを細かくしても dt の分より小さく、効きは頭打ち"),
    ]
    summary["convergence"] = {"dt": {str(dt): conv[dt] for dt in dt_list}, "order": orders,
                              "samples": {str(n_s): v for n_s, v in tbl.items()}}

    # 🎲 穴の大きさを A×10^-N で振っても線形（dt は T/8000 に合わせる）
    sweep = {}
    for item in random_scales(seed, 3, n_min=3, n_max=3, digits=2):
        d = item["value"]
        kk = clock_rate(h0, r_top, d, CLOCK["cd"])
        net = water_clock("clepsydra", orifice_d=d)
        dv = _linear_dev(net.run(0.9 * h0 / kk, trace=["clock"], sample_every=1)["trace"]["clock"], net.dt, h0, kk, 0.9 * h0 / kk)
        sweep[item["label"]] = {"t_empty": h0 / kk, "linear_dev": dv}
        checks.append(check(f"sweep:{item['label']}", "🎲", dv <= 2e-4, value=dv, tol=2e-4,
                            note=f"穴の直径 {item['label']} m でも直線（T = {h0 / kk:.0f} s）"))
    summary["sweep"] = sweep

    # 🧪 出口が鋭い穴でなく細い管（層流、Q ∝ h）なら、4 次の放物面でも一定には下がらない（Høyrup 1998 の注意）。
    #    A(h)·dh/dt = −c·h → h = (√h0 − κt/2)²: 筒 + 穴と同じ形で、目盛りの比は (2n−1):1 に戻る
    r_tube, l_tube = 0.4e-3, 0.05
    net = Network(4.0).tank("sea", 1.0).tank("clock", h0, head=h0, profile=clepsydra_area(h0, r_top))
    pipe = net.pipe("clock", "sea", radius=r_tube, length=l_tube)
    kappa = pipe["conductance"] * math.sqrt(h0) / area_cyl
    t_tube = 2.0 * math.sqrt(h0) / kappa
    res = net.run(t_tube * 1.02, trace=["clock"], sample_every=1)
    marks_t = _marks(res["trace"]["clock"], net.dt, t_tube, n_marks, h0)
    gaps_t = [a_ - b_ for a_, b_ in zip(marks_t, marks_t[1:])]
    checks += [
        close("tube:gap-ratio", "🧪", gaps_t[0] / gaps_t[-1], ratio_exact, 3e-2,
              note="細い管の出口（Q ∝ h）だと 🏺 でも目盛りが (2n−1):1 に詰まる（筒 + 穴と同じ形）"),
        check("tube:laminar", "〰️", res["links"][0]["laminar"], value=res["links"][0]["reynolds_max"], expected=RE_LAMINAR,
              note="管の中は層流（Poiseuille の前提）"),
        conserved("tube:conserve", res),
    ]
    runs["tube"] = res["series"]
    summary["tube"] = {"t_empty": t_tube, "gap_ratio": gaps_t[0] / gaps_t[-1], "reynolds_max": res["links"][0]["reynolds_max"]}
    return {"params": {**CLOCK, "dt_list": list(dt_list), "samples_list": list(samples_list), "seed": seed},
            "summary": summary, "runs": runs, "checks": checks}


# --- E18 斜面に降る雨（運動波の閉じた式との精度比べ） --------------------------------

PLANE = {"length": 30.0, "width": 1.0, "slope": 0.05, "max_depth": 0.01}  # 幅 1 m・長さ 30 m・勾配 5 % の板
RAIN = {"manning": 1e-4, "laminar": 3e-6}  # 雨の強さ [m/s]（360 mm/h・10.8 mm/h）。Re > 2000 / Re < 500 になる強さ
RAIN_DT = {"manning": 0.1, "laminar": 0.5}
RAIN_CELLS = (10, 20, 40, 80)
SKY = 10.0  # 空の水 [m³]（降らせる分より十分多く）


def rain_plane(intensity: float, cells: int, *, law: str = "manning", dt: float | None = None, fluid: str = "water",
               t_rain: float | None = None, t_end: float | None = None, sample_every: int = 1,
               mu: float | None = None, mu_scale: float | None = None) -> tuple[Network, dict, dict]:
    """空（sky）→ 雨 → 板（cells 区間、長方形断面）→ 受け皿（basin）。t_rain まで降らせ t_end まで流す。
    戻り値は (回路, run の結果, 閉じた式の定常)。既定では t_rain = 3·t_e、t_end = 6·t_e。"""
    sec = Section("rect", PLANE["width"])
    net = Network(RAIN_DT[law] if dt is None else dt, fluid, mu=mu, mu_scale=mu_scale)
    net.tank("sky", 1.0, volume=SKY).tank("basin", 1.0)
    net.channel("plane", "basin", length=PLANE["length"], cells=cells, slope=PLANE["slope"], section=sec, law=law,
                max_depth=PLANE["max_depth"])
    eq = overland_equilibrium(sec, intensity, PLANE["length"], PLANE["slope"], law, nu=net.fluid.nu,
                              h_max=PLANE["max_depth"])
    if eq["h_e"] >= 0.95 * PLANE["max_depth"]:
        raise ValueError(f"rain {intensity:g} m/s fills the plane past the rating table ({eq['h_e']:.4g} m)")
    t_rain = 3.0 * eq["t_e"] if t_rain is None else t_rain
    t_end = 2.0 * t_rain if t_end is None else t_end
    net.rain("plane", intensity, "sky", until=t_rain)
    res = net.run(t_end, sample_every=sample_every)
    eq["t_rain"] = t_rain
    return net, res, eq


def _outlet_series(res: dict) -> tuple[list[float], list[float]]:
    """出口（板の最後の区間 → 受け皿）の流量の時系列 (t, Q)。"""
    idx = next(ln["index"] for ln in res["links"] if ln["kind"] == "rating" and ln["dst"] == "basin")
    return [row["t"] for row in res["series"]], [row[f"q:{idx}"] for row in res["series"]]


def _at(ts: list[float], qs: list[float], t: float) -> float:
    """時系列を時刻 t で線形補間。"""
    for (t0, q0), (t1, q1) in zip(zip(ts, qs), zip(ts[1:], qs[1:])):
        if t0 <= t <= t1:
            return q0 if t1 == t0 else q0 + (q1 - q0) * (t - t0) / (t1 - t0)
    return qs[-1]


def _cross(ts: list[float], qs: list[float], level: float, rising: bool = True, after: float = 0.0) -> float:
    """流量が初めて level を越える（rising）／下回る時刻（線形補間）。届かなければ inf。"""
    for (t0, q0), (t1, q1) in zip(zip(ts, qs), zip(ts[1:], qs[1:])):
        if t1 < after:
            continue
        hit = (q0 < level <= q1) if rising else (q0 > level >= q1)
        if hit:
            return t0 + (t1 - t0) * (level - q0) / (q1 - q0)
    return math.inf


def _last_rain_row(res: dict, t_rain: float) -> dict:
    """雨がやむ直前の記録行（定常を読む）。"""
    return [row for row in res["series"] if row["t"] <= t_rain - 0.5 * res["dt"]][-1]


def rain_accuracy(law: str, intensity: float, cells: int, fluid: str = "water", dt: float | None = None) -> dict:
    """1 回分の精度: 定常・95 % 到達時刻・t_e/2 の流量・立ち上がりの形（0.2〜0.8 t_e）・減衰の、閉じた式との相対誤差。"""
    net, res, eq = rain_plane(intensity, cells, law=law, fluid=fluid, dt=dt)
    sec = Section("rect", PLANE["width"])
    kw = {"nu": net.fluid.nu}
    ts, qs = _outlet_series(res)
    q_e, t_e, t_rain = eq["Q_e"], eq["t_e"], eq["t_rain"]
    steady = _last_rain_row(res, t_rain)
    idx = next(ln["index"] for ln in res["links"] if ln["kind"] == "rating" and ln["dst"] == "basin")
    h95 = normal_depth(sec, 0.95 * q_e, PLANE["slope"], law, h_max=PLANE["max_depth"], **kw)
    t95_exact, t95 = h95 / intensity, _cross(ts, qs, 0.95 * q_e)
    q_half_exact = overland_rise(sec, intensity, 0.5 * t_e, PLANE["slope"], law, **kw)
    q_half = _at(ts, qs, 0.5 * t_e)
    rise = []
    for f in (0.2, 0.35, 0.5, 0.65, 0.8):
        exact = overland_rise(sec, intensity, f * t_e, PLANE["slope"], law, **kw)
        rise.append(rel_err(_at(ts, qs, f * t_e), exact))
    rec = []
    for f in (0.5, 1.0, 2.0):
        exact = overland_recession(sec, intensity, PLANE["length"], PLANE["slope"], law, f * t_e,
                                   h_max=PLANE["max_depth"], **kw)
        rec.append({"t_after": f * t_e, "engine": _at(ts, qs, t_rain + f * t_e), "exact": exact,
                    "rel_err": rel_err(_at(ts, qs, t_rain + f * t_e), exact)})
    t05_fall = _cross(ts, qs, 0.05 * q_e, rising=False, after=t_rain) - t_rain
    ch = res["channels"]["plane"]
    drained_q = res["volumes_q"]["basin"] + sum(res["volumes_q"][c] for c in ch_cells(res))
    return {
        "cells": cells, "fluid": fluid, "intensity": intensity, "dt": res["dt"], "t_e": t_e, "t_rain": t_rain,
        "Q_e": q_e, "h_e": eq["h_e"], "outflow_steady": steady[f"q:{idx}"], "depth_steady": steady["h:plane"][-1],
        "t95": t95, "t95_exact": t95_exact, "t95_rel_err": rel_err(t95, t95_exact),
        "q_half": q_half, "q_half_exact": q_half_exact, "q_half_rel_err": rel_err(q_half, q_half_exact),
        "rise_max_rel_err": max(rise), "recession": rec, "rec_max_rel_err": max(r["rel_err"] for r in rec),
        "t05_fall": t05_fall, "courant": ch["courant"], "overflow": ch["overflow"],
        "rained_q": ch["rained_q"], "drained_q": drained_q, "conservation_error_q": res["conservation_error_q"],
        "volumes_q": res["volumes_q"], "series": res["series"],
    }


def ch_cells(res: dict) -> list[str]:
    return [k for k in res["volumes_q"] if k.startswith("plane:")]


def _steady_checks(tag: str, icon: str, a: dict) -> list[dict]:
    return [
        close(f"{tag}:outflow", icon + "🔁", a["outflow_steady"], a["Q_e"], 1e-3, note="定常の流出 = i·B·L（雨がやむ直前）"),
        close(f"{tag}:depth", icon + "📏", a["depth_steady"], a["h_e"], 1e-3, note="出口の水深 = h_e（Q(h_e) = i·B·L）"),
        check(f"{tag}:courant", "🧭", a["courant"] <= 1.0, value=a["courant"], expected=1.0, note="Courant ≤ 1"),
        check(f"{tag}:rained=drained", "🌧️", a["rained_q"] == a["drained_q"], value=a["drained_q"],
              expected=a["rained_q"], tol=0, note="降った量 = 受け皿 + 板に残る水（ulp まで）"),
        check(f"{tag}:conserve", "⚖️", a["conservation_error_q"] == 0, value=a["conservation_error_q"], expected=0,
              tol=0, note="Σ水槽 + 管の中（Q 整数、ulp）"),
    ]


def e18_rain_runoff(cells_list=RAIN_CELLS, base_cells: int = 20, seed: int = 18) -> dict:
    """🌧️ 板に一様な雨。運動波の閉じた式（定常 Q_e = i·B·L、h_e、t_e = h_e/i、立ち上がり Q(i·t)、減衰の特性線）と比べ、
    区間を細かくすると誤差が縮むこと（収束）を精度表に残す。

    誤差の出どころは 2 つ: 区間に分ける（上流側の区間ほど水がたまりにくく、先端がなまる → t95・減衰の誤差 ∝ Δx）と、
    ステップに分ける（雨はそのステップの流出を決めた後に入るので、水深が dt だけ遅れる → Q(t_e/2) の誤差 ≈ m·dt/(t_e/2)）。
    Manning は粘度によらず、層流は t_e ∝ ν^(1/3)。"""
    checks, runs, summary = [], {}, {}
    sec = Section("rect", PLANE["width"])
    w = get_fluid("water")
    table, by_cells = [], {}
    for n in cells_list:
        a = rain_accuracy("manning", RAIN["manning"], n)
        by_cells[n] = a
        table.append({k: a[k] for k in ("cells", "dt", "t95_rel_err", "q_half_rel_err", "rise_max_rel_err",
                                        "rec_max_rel_err", "courant")})
        checks += _steady_checks(f"manning:{n}", "", a)
        runs[f"manning:{n}"] = a["series"][:: max(1, len(a["series"]) // 400)]
    fine = by_cells[max(cells_list)]
    base = by_cells[base_cells]
    errs = [row["t95_rel_err"] for row in table]
    errs_r = [row["rec_max_rel_err"] for row in table]
    ratios = [y / x for x, y in zip(errs, errs[1:])]
    half_dt = rain_accuracy("manning", RAIN["manning"], base_cells, dt=0.5 * RAIN_DT["manning"])
    table.append({k: half_dt[k] for k in ("cells", "dt", "t95_rel_err", "q_half_rel_err", "rise_max_rel_err",
                                          "rec_max_rel_err", "courant")})
    summary["accuracy"] = table
    alpha, m = overland_power_law(PLANE["slope"], "manning", nu=w.nu)
    nf = max(cells_list)
    checks += [
        close(f"manning:{nf}:t95", "⏱️", fine["t95"], fine["t95_exact"], 0.05,
              note="95 % に届く時刻 = h95/i（最も細かい区間。粗いほど遅れる: 精度表）"),
        check(f"manning:{nf}:rise-shape", "📈", fine["rise_max_rel_err"] <= 0.02, value=fine["rise_max_rel_err"],
              tol=0.02, note="立ち上がり Q(i·t)（0.2〜0.8 t_e）との最大相対ずれ"),
        check(f"manning:{nf}:recession", "📉", fine["rec_max_rel_err"] <= 0.05, value=fine["rec_max_rel_err"],
              tol=0.05, note="減衰（特性線の陰関数解、やんで 0.5〜2 t_e 後）との最大相対ずれ"),
        check(f"manning:{base_cells}:recession-slower", "🐢", base["t05_fall"] > base["t95"], value=base["t05_fall"],
              expected=base["t95"], note="5 % まで引く時間 > 95 % まで上がる時間（減衰は立ち上がりより遅い）"),
        close("manning:t_e-power-law", "📐", overland_t_e_power_law(RAIN["manning"], PLANE["length"], alpha, m),
              base["t_e"], 2 * base["h_e"] / PLANE["width"],
              note="広い平面の t_e = (L/(α i^(m−1)))^(1/m)（R ≈ h とみなす分 2h/B までずれる）"),
        check("manning:converges-dx", "📉", all(r <= 0.7 for r in ratios) and all(y < x for x, y in zip(errs_r, errs_r[1:])),
              value=errs, note="区間を倍にするごとに t95 の誤差が 0.7 倍未満に縮む（減衰の誤差も縮む）: 1 次の収束"),
        close("manning:q-half-dt", "⏲️", base["q_half_rel_err"], m * base["dt"] / (0.5 * base["t_e"]), 0.3,
              note="Q(t_e/2) の誤差 ≈ m·dt/(t_e/2)（雨がそのステップの流出の後に入る分）。区間数には依らない"),
        check("manning:q-half-halves", "⏲️", half_dt["q_half_rel_err"] < 0.6 * base["q_half_rel_err"],
              value=half_dt["q_half_rel_err"], expected=base["q_half_rel_err"], note="dt を半分にすると Q(t_e/2) の誤差も半分に"),
    ]
    # 🍶 Manning は粘度によらない
    e = rain_accuracy("manning", RAIN["manning"], base_cells, fluid="ethanol")
    checks.append(check("manning:fluid-free", "🟰", e["volumes_q"] == base["volumes_q"],
                        note="Manning の板は 🍶 でも各区間・受け皿の体積がビット単位で同じ"))
    # 〰️ 層流: 水とエタノール
    lam = {}
    for key in ("water", "ethanol"):
        fl = get_fluid(key)
        b = rain_accuracy("laminar", RAIN["laminar"], base_cells, fluid=key)
        lam[key] = b
        re = channel_reynolds(sec, b["h_e"], b["Q_e"], fl.nu)
        checks += _steady_checks(f"laminar:{key}", fl.icon, b)
        checks += [
            check(f"laminar:{key}:regime", "〰️", re < 500.0, value=re, expected=500.0, note="Re = v·R/ν < 500（層流の前提）"),
            check(f"laminar:{key}:rise-shape", "📈", b["rise_max_rel_err"] <= 0.05, value=b["rise_max_rel_err"], tol=0.05,
                  note=f"立ち上がり Q(i·t) との最大相対ずれ（{base_cells} 区間）"),
        ]
        runs[f"laminar:{key}"] = b["series"][:: max(1, len(b["series"]) // 400)]
        summary[f"laminar:{key}"] = {k: b[k] for k in ("t_e", "h_e", "t95", "t95_exact", "t95_rel_err",
                                                       "rise_max_rel_err", "rec_max_rel_err")}
    ratio = (get_fluid("ethanol").nu / w.nu) ** (1.0 / 3.0)
    checks += [
        close("laminar:t_e-ratio", "🍶⏱️", lam["ethanol"]["t_e"] / lam["water"]["t_e"], ratio, 1e-3,
              note="閉じた式: t_e ∝ ν^(1/3)（q = g·S·h³/(3ν)）"),
        close("laminar:t95-ratio", "🍶🏎️", lam["ethanol"]["t95"] / lam["water"]["t95"], ratio, 2e-2,
              note="回した 95 % 到達時刻の比も ν^(1/3) に沿う"),
    ]
    # 🎲 雨の強さを A×10^-5 で振る（Manning、定常だけ）
    sweep = []
    for item in random_scales(seed, 3, n_min=5, n_max=5):
        a = rain_accuracy("manning", item["value"], base_cells)
        sweep.append({"label": item["label"], "Q_e": a["Q_e"], "outflow": a["outflow_steady"], "h_e": a["h_e"],
                      "depth": a["depth_steady"], "t95_rel_err": a["t95_rel_err"]})
        checks += _steady_checks(f"sweep:{item['label']}", "🎲", a)[:2] + [
            check(f"sweep:{item['label']}:conserve", "⚖️", a["conservation_error_q"] == 0,
                  value=a["conservation_error_q"], expected=0, tol=0)]
    summary["sweep"] = sweep
    summary["manning"] = {k: base[k] for k in ("t_e", "h_e", "Q_e", "t95", "t95_exact", "q_half", "q_half_exact",
                                              "rise_max_rel_err", "recession", "t05_fall")}
    summary["manning_fine"] = {k: fine[k] for k in ("cells", "t95", "t95_exact", "rise_max_rel_err", "rec_max_rel_err")}
    summary["power_law"] = {"alpha": alpha, "m": m, "t_e": overland_t_e_power_law(RAIN["manning"], PLANE["length"], alpha, m)}
    return {"params": {"plane": PLANE, "rain": RAIN, "dt": RAIN_DT, "cells": list(cells_list), "base_cells": base_cells,
                       "seed": seed, "sky": SKY},
            "summary": summary, "runs": runs, "checks": checks}


# --- E23 慣性が要る題材（このエンジンでは再現しない） --------------------------------

JOUKOWSKY_C = 1400.0  # 管の中の水を走る圧力波の速さ [m/s]（自由な水で約 1480、管の伸びで下がる。wermac: 21 °C で 1460）
FRONT_FRACTION = 0.1  # 先端の読み方: 水深が h0 のこの割合を越えた時刻（Ritter の先端は深さ 0 の尖りなので同じ等深線で比べる）


@register_law("level-pipe")
def _law_level_pipe(net, ln, n):
    """水位差で流れる細管（U 字管の底）: Hagen–Poiseuille Q = c·(h_src − h_dst)。逆向きは対の管が受け持つ。
    等しい面積 A の 2 槽の差 x は dx/dt = −x/τ、τ = A/(2c) なので、その厳密離散化（1 − e^(−dt/τ)）で
    A/2·x の分を動かす。慣性の項（ẍ）はこのエンジンには無い。"""
    dh = net.head_q(ln["src"]) - net.head_q(ln["dst"])
    if dh <= 0:
        return 0
    return qmul(qmul(dh, ln["half_area_q"]), ln["k"])


def utube_theory(area: float, radius: float, length: float, limb: float, nu: float, g: float = G) -> dict:
    """U 字管（足の面積 A、底の管の半径 r・長さ L、足の水の高さ limb）の閉じた式。

    変位 y（片方の水面の上がり）の運動方程式（エネルギー法、底の管は Poiseuille の損失）:
        L_eff·ÿ + b·ẏ + 2g·y = 0、 L_eff = 2·limb + L·A/a、 b = 8νLA/(a r²)、 a = πr²
    一様な管（A = a）なら L_eff は水の全長で、T0 = 2π√(L_eff/2g) = 2π√(半分の長さ/g)（schoolphysics）。
    減衰比 ζ = b/(2√(2g·L_eff))。ζ < 1 なら実物は揺れる（減衰振動）、ζ ≥ 1 なら揺れずに落ち着く。
    エンジンは ÿ の項を持たない b·ẏ + 2g·y = 0 なので τ_engine = b/(2g) = A/(2c) の指数関数だけ。
    """
    a = math.pi * radius ** 2
    c = math.pi * radius ** 4 * g / (8.0 * nu * length)
    l_eff = 2.0 * limb + length * area / a
    b = 8.0 * nu * length * area / (a * radius ** 2)
    omega0 = math.sqrt(2.0 * g / l_eff)
    zeta = b / (2.0 * math.sqrt(2.0 * g * l_eff))
    out = {"area": area, "radius": radius, "length": length, "limb": limb, "nu": nu, "pipe_area": a,
           "conductance": c, "tau_engine": area / (2.0 * c), "L_eff": l_eff, "b": b, "omega0": omega0,
           "T0": 2.0 * math.pi / omega0, "zeta": zeta, "oscillates": zeta < 1.0}
    if zeta < 1.0:
        wd = omega0 * math.sqrt(1.0 - zeta * zeta)
        out.update(T_damped=2.0 * math.pi / wd, decay_time=1.0 / (zeta * omega0), tau_slow=None, tau_fast=None)
    else:
        disc = math.sqrt(b * b - 8.0 * g * l_eff)
        s_slow, s_fast = (-b + disc) / (2.0 * l_eff), (-b - disc) / (2.0 * l_eff)
        out.update(T_damped=None, decay_time=None, tau_slow=-1.0 / s_slow, tau_fast=-1.0 / s_fast)
    return out


def utube(area: float, radius: float, length: float, limb: float, delta_h: float, dt: float,
          fluid: str = "water") -> tuple[Network, dict]:
    """左右の槽（面積 A）を底の細管でつないだ U 字管。左を limb + Δh/2、右を limb − Δh/2 から始める。"""
    if abs(delta_h) / 2.0 > limb:
        raise ValueError("delta_h / 2 must not exceed limb (a limb would run dry)")
    net = Network(dt, fluid).tank("left", area, head=limb + delta_h / 2.0).tank("right", area, head=limb - delta_h / 2.0)
    th = utube_theory(area, radius, length, limb, net.fluid.nu, net.g)
    k = lag_gain(dt, th["tau_engine"])
    for src, dst in (("left", "right"), ("right", "left")):
        net.link("level-pipe", src, dst, half_area_q=to_q(area / 2.0), k=k, tau=th["tau_engine"])
    return net, th


def utube_crossover_radius(length: float, limb: float, nu: float, g: float = G, lo: float = 1e-5, hi: float = 0.1) -> float:
    """一様な U 字管（A = a）で ζ = 1 になる半径。これより太ければ実物は揺れる。"""
    return bisection(lambda r: utube_theory(math.pi * r ** 2, r, length, limb, nu, g)["zeta"] > 1.0, lo, hi)


def ritter_depth(x: float, t: float, h0: float, g: float = G) -> float:
    """Ritter（1892）の解: 平らで摩擦のない床、下流は乾いた床。c0 = √(g h0)、−c0·t ≤ x ≤ 2c0·t で
    h = (2c0 − x/t)² / 9g。上流（x < −c0 t）は h0、先端（x > 2c0 t）は 0。ダム地点 x = 0 では 4/9·h0 のまま。"""
    if t <= 0:
        return h0 if x < 0 else 0.0
    c0 = math.sqrt(g * h0)
    if x <= -c0 * t:
        return h0
    if x >= 2.0 * c0 * t:
        return 0.0
    return (2.0 * c0 - x / t) ** 2 / (9.0 * g)


def ritter_isoline_speed(fraction: float, h0: float, g: float = G) -> float:
    """Ritter の解で水深が fraction·h0 の点が進む速さ: (2 − 3√fraction)·c0（fraction → 0 で先端の 2c0）。"""
    return (2.0 - 3.0 * math.sqrt(fraction)) * math.sqrt(g * h0)


def dam_break(h0: float = 0.1, reservoir: float = 2.0, dry: float = 4.0, cells_per_m: float = 5.0,
              width: float = 0.2, slope: float = 1e-3, n_manning: float = 0.010, dt: float = 0.02,
              fluid: str = "water") -> tuple[Network, dict]:
    """ダム決壊: 長方形の樋（幅 width）の上流 reservoir [m] を水深 h0 で満たし、下流 dry [m] は乾いた床。
    エンジンの樋は等流の表（運動波）なので勾配が要る（slope > 0）。Ritter は勾配 0・摩擦 0 の解。"""
    cells = int(round((reservoir + dry) * cells_per_m))
    n_res = int(round(reservoir * cells_per_m))
    net = _sea(Network(dt, fluid))
    ch = net.channel("bed", "sea", length=reservoir + dry, cells=cells, slope=slope,
                     section={"shape": "rect", "size": width}, law="manning", n_manning=n_manning, max_depth=2.0 * h0)
    lc = ch["cell_length"]
    for name in ch["cells"][:n_res]:
        net.nodes[name]["v"] = to_q(width * h0 * lc)  # 初期の水柱（ダムの上流）
    xs = [(i + 0.5) * lc - reservoir for i in range(cells)]  # セル中心。ダム地点が x = 0
    return net, {"cells": cells, "n_res": n_res, "cell_length": lc, "x": xs, "h0": h0, "reservoir": reservoir,
                 "dry": dry, "slope": slope, "width": width}


def _front_fit(xs: list[float], ts: list[float]) -> float | None:
    """最小二乗で位置 − 時刻の傾き [m/s]（点が 2 つ未満、または時刻が同じなら None）。"""
    pts = [(t, x) for x, t in zip(xs, ts) if t is not None]
    if len(pts) < 2:
        return None
    n = len(pts)
    mt, mx = sum(t for t, _ in pts) / n, sum(x for _, x in pts) / n
    sxx = sum((t - mt) ** 2 for t, _ in pts)
    return None if sxx == 0 else sum((t - mt) * (x - mx) for t, x in pts) / sxx


def dam_break_report(h0: float, slope: float, dt: float = 0.02, t_end: float = 30.0, fluid: str = "water") -> dict:
    net, info = dam_break(h0=h0, slope=slope, dt=dt, fluid=fluid)
    res = net.run(t_end, sample_every=1)
    xs, n_res = info["x"], info["n_res"]
    c0 = math.sqrt(net.g * h0)
    thr = FRONT_FRACTION * h0
    arrivals = []
    for j in range(n_res, info["cells"]):
        t_j = next((row["t"] for row in res["series"] if row["h:bed"][j] >= thr), None)
        arrivals.append(t_j)
    v_engine = _front_fit(xs[n_res:], arrivals)
    v_ritter = ritter_isoline_speed(FRONT_FRACTION, h0, net.g)
    t_obs = info["reservoir"] / c0  # 戻りの波がダム地点に届く前（Ritter の 4/9 がまだ成り立つ）
    row = min(res["series"], key=lambda r: abs(r["t"] - t_obs))
    prof = row["h:bed"]
    ritter = [ritter_depth(x, row["t"], h0, net.g) for x in xs]
    rms = math.sqrt(sum((a - b) ** 2 for a, b in zip(prof, ritter)) / len(xs)) / h0
    dam_engine = 0.5 * (prof[n_res - 1] + prof[n_res])
    dam_ritter = 0.5 * (ritter[n_res - 1] + ritter[n_res])
    last = arrivals[-1]
    # 貯水側の各セルが 0.9·h0 を初めて割る時刻（奥の壁 → 堤の順）。Ritter は堤の近くから、エンジンは奥から
    drop_engine = [next((row["t"] for row in res["series"] if row["h:bed"][j] < 0.9 * h0), None) for j in range(n_res)]
    drop_ritter = [abs(x) / ((3.0 * math.sqrt(0.9) - 2.0) * c0) for x in xs[:n_res]]
    reversed_order = (None not in drop_engine and all(a <= b for a, b in zip(drop_engine, drop_engine[1:]))
                      and all(a > b for a, b in zip(drop_ritter, drop_ritter[1:])))
    return {
        "h0": h0, "slope": slope, "c0": c0, "front_fraction": FRONT_FRACTION, "arrivals": arrivals,
        "drawdown_engine": drop_engine, "drawdown_ritter": drop_ritter, "drawdown_reversed": reversed_order,
        "v_engine": v_engine, "v_ritter_isoline": v_ritter, "v_ritter_tip": 2.0 * c0,
        "front_ratio": None if v_engine is None else v_engine / v_ritter,
        "tip_ratio": None if v_engine is None else v_engine / (2.0 * c0),
        "t_end_engine": last, "t_end_ritter": info["dry"] / v_ritter,
        "t_obs": row["t"], "dam_depth_engine": dam_engine, "dam_depth_ritter": dam_ritter,
        "dam_ratio": dam_engine / dam_ritter, "profile_rms_over_h0": rms,
        "courant": res["channels"]["bed"]["courant"], "res": res, "info": info,
    }


def joukowsky(delta_v: float, rho: float, c: float = JOUKOWSKY_C, g: float = G) -> dict:
    """弁を急に閉めたときの圧力の跳ね上がり Δp = ρ·c·Δv（Joukowsky）。水頭にすると Δp/(ρg) = c·Δv/g。"""
    dp = rho * c * delta_v
    return {"delta_v": delta_v, "rho": rho, "c": c, "delta_p": dp, "delta_p_bar": dp / 1e5, "head": c * delta_v / g}


def e23_inertia_limits(seed: int = 0, dt_frac: float = 0.02) -> dict:
    """このエンジン（ストックとフローの準静的な模型）が、慣性の要る題材をどれだけ外すかを数で記す。

    直さない（慣性を足せば別のエンジンになる）。検証は「外れ方が記録どおりか」を確かめる。
    1. U 字管: 実物は L_eff·ÿ + b·ẏ + 2g·y = 0。ζ < 1 なら揺れる。エンジンは指数関数で落ち着くだけ
    2. ダム決壊: Ritter の先端 2√(g h0)、ダム地点 4/9·h0。エンジンの樋は運動波（勾配が要る）
    3. 水撃: 弁を急に閉めると Δp = ρ·c·Δv の衝撃。エンジンは流れが止まるだけ（圧力が無い）
    """
    checks, runs, summary = [], {}, {}
    w = get_fluid("water")
    geometries = {
        # 名前: (A, r, L, limb). manometer は一様な管（A = a）、capillary は E02/E03 の槽と細管
        "manometer": (math.pi * 0.003 ** 2, 0.003, 0.3, 0.3),
        "capillary": (0.01, 0.75e-3, 0.5, 0.1),
    }
    expect = {"manometer": True, "capillary": False}
    for name, (area, radius, length, limb) in geometries.items():
        th = utube_theory(area, radius, length, limb, w.nu)
        checks.append(check(f"{name}:kind", "🔔" if th["oscillates"] else "🫗", th["oscillates"] == expect[name],
                            value=th["zeta"], expected=1.0,
                            note=("ζ < 1: 実物は揺れる（エンジンは種類ごと外す）" if expect[name]
                                  else "ζ ≥ 1: 実物も揺れない（エンジンで足りる）")))
        if not th["oscillates"]:
            checks.append(close(f"{name}:slow-root", "🐢", th["tau_engine"], th["tau_slow"], 1.0 / th["zeta"] ** 2,
                                note="エンジンの τ = 実物の遅い根 1/|s_slow|（ずれは 1/ζ² のオーダー）"))
            checks.append(check(f"{name}:fast-root", "⚡", th["tau_fast"] < 0.01 * th["tau_engine"],
                                value=th["tau_fast"], expected=th["tau_engine"],
                                note="慣性の過渡 1/|s_fast| は τ の 1 % 未満（無視して差し支えない）"))
        else:
            checks.append(check(f"{name}:kind-mismatch", "🚧", th["T_damped"] > th["tau_engine"] * 10,
                                value=th["T_damped"], expected=th["tau_engine"],
                                note="実物の周期 ≫ エンジンの τ: 揺れ続けるものを一瞬で落ち着かせてしまう"))
        traces = {}
        sweep = random_scales(seed + 1, 3, n_min=2, n_max=2)
        dh_max = max(i["value"] for i in sweep)
        u_peak = (area / th["pipe_area"]) * (dh_max / 2.0) * (th["omega0"] if th["oscillates"] else 1.0 / th["tau_slow"])
        re_peak = u_peak * 2.0 * radius / w.nu
        checks.append(check(f"{name}:laminar", "〰️", re_peak < RE_LAMINAR, value=re_peak, expected=RE_LAMINAR,
                            note="底の管の最大流速の Re < 2300（Poiseuille の減衰で見積もってよい）"))
        for item in sweep:
            dh = item["value"]
            dt = th["tau_engine"] * dt_frac
            net, _th = utube(area, radius, length, limb, dh, dt)
            res = net.run(8.0 * th["tau_engine"], trace=["left", "right"])
            diff = [a - b for a, b in zip(res["trace"]["left"], res["trace"]["right"])]
            monotone = all(b <= a + 1e-15 for a, b in zip(diff, diff[1:])) and all(d >= -1e-15 for d in diff)
            tag = f"{name}:{item['label']}"
            checks += [
                check(f"{tag}:monotone", "📉", monotone, value=min(diff), expected=0.0,
                      note="水位差は符号を変えず単調に減る（揺れは一度も出ない）"),
                check(f"{tag}:settle", "🎯", abs(diff[-1]) <= 1e-3 * dh, value=diff[-1], tol=1e-3 * dh,
                      note="8τ 後に差は 0.1 % 未満"),
                conserved(f"{tag}:conserve", res),
            ]
            if dh == dh_max:
                traces = {"t": [(n + 1) * dt for n in range(len(diff))], "diff": diff, "delta_h": dh}
        runs[f"utube:{name}"] = [{"t": t, "diff": d} for t, d in zip(traces["t"][::25], traces["diff"][::25])]
        summary[f"utube:{name}"] = dict(th) | {"no_oscillation_reproduced": True, "reynolds_peak": re_peak,
                                               "delta_h_max": dh_max}
    # 一様な U 字管の太さで、実物が揺れるかどうかの境目
    r_star = utube_crossover_radius(0.3, 0.3, w.nu)
    radii = random_scales(seed + 2, 4, n_min=3, n_max=3)  # 1〜10 mm
    table = []
    for item in radii:
        r = item["value"]
        th = utube_theory(math.pi * r * r, r, 0.3, 0.3, w.nu)
        table.append({"radius": r, "zeta": th["zeta"], "oscillates": th["oscillates"], "tau_engine": th["tau_engine"],
                      "T0": th["T0"]})
    checks.append(check("crossover", "📐", all((row["radius"] > r_star) == row["oscillates"] for row in table),
                        value=r_star, note="一様な管（L 0.3 m・足 0.3 m）: 半径がこれより太ければ実物は揺れる"))
    summary["utube:crossover_radius"] = r_star
    summary["utube:radius_table"] = table

    # 🌊 ダム決壊（運動波の樋 vs Ritter）
    for slope in (1e-3, 1e-2):
        rep = dam_break_report(0.1, slope)
        tag = f"dam:S{slope:g}"
        checks += [
            check(f"{tag}:front-off", "🐌" if rep["front_ratio"] < 1 else "🐇",
                  rep["front_ratio"] is not None and abs(rep["front_ratio"] - 1.0) > 0.2, value=rep["front_ratio"],
                  expected=1.0, tol=0.2,
                  note=f"水深 {FRONT_FRACTION:g}·h0 の等深線の速さ: エンジン / Ritter（{rep['v_ritter_isoline']:.3f} m/s）。"
                       "緩い勾配では遅く、急な勾配では速い（どちらも 20 % 以上外す）"),
            check(f"{tag}:tip-ratio", "📏", rep["tip_ratio"] is not None and rep["tip_ratio"] < 1.0,
                  value=rep["tip_ratio"], expected=1.0, note="先端 2√(g h0) には届かない"),
            check(f"{tag}:dam-depth", "🏗️", abs(rep["dam_ratio"] - 1.0) > 0.2, value=rep["dam_depth_engine"],
                  expected=rep["dam_depth_ritter"], tol=0.2, note="t = L_res/c0 のダム地点の水深（Ritter は 4/9·h0 に留まる）"),
            check(f"{tag}:negative-wave", "↩️", rep["drawdown_reversed"], value=rep["drawdown_engine"][:3],
                  expected=rep["drawdown_ritter"][:3],
                  note="水位が下がり始める順: Ritter は堤の近くから奥へ（負の波が上流へ戻る）、エンジンは奥の壁から堤へ（運動波は下流にしか伝わらない）"),
            check(f"{tag}:profile", "📉", rep["profile_rms_over_h0"] > 0.05, value=rep["profile_rms_over_h0"],
                  expected=0.05, note="t = L_res/c0 での水面形の RMS ずれ / h0（放物線の水面形を出せない）"),
            check(f"{tag}:courant", "🧭", rep["courant"] <= 1.0, value=rep["courant"], expected=1.0, note="Courant ≤ 1"),
            conserved(f"{tag}:conserve", rep["res"]),
        ]
        runs[tag] = [{"t": r["t"], "h:bed": r["h:bed"]} for r in rep["res"]["series"][:: max(1, len(rep["res"]["series"]) // 8)]]
        summary[tag] = {k: v for k, v in rep.items() if k not in ("res", "info", "arrivals")}
    s1, s2 = summary["dam:S0.001"]["front_ratio"], summary["dam:S0.01"]["front_ratio"]
    checks.append(check("dam:slope-dependence", "🚧", s2 > 2.0 * s1, value=[s1, s2],
                        note="エンジンの先端の速さは勾配で 2 倍以上変わる（Ritter には勾配という量が無い）"))
    # 🟦 ドット絵のダム決壊: 先端は規則で 1 ティックに 1 + 横の回数 マスまで（物理ではない）
    fine = pw.run_pixels("dam-break", "water", seed=seed, frame_every=1)
    arr = pw.pixel_arrivals(fine)
    dry = [t for t in arr if t is not None and t > 0]
    front = [sum(1 for t in dry if t <= tick) for tick in range(fine["ticks"] + 1)]
    max_adv = max(b - a for a, b in zip(front, front[1:]))
    mean_speed = len(dry) / max(dry) if dry else 0.0
    cap = 1 + math.ceil(fine["lateral_rate"])
    checks += [
        check("pixel:cap", "🟦", max_adv <= cap, value=max_adv, expected=cap,
              note="ドット絵の先端は 1 ティックに「斜め 1 + 横の回数」マスまで（規則の上限、物理ではない）"),
        check("pixel:mean", "🟦", mean_speed <= 1.0 + fine["lateral_rate"], value=mean_speed,
              expected=1.0 + fine["lateral_rate"], note="平均の先端の速さ [マス/ティック] ≤ 1 + 横の頻度"),
    ]
    summary["pixel:dam-break"] = {"max_advance_cells_per_tick": max_adv, "mean_cells_per_tick": mean_speed,
                                  "cap": cap, "lateral_rate": fine["lateral_rate"], "note": "g も m も無い盤。速さは規則"}
    # 💥 水撃
    q, area_pipe = 1e-4, math.pi * 0.005 ** 2  # 10 mm の管に 0.1 L/s → 約 1.27 m/s
    dv = q / area_pipe
    jk = joukowsky(dv, w.rho)
    dt = 0.01
    net = Network(dt).tank("up", 0.05, head=0.5).tank("down", 0.05)
    net.pump("up", "down", rate=q, until=1.0)
    res = net.run(2.0, sample_every=1)
    flows = [r["q:0"] for r in res["series"]]
    stop = next(r["t"] for r, f in zip(res["series"], flows) if f == 0.0)
    checks += [
        check("hammer:stops-in-one-step", "🛑", abs(stop - (1.0 + dt)) <= dt / 2 + 1e-12, value=stop, expected=1.0 + dt,
              note="弁を閉めた次のステップで流量 0（圧力の跳ね上がりは起きない＝量が無い）"),
        check("hammer:ignored-surge", "💥", jk["delta_p"] > 1e6, value=jk["delta_p"], expected=1e6,
              note=f"無視される衝撃 Δp = ρ·c·Δv（Δv = {dv:.2f} m/s, c = {JOUKOWSKY_C:g} m/s）= 水頭 {jk['head']:.0f} m"),
        conserved("hammer:conserve", res),
    ]
    summary["hammer"] = jk | {"stop_time": stop}
    return {"params": {"seed": seed, "dt_frac": dt_frac, "geometries": {k: list(v) for k, v in geometries.items()},
                       "front_fraction": FRONT_FRACTION, "joukowsky_c": JOUKOWSKY_C},
            "summary": summary, "runs": runs, "checks": checks}


# --- E24 ワインの蒸留（⚗️ 単式蒸留・コニャック式の 2 回蒸留） ------------------------------

STILL = {"wine_abv": 12.0, "wine_volume": 0.025, "power": 5000.0, "brouillis_cut": 5.0, "heads_fraction": 0.01,
         "heart_cut": 60.0, "secondes_cut": 5.0}
SPECIES = ("ethanol", "water")


def _vol(parts: dict) -> float:
    return mix_volume(parts["ethanol"], parts["water"])


def still_network(ethanol: float, water: float, *, dt: float = 1.0, power: float = STILL["power"],
                  receivers=("receiver",), cuts=None, fluid: str = "water") -> Network:
    """⚗️ 釜（pot）と受け器。中身はエタノール・水のモル数。cuts を渡すと受け器を切り替える（前留・本留・後留）。"""
    net = Network(dt, fluid, species=SPECIES, volume_of=_vol)
    net.tank("pot", 0.1, parts={"ethanol": ethanol, "water": water})
    for r in receivers:
        net.tank(r, 0.05)
    if cuts:
        net.boil("pot", power=power, cuts=cuts)
    else:
        net.boil("pot", receivers[0], power=power)
    return net


def _fraction(res: dict, name: str) -> dict:
    p = res["parts_q"][name]
    e, w = from_q(p["ethanol"]), from_q(p["water"])
    return {"ethanol": e, "water": w, "litres": mix_volume(e, w) * 1e3, "abv": abv(e, w) if e + w > 0 else 0.0}


def _parts_exact(res: dict, start: dict) -> bool:
    """受け器と釜の成分の合計が、入れた量と Q 整数で一致するか。"""
    return all(sum(nd[sp] for nd in res["parts_q"].values()) == start[sp] for sp in SPECIES)


def e24_distillation(dt_list=(4.0, 1.0, 0.25), seed: int = 24) -> dict:
    """⚗️ ワインを釜で沸かす単式蒸留を、Rayleigh 式（ln W0/W = ∫ dx/(y − x)、Rayleigh 卿がエタノールと水で調べた形）と
    比べる。気液平衡は Perry's の表（101.3 kPa）、度数は 20 °C の密度表で換算（混ぜると縮む）。
    そのうえでコニャック式の 2 回蒸留（初留 → 本留で前留・ハート・セコンド・後留に切り分け）を回し、
    公表されている度数（初留 28〜32 %、ハート約 70 %）と照らす。"""
    checks, runs, summary = [], {}, {}
    ne, nw = moles_for_abv(STILL["wine_abv"], STILL["wine_volume"])
    w0, x0 = ne + nw, mole_fraction(ne, nw)

    # 📐 Rayleigh 式: 釜に残る量と組成の関係。dt を 1/4 にすると誤差も 1/4（前進 Euler の 1 次）
    conv = {}
    for dt in dt_list:
        net = still_network(ne, nw, dt=dt)
        res = net.run(3600.0, sample_every=max(1, int(60 / dt)))
        pot = _fraction(res, "pot")
        x = mole_fraction(pot["ethanol"], pot["water"])
        w_ray = w0 * math.exp(-rayleigh_log(x, x0))
        conv[dt] = abs(pot["ethanol"] + pot["water"] - w_ray) / w_ray
        checks.append(check(f"rayleigh:dt{dt:g}:species", "⚖️", res["parts_error_q"] == {"ethanol": 0, "water": 0},
                            value=res["parts_error_q"], note="エタノールと水がそれぞれ 0 ulp で保存"))
        if dt == 1.0:
            runs["rayleigh"] = res["series"]
            rows = res["series"]
            xs = [mole_fraction(r["parts:pot"]["ethanol"], r["parts:pot"]["water"]) for r in rows]
            temps = [bubble_t(x_) for x_ in xs]
            running = [abv(vle_y(x_), 1.0 - vle_y(x_)) for x_ in xs]
            summary["rayleigh"] = {"x_end": x, "pot_left_mol": pot["ethanol"] + pot["water"], "rayleigh_mol": w_ray,
                                   "t_start": temps[0], "t_end": temps[-1], "running_abv_start": running[0],
                                   "running_abv_end": running[-1]}
    errs = [conv[dt] for dt in dt_list]
    orders = [math.log(errs[i] / errs[i + 1]) / math.log(dt_list[i] / dt_list[i + 1]) for i in range(len(errs) - 1)]
    checks += [
        check("rayleigh:error", "📐", errs[-1] <= 1e-4, value=errs[-1], tol=1e-4,
              note=f"釜に残る量 W と Rayleigh 式の W のずれ（dt = {dt_list[-1]:g} s）"),
        check("rayleigh:order", "📐", all(0.8 <= o <= 1.2 for o in orders), value=[round(o, 3) for o in orders],
              expected=1.0, note="dt を 1/4 にするとずれも約 1/4（前進 Euler の 1 次）"),
        check("pot:temperature-rises", "🌡️", all(a <= b + 1e-12 for a, b in zip(temps, temps[1:])),
              value=[round(temps[0], 2), round(temps[-1], 2)], note="釜の沸点はエタノールが抜けるほど上がる（表の沸点）"),
        check("running:abv-falls", "📉", all(a >= b - 1e-12 for a, b in zip(running, running[1:])),
              value=[round(running[0], 2), round(running[-1], 2)], note="流れ出る留液の度数は下がっていく"),
    ]
    summary["convergence"] = {"dt": {str(dt): conv[dt] for dt in dt_list}, "order": orders}

    # 🎲 火力を A×10^3 W で振っても Rayleigh の関係は同じ（時間だけが 1/火力で縮む）
    sweep = {}
    for item in random_scales(seed, 3, n_min=-3, n_max=-3, digits=2):
        p_w = item["value"]
        net = still_network(ne, nw, dt=1.0, power=p_w)
        t_end = 3600.0 * STILL["power"] / p_w
        res = net.run(t_end, sample_every=max(1, int(t_end / 50)))
        pot = _fraction(res, "pot")
        x = mole_fraction(pot["ethanol"], pot["water"])
        err = abs(pot["ethanol"] + pot["water"] - w0 * math.exp(-rayleigh_log(x, x0))) / (pot["ethanol"] + pot["water"])
        sweep[item["label"]] = {"t_end": t_end, "rayleigh_error": err, "x_end": x}
        checks.append(check(f"sweep:{item['label']}", "🎲", err <= 2e-3, value=err, tol=2e-3,
                            note=f"火力 {item['label']} W でも Rayleigh 式どおり（時間は {t_end:.0f} s）"))
    summary["power_sweep"] = sweep

    # 🫧 共沸点: y = x なので、蒸留しても組成が変わらない
    na, nwa = AZEOTROPE_X * 300.0, (1.0 - AZEOTROPE_X) * 300.0
    net = still_network(na, nwa, dt=1.0)
    res = net.run(600.0)  # 600 s で釜の約 1/4 を蒸発させる（空にはしない）
    pot, rec = _fraction(res, "pot"), _fraction(res, "receiver")
    checks += [
        close("azeotrope:pot", "🫧", mole_fraction(pot["ethanol"], pot["water"]), AZEOTROPE_X, 1e-9,
              note="共沸点（モル分率 0.8943）の液は蒸留しても組成が変わらない"),
        close("azeotrope:distillate", "🫧", rec["abv"], abv(na, nwa), 1e-6, note="留液の度数 = 釜の度数（これ以上は濃くならない）"),
    ]
    summary["azeotrope"] = {"abv": abv(na, nwa), "distillate_abv": rec["abv"]}

    # 🍇 初留（première chauffe）: ワインを、流れ出る度数が brouillis_cut を切るまで受ける
    start = {sp: to_q(v) for sp, v in (("ethanol", ne), ("water", nw))}
    net = still_network(ne, nw, dt=1.0, receivers=("brouillis", "tails1"),
                        cuts=[{"to": "brouillis", "abv_below": STILL["brouillis_cut"]}, {"to": "tails1"}])
    res1 = net.run(6000.0, sample_every=60)
    b = _fraction(res1, "brouillis")
    t1 = res1["links"][0]["cuts"][0].get("ended")
    checks += [
        check("premiere:brouillis-abv", "🍇", 28.0 <= b["abv"] <= 32.0, value=b["abv"], expected=30.0, tol=2.0,
              note=f"初留の度数が公表の 28〜32 %（流れ出る度数 {STILL['brouillis_cut']:g} % で切った場合）"),
        check("premiere:exact", "⚖️", _parts_exact(res1, start), note="受け器 + 釜 = 入れたワイン（成分ごとに Q 整数で一致）"),
    ]
    runs["premiere"] = res1["series"]
    cut_table = {}
    for cut in (2.0, 3.0, 5.0, 8.0):
        r = still_network(ne, nw, dt=1.0, receivers=("b", "t"), cuts=[{"to": "b", "abv_below": cut}, {"to": "t"}]).run(
            6000.0, sample_every=600)
        fb = _fraction(r, "b")
        cut_table[f"{cut:g}"] = {"abv": fb["abv"], "litres": fb["litres"]}
    summary["premiere"] = {"brouillis": b, "cut_time": t1, "volume_ratio": b["litres"] / (STILL["wine_volume"] * 1e3),
                           "cut_table": cut_table}

    # 🥃 本留（bonne chauffe）: 初留を釜へ戻し、前留（量で 1 %）→ ハート（60 % を切るまで）→ セコンド → 後留
    heads_v = STILL["heads_fraction"] * b["litres"] * 1e-3
    net = still_network(b["ethanol"], b["water"], dt=1.0, receivers=("heads", "heart", "secondes", "tails"),
                        cuts=[{"to": "heads", "volume": heads_v}, {"to": "heart", "abv_below": STILL["heart_cut"]},
                              {"to": "secondes", "abv_below": STILL["secondes_cut"]}, {"to": "tails"}])
    start2 = {sp: to_q(b[sp]) for sp in SPECIES}
    res2 = net.run(6000.0, sample_every=60)
    fr = {k: _fraction(res2, k) for k in ("heads", "heart", "secondes", "tails")}
    checks += [
        check("bonne:heart-abv", "🥃", 68.0 <= fr["heart"]["abv"] <= 72.0, value=fr["heart"]["abv"], expected=70.0,
              tol=2.0, note="ハートの度数が公表の 68〜72 %（60 % を切ったらセコンドへ、の切り方そのまま）"),
        check("bonne:order", "🪜", fr["heads"]["abv"] > fr["heart"]["abv"] > fr["secondes"]["abv"] > fr["tails"]["abv"],
              value=[round(fr[k]["abv"], 1) for k in ("heads", "heart", "secondes", "tails")],
              note="前留 > ハート > セコンド > 後留 の順に薄くなる"),
        check("bonne:exact", "⚖️", _parts_exact(res2, start2), note="4 つの受け器 + 釜 = 入れた初留（成分ごとに Q 整数で一致）"),
        check("bonne:species", "⚖️", res2["parts_error_q"] == {"ethanol": 0, "water": 0}, value=res2["parts_error_q"],
              note="成分ごとに 0 ulp"),
    ]
    runs["bonne"] = res2["series"]
    yield_ratio = STILL["wine_volume"] * 1e3 / fr["heart"]["litres"]
    summary["bonne"] = {k: v for k, v in fr.items()} | {"cuts": res2["links"][0]["cuts"],
                                                         "wine_litres_per_heart_litre": yield_ratio}
    return {"params": {**STILL, "dt_list": list(dt_list), "seed": seed}, "summary": summary, "runs": runs,
            "checks": checks}


# --- E25 天使の取り分（👼 樽の木を抜けていく水とエタノール） -----------------------------

YEAR = 365.25 * 86400.0  # [s]
CASK = {"volume": 0.2, "fill_abv": 63.5}  # 200 L の樽、スコッチの標準の詰め度数 63.5 %（whiskipedia）
SCOTLAND = {"rh": 0.88, "loss": 0.02, "abv_drop": 0.5}  # 湿度 88 %（puninwine）、年 2 %（各社）、度数 −0.5 点/年
CLIMATES = {  # rh・loss は入力（🎯 = 較正に使う出典の数字、📝 = 置いた仮定）。f は温度の倍率（loss から合わせる）
    "scotland": {"icon": "🏴", "rh": 0.88, "loss": 0.02, "src": "🎯 湿度 88 %・年 2 %（puninwine / alcoholprofessor）"},
    "kentucky-top": {"icon": "🔥", "rh": 0.45, "loss": 0.04,
                     "src": "🎯 年 4 %（alcoholprofessor）/ 📝 湿度 45 %（上の階は暑く乾く: vinepair、数値は仮定）"},
    "kentucky-ground": {"icon": "🧱", "rh": 0.75, "loss": 0.03,
                        "src": "📝 湿度 75 %・年 3 %（下の階は湿る: vinepair、数値は仮定）"},
    "taiwan": {"icon": "🌴", "rh": 0.80, "loss": 0.11, "src": "🎯 年 10〜12 %（thewhiskeywash / milroys）/ 📝 湿度 80 %"},
}
SCOTCH_MIN_ABV = 40.0  # スコッチとして売れる下限（alcoholprofessor / puninwine）


def angel_rates(fill_abv: float = CASK["fill_abv"], volume: float = CASK["volume"], target: dict = SCOTLAND) -> dict:
    """スコットランドの 1 年目（体積 −loss、度数 −abv_drop 点）に合う k_e・k_w [1/年] を厳密に解く（🎯 較正）。

    1 年後の中身は moles_for_abv(度数, 体積) で決まるので、k_s = ln(n_s0 / n_s1) がそのまま出る。"""
    ne0, nw0 = moles_for_abv(fill_abv, volume)
    ne1, nw1 = moles_for_abv(fill_abv - target["abv_drop"], volume * (1.0 - target["loss"]))
    return {"k_e": math.log(ne0 / ne1), "k_w": math.log(nw0 / nw1), "rh_ref": target["rh"]}


def water_rate(k_w_ref: float, rh: float, rh_ref: float) -> float:
    """水の抜ける速さは空気の乾き（1 − RH）に比例する: k_w(RH) = k_w,ref·(1 − RH)/(1 − RH_ref)。
    エタノールは空気に無いので湿度によらない（k_e 一定）。"""
    return k_w_ref * max(1.0 - rh, 0.0) / (1.0 - rh_ref)


def crossover_rh(rates: dict) -> float:
    """度数が動かない湿度 RH*: k_w(RH*) = k_e → RH* = 1 − (1 − RH_ref)·k_e/k_w,ref。温度の倍率にはよらない。"""
    return 1.0 - (1.0 - rates["rh_ref"]) * rates["k_e"] / rates["k_w"]


def cask_state(t_years: float, rates: dict, rh: float, f: float = 1.0, fill_abv: float = CASK["fill_abv"],
               volume: float = CASK["volume"]) -> tuple[float, float]:
    """閉じた式 n_s(t) = n_s0·e^(−f·k_s·t) の (エタノール, 水) [mol]。"""
    ne0, nw0 = moles_for_abv(fill_abv, volume)
    kw = water_rate(rates["k_w"], rh, rates["rh_ref"])
    return ne0 * math.exp(-f * rates["k_e"] * t_years), nw0 * math.exp(-f * kw * t_years)


def fit_temperature(rates: dict, rh: float, loss: float, volume: float = CASK["volume"]) -> float:
    """1 年目の体積の減り loss に合う温度の倍率 f（二分法、🎯 較正）。"""
    v0 = mix_volume(*cask_state(0.0, rates, rh, volume=volume))
    lo, hi = 1e-3, 1e3
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        lost = 1.0 - mix_volume(*cask_state(1.0, rates, rh, mid, volume=volume)) / v0
        lo, hi = (mid, hi) if lost < loss else (lo, mid)
    return math.sqrt(lo * hi)


def years_to_abv(limit: float, rates: dict, rh: float, f: float = 1.0, t_max: float = 500.0) -> float | None:
    """度数が limit を切るまでの年数（閉じた式を二分法）。切らなければ None。"""
    if abv(*cask_state(t_max, rates, rh, f)) > limit:
        return None
    return bisection(lambda t: abv(*cask_state(t, rates, rh, f)) > limit, 0.0, t_max)


@register_law("soak")
def _law_soak(net, ln, n):
    """🪵 木が中身を吸う: 木（dst）の量が capacity に届くまで、1 ステップに (capacity − 今の量)·(1 − e^(−dt/τ))。
    成分モードでは樽の中身の比のまま吸う。"""
    room = ln["capacity_q"] - net.nodes[ln["dst"]]["v"]
    return qmul(room, ln["gain"]) if room > 0 else 0


def angel_cask(rates: dict, climates: dict, *, steps_per_year: int = 365, soak: dict | None = None) -> Network:
    """樽を気候ごとに 1 つずつ並べ、どれも空（air）へ抜ける。soak = {"capacity": mol, "tau_years": 年} で木の吸い込みも。
    刻みは 1 年を steps_per_year 等分（既定 365 でほぼ 1 日。年の区切りがちょうどステップに乗る）。"""
    net = Network(YEAR / steps_per_year, species=("ethanol", "water"), volume_of=lambda p: mix_volume(p["ethanol"], p["water"]))
    ne0, nw0 = moles_for_abv(CASK["fill_abv"], CASK["volume"])
    net.tank("air", 1.0)
    for name, c in climates.items():
        net.tank(name, 0.25, parts={"ethanol": ne0, "water": nw0})
        kw = water_rate(rates["k_w"], c["rh"], rates["rh_ref"])
        net.evaporate(name, "air", rates={"ethanol": c["f"] * rates["k_e"] / YEAR, "water": c["f"] * kw / YEAR})
        if soak:
            net.tank(f"{name}:wood", 1.0)
            net.link("soak", name, f"{name}:wood", capacity_q=to_q(soak["capacity"]),
                     gain=lag_gain(net.dt, soak["tau_years"] * YEAR))
    return net


def _cask_series(res: dict, name: str) -> list[tuple[float, float, float]]:
    """(年, 体積 m³, 度数) の並び（0 年目を含む）。"""
    ne0, nw0 = moles_for_abv(CASK["fill_abv"], CASK["volume"])
    out = [(0.0, mix_volume(ne0, nw0), CASK["fill_abv"])]
    for row in res["series"]:
        p = row[f"parts:{name}"]
        out.append((row["t"] / YEAR, mix_volume(p["ethanol"], p["water"]), abv(p["ethanol"], p["water"])))
    return out


def e25_angels_share(years: int = 50, steps_per_year: int = 365, seed: int = 25) -> dict:
    """👼 天使の取り分。🎯 スコットランドの数字で 2 つの速さを合わせ（較正）、同じ式で他の気候を予測する。

    検証できるのは「式どおりに動くか」（閉じた式・成分ごとの保存・dt によらない）で、物理としての正しさではない。
    湿度で度数の上下が入れ替わる・その湿度 RH*・40 % を切る年数は、較正に使っていない出力。"""
    checks, runs, summary = [], {}, {}
    rates = angel_rates()
    rh_star = crossover_rh(rates)
    climates = {k: {**c, "f": (1.0 if k == "scotland" else fit_temperature(rates, c["rh"], c["loss"]))}
                for k, c in CLIMATES.items()}
    net = angel_cask(rates, climates, steps_per_year=steps_per_year)
    res = net.run(years * YEAR, sample_every=steps_per_year)
    v0 = CASK["volume"]
    checks.append(check("conserve:species", "⚖️", res["parts_error_q"] == {"ethanol": 0, "water": 0},
                        value=res["parts_error_q"], note="成分ごとの総量（樽 + 空）が 0 ulp"))
    for name, c in climates.items():
        ser = _cask_series(res, name)
        y1 = ser[1]
        loss1 = 1.0 - y1[1] / v0
        e_t, w_t = cask_state(years, rates, c["rh"], c["f"])
        p = res["parts_q"][name]
        e_s, w_s = from_q(p["ethanol"]), from_q(p["water"])
        rising = ser[-1][2] > CASK["fill_abv"]
        checks += [
            close(f"{name}:closed-form", c["icon"] + "📐", e_s, e_t, 1e-8, note=f"{years} 年後のエタノール = n0·e^(−f·k_e·t)"),
            close(f"{name}:closed-form-water", c["icon"] + "📐", w_s, w_t, 1e-8, note="水 = n0·e^(−f·k_w(RH)·t)"),
            close(f"{name}:loss-1y", c["icon"] + "🎯", loss1, c["loss"], 1e-6,
                  note="1 年目の体積の減り（🎯 較正どおりか。これは合わせた数なので予測ではない）"),
            check(f"{name}:abv-direction", c["icon"] + ("📈" if c["rh"] < rh_star else "📉"),
                  rising == (c["rh"] < rh_star), value=round(ser[-1][2], 3), expected=CASK["fill_abv"],
                  note=f"湿度 {c['rh']:.0%} が RH* = {rh_star:.1%} より{'低い → 度数が上がる' if c['rh'] < rh_star else '高い → 度数が下がる'}（🔮 予測）"),
        ]
        runs[name] = [{"year": round(t, 6), "liters": v * 1e3, "abv": a} for t, v, a in ser]
        summary[name] = {"rh": c["rh"], "f": c["f"], "loss_1y": loss1, "abv_1y": y1[2], "abv_end": ser[-1][2],
                         "liters_end": ser[-1][1] * 1e3, "source": c["src"]}
    sc = _cask_series(res, "scotland")
    checks.append(close("scotland:abv-1y", "🏴🎯", sc[0][2] - sc[1][2], SCOTLAND["abv_drop"], 1e-6,
                        note="1 年目の度数の下がり 0.5 点（🎯 較正どおり）"))
    # 📉 本場の向き: スコットランドは下がり、ケンタッキー（上の階）は上がる（alcoholprofessor・vinepair の向き）
    checks += [
        check("sign:scotland-falls", "🏴📉", summary["scotland"]["abv_end"] < CASK["fill_abv"],
              value=summary["scotland"]["abv_end"], note="湿ったスコットランドでは度数が下がる（出典どおりの向き）"),
        check("sign:kentucky-rises", "🔥📈", summary["kentucky-top"]["abv_end"] > CASK["fill_abv"],
              value=summary["kentucky-top"]["abv_end"], note="乾いた上の階では度数が上がる（湿度 45 % は仮定）"),
        check("sign:same-warehouse", "🧱📉", summary["kentucky-ground"]["abv_end"] < summary["kentucky-top"]["abv_end"],
              value=[summary["kentucky-ground"]["abv_end"], summary["kentucky-top"]["abv_end"]],
              note="同じ蔵でも湿った下の階ほど度数が低い（vinepair の向き）"),
    ]
    # 🎯 RH*: その湿度の樽は度数が動かない（温度の倍率を変えても）
    flat = {}
    for f in (1.0, 3.0):
        net_f = angel_cask(rates, {"flat": {"rh": rh_star, "f": f}}, steps_per_year=52)
        r_f = net_f.run(30 * YEAR, sample_every=52)
        ser = _cask_series(r_f, "flat")
        flat[f"x{f:g}"] = max(abs(a - CASK["fill_abv"]) for _t, _v, a in ser)
        checks.append(check(f"rh-star:f{f:g}", "⚖️", flat[f"x{f:g}"] < 1e-6, value=flat[f"x{f:g}"], tol=1e-6,
                            note=f"RH* = {rh_star:.2%} の樽は 30 年たっても度数が {CASK['fill_abv']} % のまま（温度 ×{f:g}）"))
    # ⏳ 40 % を切るまで（閉じた式の二分法 vs シミュレーション）
    t40 = years_to_abv(SCOTCH_MIN_ABV, rates, SCOTLAND["rh"])
    t40_sim = next((t for t, _v, a in sc if a < SCOTCH_MIN_ABV), None)
    checks.append(check("scotland:years-to-40", "⏳", t40 is not None and t40_sim is not None and 0 <= t40_sim - t40 < 1.0,
                        value=t40_sim, expected=t40, tol=1.0,
                        note="スコッチの下限 40 % を切る年（シミュレーションは 1 年ごとに見るので 1 年以内の遅れ）"))
    # 🎲 湿度を A×10^-1 で振る: 度数の上下は RH* との大小で決まる
    sweep = []
    for item in random_scales(seed, 6, n_min=1, n_max=1, digits=2):
        rh = item["value"]
        e1, w1 = cask_state(10.0, rates, rh)
        up = abv(e1, w1) > CASK["fill_abv"]
        sweep.append({"rh": rh, "abv_10y": abv(e1, w1)})
        checks.append(check(f"sweep:rh{item['label']}", "🎲", up == (rh < rh_star), value=round(abv(e1, w1), 3),
                            note=f"湿度 {rh:.0%}: RH* より{'低い → 上がる' if rh < rh_star else '高い → 下がる'}"))
    # 🎲 刻みによらない（厳密離散化）: 1 年を A×10^N 等分（N = 1〜2）しても 20 年後が閉じた式どおり
    for item in random_scales(seed + 1, 3, n_min=-2, n_max=-1, digits=2):
        spy = max(1, int(item["value"]))
        r_d = angel_cask(rates, {"scotland": climates["scotland"]}, steps_per_year=spy).run(20 * YEAR)
        p = r_d["parts_q"]["scotland"]
        val = abv(from_q(p["ethanol"]), from_q(p["water"]))
        checks.append(close(f"dt:{spy}/yr", "🧮", val, abv(*cask_state(20.0, rates, SCOTLAND["rh"])), 1e-9,
                            note=f"1 年を {spy} 等分しても 20 年後の度数が閉じた式どおり（e^(−k·dt) の厳密離散化）"))
    # 🪵 木の吸い込み（任意）: 1 年目だけ多く減り、そのあと年 2 % に落ち着く
    soak = angel_soak(rates)
    checks += [
        close("soak:first-year", "🪵🎯", soak["loss_1y"], SOAK["first_year_loss"], 1e-6,
              note="1 年目の減り 3.75 %（🎯 puninwine の 3.5〜4 % の真ん中に capacity を合わせた）"),
        check("soak:settles", "🪵🔮", abs(soak["loss_2y"] - SCOTLAND["loss"]) < 0.002, value=soak["loss_2y"],
              expected=SCOTLAND["loss"], tol=0.002, note="2 年目は木が満ちて年 2 % 近くに戻る（🔮 予測）"),
        check("soak:conserve", "⚖️", soak["parts_error_q"] == {"ethanol": 0, "water": 0}, value=soak["parts_error_q"],
              note="樽 + 木 + 空で成分ごとに 0 ulp"),
    ]
    summary.update({"rates_per_year": rates, "rh_star": rh_star, "years_to_40": t40, "years_to_40_sim": t40_sim,
                    "liters_at_40": (mix_volume(*cask_state(t40, rates, SCOTLAND["rh"])) * 1e3) if t40 else None,
                    "rh_sweep": sweep, "soak": soak, "rh_star_flat": flat})
    return {"params": {"cask": CASK, "scotland": SCOTLAND, "climates": CLIMATES, "years": years,
                       "steps_per_year": steps_per_year,
                       "soak": SOAK, "seed": seed},
            "summary": summary, "runs": runs, "checks": checks}


SOAK = {"first_year_loss": 0.0375, "tau_years": 0.25}  # 🎯 1 年目 3.5〜4 %（puninwine）の真ん中 / 📝 木が満ちるまで約 3 か月


def angel_soak(rates: dict) -> dict:
    """木の吸い込みの capacity [mol] を 1 年目の減り 3.75 % に合わせ（二分法）、2 年目の減りを出す。"""
    climate = {"scotland": {"rh": SCOTLAND["rh"], "f": 1.0}}
    v0 = CASK["volume"]

    def run(cap: float) -> dict:
        net = angel_cask(rates, climate, steps_per_year=52, soak={"capacity": cap, "tau_years": SOAK["tau_years"]})
        res = net.run(2 * YEAR, sample_every=52)
        ser = _cask_series(res, "scotland")
        return {"loss_1y": 1.0 - ser[1][1] / v0, "loss_2y": 1.0 - ser[2][1] / ser[1][1], "res": res}

    cap = bisection(lambda c: run(c)["loss_1y"] < SOAK["first_year_loss"], 0.0, 1000.0, 60)
    out = run(cap)
    ne0, nw0 = moles_for_abv(CASK["fill_abv"], CASK["volume"])
    return {"capacity_mol": cap, "capacity_liters": cap / (ne0 + nw0) * CASK["volume"] * 1e3,
            "loss_1y": out["loss_1y"], "loss_2y": out["loss_2y"], "parts_error_q": out["res"]["parts_error_q"]}


# --- E26〜E28 お金・金・レアメタル（💴 🪙 ⛏️ 同じ配管に、流す物だけ置き換える） ----------

HUB_GROUPS = ("a", "b", "c")
# 位置は 源 → 拠点 3 つ（それぞれ腕が 1 つ）→ 外（→ 源へ戻る）。時間の単位は日、量の単位は題材ごと（兆・t・kt）。
# 💴 源 = 中央銀行（十分大きい有限の初期残高から一定の速さで出す）/ 拠点 = メガ銀行 / 腕 = 同じグループの証券 / 外 = 市中
# 🪙 源 = 金庫（総量が決まっている。出る量は残りに比例）ほかは 💴 と同じ
# ⛏️ 源 = 産出国（埋蔵量が決まっている）/ 拠点 = 製錬する加工国 / 腕 = 部品にする加工国 / 外 = 市場（戻らない）
HUB_FLOWS = {
    "money": {"icon": "💴", "unit": "兆", "names": ("central-bank", "bank", "sec", "town"),
              "source_volume": float(1 << 20), "issue": 0.5, "shares": (0.5, 0.25, 0.25), "tau_feed": None,
              "delay_feed": 0.0, "tau_out": 16.0, "tau_arm": 16.0, "delay_arm": 2.0, "tau_sat": 16.0, "tau_back": 16.0,
              "dt": 0.125},
    "gold": {"icon": "🪙", "unit": "t", "names": ("vault", "bank", "sec", "town"),
             "source_volume": 128.0, "issue": None, "shares": None, "tau_feed": (8.0, 16.0, 16.0),
             "delay_feed": 0.0, "tau_out": 16.0, "tau_arm": 16.0, "delay_arm": 2.0, "tau_sat": 16.0, "tau_back": 16.0,
             "dt": 0.125},
    "rare-metal": {"icon": "⛏️", "unit": "kt", "names": ("mine", "refiner", "maker", "market"),
                   "source_volume": 128.0, "issue": None, "shares": None, "tau_feed": (60.0, 120.0, 120.0),
                   "delay_feed": 30.0, "tau_out": None, "tau_arm": 40.0, "delay_arm": 0.0, "tau_sat": 20.0,
                   "tau_back": None, "dt": 0.25},
}


def hub_names(kind: str) -> dict:
    src, hub, arm, out = HUB_FLOWS[kind]["names"]
    return {"source": src, "hubs": [f"{hub}-{g}" for g in HUB_GROUPS], "arms": [f"{arm}-{g}" for g in HUB_GROUPS],
            "outside": out}


def hub_network(kind: str, *, dt: float | None = None, inflows: str = "uniform", **over) -> Network:
    """💴 🪙 ⛏️ の配管。源 → 拠点は issue（一定の速さ、ポンプ）か tau_feed（残りに比例、弁）。
    拠点 → 腕（tau_arm、遅れ delay_arm）、拠点 → 外（tau_out）、腕 → 外（tau_sat）、外 → 源（tau_back）。
    None の管は付けない。弁は 1 ステップに 1 − e^(−dt/τ) の割合（厳密離散化）。"""
    p = {**HUB_FLOWS[kind], **over}
    nm = hub_names(kind)
    net = Network(dt=p["dt"] if dt is None else dt, inflows=inflows)
    net.tank(nm["source"], 1.0, volume=p["source_volume"])
    for hub, arm in zip(nm["hubs"], nm["arms"]):
        net.tank(hub, 1.0).tank(arm, 1.0)
    net.tank(nm["outside"], 1.0)
    for i, (hub, arm) in enumerate(zip(nm["hubs"], nm["arms"])):
        if p["issue"] is not None:
            net.pump(nm["source"], hub, rate=p["issue"] * p["shares"][i], delay=p["delay_feed"])
        else:
            net.valve(nm["source"], hub, tau=p["tau_feed"][i], delay=p["delay_feed"])
        net.valve(hub, arm, tau=p["tau_arm"], delay=p["delay_arm"])
        if p["tau_out"] is not None:
            net.valve(hub, nm["outside"], tau=p["tau_out"])
        net.valve(arm, nm["outside"], tau=p["tau_sat"])
    if p["tau_back"] is not None:
        net.valve(nm["outside"], nm["source"], tau=p["tau_back"])
    return net


def hub_columns(kind: str, net: Network) -> list[str]:
    """断面図の列: 源 ｜（遅れの管）拠点 ｜（遅れの管）腕 … ｜ 外。"""
    nm = hub_names(kind)
    piped = {(ln["src"], ln["dst"]) for ln in net.links if ln["delay_steps"]}
    cols = [nm["source"]]
    for hub, arm in zip(nm["hubs"], nm["arms"]):
        cols += [f"{nm['source']}>{hub}"] * ((nm["source"], hub) in piped) + [hub]
        cols += [f"{hub}>{arm}"] * ((hub, arm) in piped) + [arm]
    return cols + [nm["outside"]]


def hub_caps(kind: str) -> dict:
    """断面図の縦いっぱい。🪙 ⛏️ は全部の槽を最初の総量に揃える（1 段 = 総量の 1/段数）。
    💴 は源の初期残高が大きいので、源以外を「行き渡った後にいちばん多い槽」に揃える（源は走らせた中の最大のまま）。"""
    p, nm = HUB_FLOWS[kind], hub_names(kind)
    tanks = [*nm["hubs"], *nm["arms"], nm["outside"]]
    if p["issue"] is None:
        return dict.fromkeys([nm["source"], *tanks], p["source_volume"])
    steady = hub_steady(kind)
    return dict.fromkeys(tanks, max(steady[n] for n in tanks))


def steady_tank(rate_out: float, actual: float, counted: float, dt: float | None = None) -> tuple[float, float]:
    """弁だけが出る槽の行き渡った後: (槽の量, 直結の下の槽が数える出る速さ)。出口の合計の速さ rate_out、
    実際に入る速さ actual・同じステップの中で数える速さ counted。dt=None は連続の式 actual/rate_out（💴 🪙 🫀 で共通）。"""
    if dt is None:
        return actual / rate_out, actual
    big_l = rate_out * dt
    phi = -math.expm1(-big_l)
    base = actual - counted * (1.0 - phi / big_l)
    return base * dt / phi, base


def hub_steady(kind: str, *, dt: float | None = None, **over) -> dict:
    """行き渡った後の量（閉じた式）。どの槽も「入ってくる速さ × 滞在時間」= J/K。

    dt を渡すと離散の定常（エンジンの既定 outlets="exact"・inflows="uniform" の 1 ステップをそのまま解いたもの）。
    槽に実際に入る速さ J と、同じステップの中で数える速さ I（上の槽が後から増やした分は数えない。遅れの管を通った分は
    全部数える）を分けて持つ: 槽の量 = (J − I·θ)·dt/Φ、Φ = 1 − e^(−K·dt)、θ = 1 − Φ/(K·dt)。出る量は速さの比で分かれ、
    直結の管で下の槽が数えるのは (J − I·θ) の分だけ。I = J の槽（上がポンプ・遅れの管）は J/K と厳密に一致する。
    💴 は入ってくる速さが issue で決まる（有限の初期残高が足りる範囲で循環させ、新規発行はしない）。🪙 は総量 G が決まっていて、全体の流れの速さを
    G = Σ槽 + Σ管の中 から決める（取り分は滞在時間の比で、速さにはよらない）。
    """
    p = {**HUB_FLOWS[kind], **over}
    nm = hub_names(kind)
    if p["tau_back"] is None:
        raise ValueError("no steady state without a way back to the source")

    def tank(rate_out: float, actual: float, counted: float) -> tuple[float, float]:
        return steady_tank(rate_out, actual, counted, dt)

    r_arm, r_out, r_sat, r_back = (1.0 / p[k] if p[k] else 0.0 for k in ("tau_arm", "tau_out", "tau_sat", "tau_back"))
    to_arm = r_arm / (r_arm + r_out)
    shares = p["shares"] if p["issue"] is not None else [(1.0 / t) / sum(1.0 / x for x in p["tau_feed"]) for t in p["tau_feed"]]

    def solve(flow: float, feed_counted: float) -> dict:
        """全体の流れの速さ flow（源から出る合計）と、源から拠点への直結で数えられる速さの合計 feed_counted。"""
        out, town_actual, town_counted = {}, 0.0, 0.0
        for share, hub, arm in zip(shares, nm["hubs"], nm["arms"]):
            f = flow * share
            f_counted = f if p["delay_feed"] or p["issue"] is not None else feed_counted * share
            out[hub], hub_base = tank(r_arm + r_out, f, f_counted)
            a = f * to_arm  # 腕へ行く速さ
            out[arm], arm_base = tank(r_sat, a, a if p["delay_arm"] else hub_base * to_arm)
            out[f"{hub}>{arm}"] = a * p["delay_arm"]
            out[f"{nm['source']}>{hub}"] = f * p["delay_feed"]
            town_actual += f
            town_counted += hub_base * (1.0 - to_arm) + arm_base
        out[nm["outside"]], back_base = tank(r_back, town_actual, town_counted)
        out["_back"] = back_base
        return out

    if p["issue"] is not None:
        out = solve(p["issue"], p["issue"])
        feeds = [p["issue"] * s for s in shares]
    else:  # 🪙 源も槽: 入る速さ flow、数える速さは外の槽の (J − I·θ)。全部が flow に比例するので、総量から flow を決める
        k_feed = sum(1.0 / t for t in p["tau_feed"])
        unit = solve(1.0, 1.0)
        for _ in range(8):  # 源が数える速さ → 源が直結で送る数える速さ → 外の槽 … を 1 巡りずつ（すぐ落ち着く）
            source, feed_base = tank(k_feed, 1.0, unit["_back"])
            unit = solve(1.0, feed_base)
        source, _base = tank(k_feed, 1.0, unit["_back"])
        flow = p["source_volume"] / (source + sum(v for k, v in unit.items() if k != "_back"))
        out = {k: v * flow for k, v in unit.items()}
        out[nm["source"]] = source * flow
        feeds = [flow * s for s in shares]
    out.pop("_back")
    cycle = p["delay_feed"] + 1.0 / (r_arm + r_out) + to_arm * (p["delay_arm"] + 1.0 / r_sat) + 1.0 / r_back
    return {"feeds": feeds, "cycle_time": cycle, **out}


def bateman(rates, t: float, n0: float = 1.0) -> list[float]:
    """1 次の段がつながった鎖（Bateman の式）。速さが全部違うとき、時刻 t に段 j にある量は
    N_j = n0·Π_{i<j} λ_i · Σ_{i≤j} e^(−λ_i t) / Π_{l≤j, l≠i} (λ_l − λ_i)。t < 0 は 0。"""
    out = []
    for j in range(len(rates)):
        if t < 0:
            out.append(0.0)
            continue
        total = sum(math.exp(-rates[i] * t) / math.prod(rates[m] - rates[i] for m in range(j + 1) if m != i)
                    for i in range(j + 1))
        out.append(n0 * math.prod(rates[:j]) * total)
    return out


def bateman_peak_time(k1: float, k2: float) -> float:
    """鎖の 2 段目が最大になる時刻 ln(k2/k1) / (k2 − k1)。"""
    return math.log(k2 / k1) / (k2 - k1)


def _in_transit(res: dict) -> float:
    return from_q(res["total_q"] - sum(res["volumes_q"].values()))


def _hub_checks(tag: str, kind: str, res: dict, dt: float, tol_exact: float, **over) -> tuple[list[dict], dict]:
    """行き渡った後の量を、離散の定常（きつく）と連続の式（刻みの 1 次まで）に照らす。"""
    nm = hub_names(kind)
    p = {**HUB_FLOWS[kind], **over}
    exact, smooth = hub_steady(kind, dt=dt, **over), hub_steady(kind, **over)
    k_max = max(sum(1.0 / p[k] for k in ("tau_arm", "tau_out") if p[k]), 1.0 / p["tau_sat"], 1.0 / p["tau_back"],
                sum(1.0 / t for t in p["tau_feed"]) if p["tau_feed"] else 0.0)  # いちばん速く空く槽
    names = [n for n in [nm["source"], *nm["hubs"], *nm["arms"], nm["outside"]] if n in exact]
    got = {n: from_q(res["volumes_q"][n]) for n in names}
    worst_exact = max(rel_err(got[n], exact[n]) for n in names)
    worst_smooth = max(rel_err(got[n], smooth[n]) for n in names)
    pipes = sum(v for k, v in exact.items() if ">" in k)
    checks = [
        check(f"{tag}steady:discrete", "🎯", worst_exact <= tol_exact, value=worst_exact, tol=tol_exact,
              note="どの槽も 入る速さ × 滞在時間（離散の定常: エンジンの 1 ステップをそのまま解いた式）"),
        check(f"{tag}steady:closed-form", "📐", worst_smooth <= dt * k_max, value=worst_smooth, tol=dt * k_max,
              note="連続の式（速さ 1/τ）との差（dt × いちばん速い槽の速さ まで。上がポンプ・遅れの管の槽は厳密）"),
        close(f"{tag}transit", "⏳", _in_transit(res), pipes, tol_exact, note="遅れの管の中 = 流れる速さ × 遅れ"),
        conserved(f"{tag}conserve", res),
    ]
    return checks, {"got": got, "discrete": {n: exact[n] for n in names}, "closed_form": {n: smooth[n] for n in names},
                    "in_transit": _in_transit(res), "cycle_time": smooth["cycle_time"]}


def settlement(balance: float, *, pay: float = 1.0, back: float = 0.75, overlap: bool = False,
               dt: float = 1.0 / 64) -> dict:
    """🔁 銀行 a が b へ pay、b が a へ back を 1 日で払う。overlap=False は a が午前・b が午後（順番に決済）、
    True は 1 日かけて同時（差し引き）。a の手元（balance）が足りないと、足りない分は払われない。"""
    net = Network(dt=dt)
    net.tank("bank-a", 1.0, volume=balance).tank("bank-b", 1.0, volume=4.0 * pay)
    if overlap:
        net.pump("bank-a", "bank-b", rate=pay, until=1.0)
        net.pump("bank-b", "bank-a", rate=back, until=1.0)
    else:
        net.pump("bank-a", "bank-b", rate=2.0 * pay, until=0.5)
        net.pump("bank-b", "bank-a", rate=2.0 * back, start=0.5, until=1.0)
    return net.run(1.0, trace=["bank-a"])


def e26_money_flow(seed: int = 26) -> dict:
    """💴 架空の国のお金の流れ。中央銀行が一定の速さで出し、メガ銀行 → 証券（T+2 の遅れ）→ 市中 → 中央銀行へ戻る。
    物価・金利・政策は入れない（配管だけ）。量は「出す速さ × 滞在時間」に落ち着き、出す速さに比例する。"""
    kind, p = "money", HUB_FLOWS["money"]
    nm = hub_names(kind)
    dt = p["dt"]
    net = hub_network(kind)
    res = net.run(40 * p["tau_back"], trace=[nm["arms"][0]])
    checks, steady = _hub_checks("", kind, res, dt, 1e-9)
    arm_link = next(ln for ln in res["links"] if ln["src"] == nm["hubs"][0] and ln["dst"] == nm["arms"][0])
    lag = arm_link["first_arrival"] - arm_link["first_flow"]
    outside_cb = sum(v for k, v in steady["got"].items()) + steady["in_transit"]
    smooth = hub_steady(kind)
    checks += [
        check("t+2", "📅", lag == p["delay_arm"], value=lag, expected=p["delay_arm"], tol=0,
              note="銀行を出てから証券に届くまでちょうど 2 日（遅れの管）"),
        close("stock", "💴", outside_cb, p["issue"] * smooth["cycle_time"], dt / 4.0,
              note="中央銀行の外にある量 = 出す速さ × 一巡りの滞在時間（刻みの 1 次まで）"),
    ]
    summary = {"steady": steady, "t_plus": lag, "outside_central_bank": outside_cb}
    for item in random_scales(seed, 3, n_min=-1, n_max=1):  # 出す速さ 0.1〜99.9 兆/日
        q = item["value"]
        r = hub_network(kind, issue=q).run(40 * p["tau_back"])
        tag = f"issue={item['label']}:"
        more, st = _hub_checks(tag, kind, r, dt, 1e-9, issue=q)
        ratio = st["got"][nm["outside"]] / steady["got"][nm["outside"]]
        checks += more[:1] + more[3:] + [close(f"{tag}linear", "🎲📈", ratio, q / p["issue"], 1e-9,
                                               note="出す速さを変えると源以外の定常量も同じ倍率（残高が足りる範囲）")]
    # 🔁 決済: 順番に払うなら手元は払う額だけ要る。同時（差し引き）なら 差額 + 1 ステップ分で足りる
    pay, back, sdt = 1.0, 0.75, 1.0 / 64
    full = settlement(pay)
    need_net = (pay - back) + back * sdt
    netted = settlement(need_net, overlap=True)
    checks += [
        check("gross:need", "🔁", full["links"][0]["volume_moved_q"] == to_q(pay) and min(full["trace"]["bank-a"]) == 0.0,
              value=min(full["trace"]["bank-a"]), expected=0.0, tol=0,
              note="順番に決済: 手元が払う額ちょうどなら全部払えて、途中で 0 に触れる"),
        check("net:need", "🤝", netted["links"][0]["volume_moved_q"] == to_q(pay), value=need_net,
              note="同時に決済: 手元は 差額 + 1 ステップ分で全部払える（順番なら払う額 1 が要る）"),
    ]
    for item in random_scales(seed + 1, 3, n_min=1, n_max=3):  # 足りない割合 x
        x = item["value"]
        short = settlement(pay * (1.0 - x))
        paid = short["links"][0]["volume_moved_q"]
        checks.append(check(f"short={item['label']}:unpaid", "🎲🕳️", paid == to_q(pay * (1.0 - x)),
                            value=from_q(to_q(pay) - paid), expected=x * pay, tol=0,
                            note="手元が足りないと、払われるのは手元の分ちょうど（残りは未決済）"))
        checks.append(conserved(f"short={item['label']}:conserve", short))
    summary["settlement"] = {"gross_need": pay, "net_need": need_net, "saved": pay - need_net}
    return {"params": {"flow": p, "seed": seed}, "summary": summary, "runs": {"money": res["series"]}, "checks": checks}


def e27_gold_flow(seed: int = 27) -> dict:
    """🪙 同じ配管で、流す物を総量の決まった金に置き換える。源は刷れないので、出る量は金庫の残りに比例する。
    行き渡った後の取り分は滞在時間の比だけで決まり、速く回しても総量は増えない。"""
    kind, p = "gold", HUB_FLOWS["gold"]
    nm = hub_names(kind)
    dt, total = p["dt"], p["source_volume"]
    t_end = 1024.0
    res = hub_network(kind).run(t_end)
    checks, steady = _hub_checks("", kind, res, dt, 1e-9)
    # 槽 + 管の中の合計は整数（Q）で比べる（float で足すと、足し方の違いで最後の 1 桁が動く: Python 3.11 と 3.12 以降）
    held = from_q(res["total_q"])
    share = steady["closed_form"][nm["source"]] / total
    checks += [
        check("total", "🪙", res["total_q"] == to_q(total) and held == total, value=held, expected=total, tol=0,
              note="槽 + 管の中の合計は最初の総量のまま（増えも減りもしない）"),
        close("share", "⚖️", steady["got"][nm["source"]] / total, 1.0 / (1.0 + sum(1.0 / t for t in p["tau_feed"])
                                                                         * steady["cycle_time"]), dt / 4.0,
              note="金庫に残る割合 = 1 / (1 + Σ 一巡りの滞在時間/τ_feed)"),
    ]
    summary = {"steady": steady, "vault_share": share}
    keys = ("tau_feed", "tau_out", "tau_arm", "tau_sat", "tau_back", "delay_arm")
    for item in random_scales(seed, 3, n_min=0, n_max=0):  # 全部の時定数と遅れを c 倍（1〜9.99 倍ゆっくり）
        c = item["value"]
        over = {k: (tuple(c * t for t in p[k]) if isinstance(p[k], tuple) else c * p[k]) for k in keys}
        r = hub_network(kind, **over).run(c * t_end)
        got = {n: from_q(r["volumes_q"][n]) for n in steady["got"]}
        worst = max(rel_err(got[n], steady["closed_form"][n]) for n in got)
        tag = f"slow={item['label']}:"
        checks += [
            check(f"{tag}same-shares", "🎲⚖️", worst <= dt / 4.0, value=worst, tol=dt / 4.0,
                  note="回る速さを変えても取り分は同じ（滞在時間の比だけで決まる）"),
            conserved(f"{tag}conserve", r),
        ]
    for item in random_scales(seed + 1, 3, n_min=-2, n_max=-1):  # 外 → 金庫の時定数 10〜999 日
        tau = item["value"]
        r = hub_network(kind, tau_back=tau).run(max(t_end, 24 * tau))
        more, _st = _hub_checks(f"tau_back={item['label']}:", kind, r, dt, 1e-9, tau_back=tau)
        checks += more
    return {"params": {"flow": p, "seed": seed}, "summary": summary, "runs": {"gold": res["series"]}, "checks": checks}


def _metal_run(dt: float, t_end: float, **over) -> tuple[dict, dict]:
    kind = "rare-metal"
    nm = hub_names(kind)
    res = hub_network(kind, dt=dt, **over).run(t_end, trace=[nm["source"], nm["hubs"][0], nm["arms"][0]])
    return res, nm


def e28_rare_metal(dt_list=(1.0, 0.5, 0.25), seed: int = 28) -> dict:
    """⛏️ 産出国（埋蔵量が決まっている）→ 船 30 日 → 製錬する加工国 → 部品にする加工国 → 市場。
    1 次の段の鎖なので Bateman の式が答え: 産出国は指数で減り、加工国の在庫は山を作ってから減る。"""
    p = HUB_FLOWS["rare-metal"]
    reserve, ship = p["source_volume"], p["delay_feed"]
    k_feed = [1.0 / t for t in p["tau_feed"]]
    big_k, k_arm, k_sat = sum(k_feed), 1.0 / p["tau_arm"], 1.0 / p["tau_sat"]
    t_end = 240.0
    checks, summary, errors, errors_end = [], {}, [], []

    def at(trace, t, dt):
        return trace[int(round(t / dt)) - 1]

    for dt in dt_list:
        res, nm = _metal_run(dt, t_end)
        tr = res["trace"]
        hub, arm = tr[nm["hubs"][0]], tr[nm["arms"][0]]
        n0 = reserve * k_feed[0] / big_k  # 加工国 a の取り分
        s_peak = bateman_peak_time(big_k, k_arm)
        peak_ref = bateman([big_k, k_arm], s_peak, n0)[1]
        i_peak = max(range(len(hub)), key=hub.__getitem__)
        errors.append(rel_err(hub[i_peak], peak_ref))
        old_hub = _metal_run(dt, t_end, inflows="end")[0]["trace"][nm["hubs"][0]]  # 入った量をステップの終わりに足す扱い
        errors_end.append(rel_err(max(old_hub), peak_ref))
        if dt != dt_list[-1]:
            continue
        mine_worst = max(rel_err(at(tr[nm["source"]], t, dt), reserve * math.exp(-big_k * t)) for t in (30.0, 60.0, 120.0))
        maker_worst = max(rel_err(at(arm, ship + s, dt), bateman([big_k, k_arm, k_sat], s, n0)[2])
                          for s in (20.0, 40.0, 80.0, 160.0))
        feed = next(ln for ln in res["links"] if ln["dst"] == nm["hubs"][0])
        afloat = reserve * (math.exp(-big_k * (t_end - ship)) - math.exp(-big_k * t_end))
        checks += [
            check("mine:exponential", "⛏️", mine_worst <= 1e-9, value=mine_worst, tol=1e-9,
                  note="産出国の残り = 埋蔵量 × e^(−Σt/τ_feed)（30・60・120 日。出口が 3 本でも厳密: 🚿 E29）"),
            check("ship:arrival", "🚢", feed["first_arrival"] == ship, value=feed["first_arrival"], expected=ship, tol=0,
                  note="最初の船が着くのはちょうど 30 日後"),
            close("refiner:peak-time", "⛰️", (i_peak + 1) * dt, ship + s_peak, 4 * dt / (ship + s_peak),
                  note="加工国の在庫が最大になる時刻 = 船の遅れ + ln(k2/k1)/(k2 − k1)"),
            close("refiner:peak", "📦", hub[i_peak], peak_ref, dt / p["tau_arm"],
                  note="そのときの在庫（Bateman の 2 段目、刻みの 1 次まで）"),
            check("maker:bateman", "🔩", maker_worst <= dt / p["tau_sat"], value=maker_worst, tol=dt / p["tau_sat"],
                  note="部品の加工国の在庫（Bateman の 3 段目、着いてから 20・40・80・160 日）"),
            close("afloat", "🌊", _in_transit(res), afloat, 1e-9, note="船の上にある量 = 過去 30 日に出た分"),
            conserved("conserve", res),
        ]
        summary = {"peak_day": (i_peak + 1) * dt, "peak_day_closed_form": ship + s_peak, "peak": hub[i_peak],
                   "peak_closed_form": peak_ref, "mine_left": at(tr[nm["source"]], t_end, dt),
                   "market": from_q(res["volumes_q"][nm["outside"]]), "afloat": _in_transit(res)}
        runs = {"rare-metal": res["series"]}
    orders = [errors_end[i] / errors_end[i + 1] for i in range(len(errors_end) - 1)]
    checks += [
        check("end:order", "📏", all(1.6 <= r <= 2.4 for r in orders), value=orders, expected=2.0,
              note="入った量をステップの終わりに足す扱い（inflows=\"end\"）: 刻みを半分にすると山の高さの誤差も半分（1 次）"),
        check("uniform:better", "🚰", all(a < b / 10.0 for a, b in zip(errors, errors_end)), value=errors,
              expected=errors_end, note="既定（一定の速さで入ると見なす）: どの刻みでも山の高さの誤差が 1/10 以下"),
    ]
    summary["peak_errors"] = dict(zip(map(str, dt_list), errors))
    summary["peak_errors_end"] = dict(zip(map(str, dt_list), errors_end))
    for item in random_scales(seed, 3, n_min=-1, n_max=-1):  # 製錬の時定数 10〜99.9 日
        tau = item["value"]
        if abs(1.0 / tau - big_k) < 1e-3:  # 速さが重なると Bateman の式が割れる
            continue
        dt = dt_list[-1]
        res, nm = _metal_run(dt, 480.0, tau_arm=tau)
        hub = res["trace"][nm["hubs"][0]]
        i_peak = max(range(len(hub)), key=hub.__getitem__)
        s_peak = bateman_peak_time(big_k, 1.0 / tau)
        ref = bateman([big_k, 1.0 / tau], s_peak, reserve * k_feed[0] / big_k)[1]
        tag = f"tau_arm={item['label']}:"
        checks += [
            close(f"{tag}peak", "🎲📦", hub[i_peak], ref, dt / min(tau, 1.0 / big_k), note="山の高さ（Bateman）"),
            close(f"{tag}peak-time", "🎲⛰️", (i_peak + 1) * dt, ship + s_peak, 8 * dt / (ship + s_peak), note="山の時刻"),
            conserved(f"{tag}conserve", res),
        ]
    return {"params": {"flow": p, "dt_list": list(dt_list), "seed": seed}, "summary": summary, "runs": runs,
            "checks": checks}


# --- E29 出口の数（🚿 同じ槽から弁が n 本出るとき: 別々に取るとずれる、槽ごとに厳密に解けばずれない） ----------

OUTLETS = {"volume": 128.0, "tau": 30.0, "dt": 0.25, "t_end": 120.0, "counts": (1, 2, 3, 4, 8)}
# 刻みをわざと粗くして、ずれを見える大きさに。列: 別々に取る 1・2・4・8 本 ｜ 厳密に解く 8 本
OUTLET_SCENE = {"tau": 8.0, "dt": 4.0, "t_end": 32.0, "counts": (1, 2, 4, 8, 8),
                "modes": ("independent",) * 4 + ("exact",)}


def outlet_factor(rates, dt: float) -> float:
    """弁が別々に取る槽（outlets="independent"）の、1 ステップ後に残る割合 1 − Σ(1 − e^(−k_i·dt))。
    どの弁もステップの最初の量を見て取るので、合わせると e^(−Σk_i·dt) より少し多く取る。"""
    return 1.0 - sum(1.0 - math.exp(-k * dt) for k in rates)


def outlet_error(rates, dt: float, t: float) -> float:
    """別々に取るとき、時刻 t の残りが厳密な e^(−Σk·t) から外れる割合の、刻みの 1 次の見積もり −t·dt·(K² − Σk_i²)/2。
    合計の速さ K を n 等分すると K² − Σk_i² = K²·(1 − 1/n)。1 本なら 0、本数を増やすと頭打ち。"""
    big_k = sum(rates)
    return -t * dt * (big_k * big_k - sum(k * k for k in rates)) / 2.0


def outlet_network(shares_list, *, tau: float, dt: float, volume: float = OUTLETS["volume"],
                   pump: bool = False, outlets: str = "exact", modes=None, pump_extra: float = 0.0) -> Network:
    """槽 pool-1, pool-2, … を並べ、それぞれ shares（合計 1）の割合で出口を付ける。合計の時定数はどれも tau。
    pump=True は弁の代わりに、一定の速さ volume/tau·share のポンプ（残りによらない）。
    outlets は出口の取り方（"exact" = 槽ごとに厳密、"independent" = 別々）。modes で槽ごとに変えられる。
    pump_extra > 0 なら、弁に加えて一定の速さ volume/tau·pump_extra のポンプも付ける（弁とポンプが混ざる槽）。"""
    net = Network(dt=dt, outlets=outlets)
    for i, shares in enumerate(shares_list, 1):
        net.tank(f"pool-{i}", 1.0, volume=volume, outlets=None if modes is None else modes[i - 1])
        for j, share in enumerate(shares, 1):
            net.tank(f"out-{i}-{j}", 1.0)
            if pump:
                net.pump(f"pool-{i}", f"out-{i}-{j}", rate=volume / tau * share)
            else:
                net.valve(f"pool-{i}", f"out-{i}-{j}", tau=tau / share)
        if pump_extra:
            net.tank(f"pumped-{i}", 1.0)
            net.pump(f"pool-{i}", f"pumped-{i}", rate=volume / tau * pump_extra)
    return net


def e29_outlet_count(seed: int = 29) -> dict:
    """🚿 出口の数。弁を別々に取らせる（outlets="independent"）と、n 本で合わせた減り方が刻みの 1 次ずれる（速く空く）:
    ずれは (1 − 1/n) に比例して頭打ち、刻みを半分にすると半分。エンジンの既定（outlets="exact"）は槽ごとに
    1 ステップを厳密に解くので、何本でも、ポンプが混ざっても、ずれない。命題の判定は Aoi のサイクル3。"""
    p = OUTLETS
    vol, tau, dt, t_end = p["volume"], p["tau"], p["dt"], p["t_end"]
    equal = [[1.0 / n] * n for n in p["counts"]]
    checks, errors = [], {}
    res = outlet_network(equal, tau=tau, dt=dt, outlets="independent").run(t_end)
    exact = vol * math.exp(-t_end / tau)
    steps = int(round(t_end / dt))
    for i, (n, shares) in enumerate(zip(p["counts"], equal), 1):
        got = from_q(res["volumes_q"][f"pool-{i}"])
        rates = [s / tau for s in shares]
        errors[n] = got / exact - 1.0
        checks.append(close(f"n={n}:discrete", "🎯", got, vol * outlet_factor(rates, dt) ** steps, 1e-9,
                            note="別々に取る: 残り = 最初の量 × (1 − Σ(1 − e^(−k_i·dt)))^ステップ数"))
        if n == 1:
            checks.append(check("n=1:exact", "1️⃣", abs(errors[n]) <= 1e-9, value=errors[n], expected=0.0, tol=1e-9,
                                note="出口が 1 本なら厳密（e^(−t/τ) どおり）"))
        else:
            checks.append(close(f"n={n}:first-order", "📐", errors[n], outlet_error(rates, dt, t_end), 0.02,
                                note="別々に取る: ずれ = −t·dt·K²·(1 − 1/n)/2（速く空く側）"))
            checks.append(close(f"n={n}:ratio", "🚿", errors[n] / errors[2], (1.0 - 1.0 / n) / 0.5, 0.01,
                                note="2 本のときとの比 = (1 − 1/n)/(1/2)"))
    checks.append(conserved("conserve", res))
    half = outlet_network(equal, tau=tau, dt=dt / 2, outlets="independent").run(t_end)
    order = {n: errors[n] / (from_q(half["volumes_q"][f"pool-{i}"]) / exact - 1.0)
             for i, n in enumerate(p["counts"], 1) if n > 1}
    checks.append(check("order", "📏", all(1.9 <= r <= 2.1 for r in order.values()), value=list(order.values()),
                        expected=2.0, note="別々に取る: 刻みを半分にするとずれも半分（1 次）"))
    limit = -t_end * dt / (2.0 * tau * tau)
    checks.append(check("ceiling", "🧢", all(limit * 1.02 <= e <= 0.0 for e in errors.values()), value=min(errors.values()),
                        expected=limit, note="この刻みでは −t·dt·K²/2 を超えない（粗い刻み・多い本数では超える: Aoi の O5・O16）"))
    # 🎲 等分でない出口: 効くのは本数そのものではなく K² − Σk_i²（1 本に偏るほど小さい）
    uneven = []
    for item in random_scales(seed, 3, n_min=1, n_max=2):  # 1 本目の割合 0.01〜0.999
        a = item["value"]
        shares = [a, (1.0 - a) / 2.0, (1.0 - a) / 2.0]
        uneven.append(shares)
        r = outlet_network([shares], tau=tau, dt=dt, outlets="independent").run(t_end)
        got = from_q(r["volumes_q"]["pool-1"]) / exact - 1.0
        rates = [s / tau for s in shares]
        tag = f"share={item['label']}:"
        checks += [close(f"{tag}first-order", "🎲📐", got, outlet_error(rates, dt, t_end), 0.02,
                         note="別々に取る: 等分でなくても −t·dt·(K² − Σk_i²)/2"),
                   conserved(f"{tag}conserve", r)]
    # ✅ 既定（槽ごとに 1 ステップを厳密に解く）: 何本でも、等分でなくても、粗い刻みでも、ずれない
    fixed = outlet_network(equal + uneven, tau=tau, dt=dt).run(t_end)
    worst = max(abs(from_q(v) / exact - 1.0) for k, v in fixed["volumes_q"].items() if k.startswith("pool-"))
    coarse_dt = tau / 2.0
    coarse = outlet_network(equal, tau=tau, dt=coarse_dt).run(4 * tau)
    coarse_worst = max(abs(from_q(v) / (vol * math.exp(-4.0)) - 1.0) for k, v in coarse["volumes_q"].items()
                       if k.startswith("pool-"))
    arrive = {k: from_q(v) for k, v in fixed["volumes_q"].items() if k.startswith(f"out-{len(equal) + 1}-")}
    want = [s * (vol - exact) for s in uneven[0]]
    split_worst = max(rel_err(arrive[f"out-{len(equal) + 1}-{j}"], w) for j, w in enumerate(want, 1))
    checks += [
        check("exact:any-count", "✅", worst <= 1e-9, value=worst, tol=1e-9,
              note="既定の取り方: 1〜8 本の等分と、偏った 3 本 × 3 組で、残りは e^(−t/τ) どおり"),
        check("exact:coarse", "✅🪨", coarse_worst <= 1e-9, value=coarse_worst, tol=1e-9,
              note="刻みが時定数の半分でも厳密"),
        check("exact:split", "✅🔱", split_worst <= 1e-9, value=split_worst, tol=1e-9,
              note="出た量は速さの比で分かれる（偏った 3 本）"),
        conserved("exact:conserve", fixed),
    ]
    # ✅💴 弁とポンプが混ざる槽: dV/dt = −K·V − q の解 (V0 + q/K)·e^(−K·t) − q/K
    share_q, t_mix = 1.0 / 8.0, 2.0 * tau
    q = vol / tau * share_q
    mix_exact = (vol + q * tau) * math.exp(-t_mix / tau) - q * tau
    mixed = {m: outlet_network([[0.5, 0.5]], tau=tau, dt=dt, outlets=m, pump_extra=share_q).run(t_mix)
             for m in OUTLET_MODES}
    mix_err = {m: from_q(r["volumes_q"]["pool-1"]) / mix_exact - 1.0 for m, r in mixed.items()}
    checks += [
        check("exact:with-pump", "✅💴", abs(mix_err["exact"]) <= 1e-9, value=mix_err["exact"], expected=0.0, tol=1e-9,
              note="弁 2 本 + ポンプの槽も厳密"),
        check("independent:with-pump", "💴", mix_err["independent"] < -1e-3, value=mix_err["independent"],
              note="別々に取ると、弁とポンプのあいだでもずれる"),
        conserved("exact:with-pump:conserve", mixed["exact"]),
    ]
    # 💴 一定量で出すポンプだけ: 取る量が残りによらないので、何本に分けても同じ（槽が足りているあいだ）
    pumped = outlet_network(equal, tau=tau, dt=dt, pump=True).run(tau / 2.0)
    left = [pumped["volumes_q"][f"pool-{i}"] for i in range(1, len(equal) + 1)]
    spread = max(left) - min(left)
    checks += [
        check("pump:same", "💴", spread <= max(p["counts"]) * steps, value=spread, expected=0, tol=max(p["counts"]) * steps,
              note="ポンプは何本でも残りが同じ（差は切り捨ての ulp だけ）"),
        close("pump:half", "⏱️", from_q(left[0]), vol / 2.0, 1e-12, note="τ/2 でちょうど半分（速さ = 最初の量/τ）"),
        conserved("pump:conserve", pumped),
    ]
    summary = {"independent_errors": {str(n): e for n, e in errors.items()}, "ceiling": limit,
               "tau_effective": {str(n): -dt / math.log(outlet_factor([1.0 / (n * tau)] * n, dt)) for n in p["counts"]},
               "order": {str(n): r for n, r in order.items()}, "exact_worst": worst, "exact_coarse_worst": coarse_worst,
               "with_pump": mix_err}
    return {"params": {"outlets": p, "seed": seed}, "summary": summary, "runs": {"outlets": res["series"]}, "checks": checks}


# --- E30 血液の循環（🫀 心臓 → 器官 → 心臓） ----------------------------------------

# 安静時の大人 1 人。時間は秒、量は L。体積は 血液量 × 割合（Young 2010 の表 1.2 = Guyton の数字、残り 7 % が心臓）、
# 器官へ行く割合は心拍出量の内訳（Columbia 大 Regional Circulations の表）。器官の割合は細動脈・毛細血管の 7 % を
# 流れの比で分けたもの（どの器官も滞在時間は同じ）。scene は draw の (秒, 何ステップごとに描くか)。
# stop は draw で拍動を止める秒（🛑 復路: 行き渡った後に止めると、弁で逆流しないので血が全部心臓に戻り、
# 最初の絵 = 全部心臓 に戻って一巡する。始まりの「全部心臓」も模型の仮置きなので、終わりも同じ仮置きに戻す）。
CIRCULATION = {
    "blood": 5.0, "output": 5.0 / 60.0, "dt": 0.25, "scene": (288.0, 2), "stop": 128.0,
    "flows": {"brain": 0.15, "coronary": 0.05, "kidney": 0.22, "gut-liver": 0.25, "muscle": 0.17, "skin": 0.08,
              "other": 0.08},
    "volumes": {"artery": 0.13, "organs": 0.07, "vein": 0.64, "lung": 0.09},
}


def circulation_taus() -> dict:
    """血管の時定数（ふつうの状態で決め、ポンプの速さや血液量を変えても血管は同じ）。動脈 → 器官の弁は
    1/τ = 血管の通りやすさで、比が器官へ行く割合。器官・静脈・肺は 体積 / 心拍出量（滞在時間）。"""
    p = CIRCULATION
    f, v, vol = p["output"], p["blood"], p["volumes"]
    k_art = f / (v * vol["artery"])
    return {"artery": {o: 1.0 / (k_art * s) for o, s in p["flows"].items()}, "organ": v * vol["organs"] / f,
            "vein": v * vol["vein"] / f, "lung": v * vol["lung"] / f}


def circulation_columns() -> list[str]:
    """断面図の列: 心臓 ｜ 動脈 ｜ 器官 7 つ ｜ 静脈 ｜ 肺。"""
    return ["heart", "artery", *CIRCULATION["flows"], "vein", "lung"]


def circulation_network(*, blood: float | None = None, output: float | None = None, dt: float | None = None,
                        stop: float | None = None) -> Network:
    """🫀 心臓（ポンプ）→ 動脈 → 器官 7 つ（並列）→ 静脈 → 肺 → 心臓。最初は血液が全部心臓にある。
    心臓は一定の速さ output [L/s] で送るが、空になればそのときあるだけ（戻ってきた分）しか送れない。
    stop [s] を渡すとそこで拍動を止める（🛑 弁は開いたままで逆流しないので、送った血は心臓へ戻り切る）。"""
    p, tau = CIRCULATION, circulation_taus()
    net = Network(dt=p["dt"] if dt is None else dt)
    net.tank("heart", 1.0, volume=p["blood"] if blood is None else blood)
    for name in circulation_columns()[1:]:
        net.tank(name, 1.0)
    net.pump("heart", "artery", rate=p["output"] if output is None else output, until=stop)
    for organ in p["flows"]:
        net.valve("artery", organ, tau=tau["artery"][organ])
        net.valve(organ, "vein", tau=tau["organ"])
    net.valve("vein", "lung", tau=tau["vein"])
    net.valve("lung", "heart", tau=tau["lung"])
    return net


def circulation_steady(*, blood: float | None = None, output: float | None = None, dt: float | None = None) -> dict:
    """行き渡った後（閉じた式）。どの槽も 入る速さ × 滞在時間（中心容積の原理 V = F·T）、心臓は残り。
    dt を渡すと離散の定常（💴 🪙 と同じ steady_tank）。血液が少ないかポンプが速くて心臓が空になると、送れるのは
    戻ってくる分だけ: F = 血液量 / (T_残り + dt)（空の心臓に 1 ステップで届いた量を次のステップで送る。連続なら dt → 0）。"""
    p, tau = CIRCULATION, circulation_taus()
    blood = p["blood"] if blood is None else blood
    want = p["output"] if output is None else output
    k_art = sum(1.0 / t for t in tau["artery"].values())
    unit = {}  # 送る速さ 1 L/s あたりの量
    unit["artery"], art_base = steady_tank(k_art, 1.0, 1.0, dt)  # 上がポンプなので厳密に 1/K
    vein_counted = 0.0
    for organ, t_art in tau["artery"].items():
        share = (1.0 / t_art) / k_art
        unit[organ], base = steady_tank(1.0 / tau["organ"], share, art_base * share, dt)
        vein_counted += base
    unit["vein"], vein_base = steady_tank(1.0 / tau["vein"], 1.0, vein_counted, dt)
    unit["lung"], _base = steady_tank(1.0 / tau["lung"], 1.0, vein_base, dt)
    rest = sum(unit.values())
    flow = min(want, blood / (rest + (dt or 0.0)))
    out = {k: v * flow for k, v in unit.items()}
    out["heart"] = blood - flow * rest
    return {"output": flow, "limited": flow < want, "rest_time": rest, **out}


def node_rates(res: dict) -> dict:
    """最後のステップの、槽ごとの 入る速さ・出る速さ [L/s]。"""
    rates: dict = {}
    for ln in res["links"]:
        q = res["final"][f"q:{ln['index']}"]
        rates.setdefault(ln["dst"], {"in": 0.0, "out": 0.0})["in"] += q
        rates.setdefault(ln["src"], {"in": 0.0, "out": 0.0})["out"] += q
    return rates


def _circulation_checks(tag: str, res: dict, st: dict) -> list[dict]:
    names = circulation_columns()
    worst = max(rel_err(from_q(res["volumes_q"][n]), st[n]) for n in names if st[n] > 0)
    pumped = node_rates(res)["heart"]["out"]
    return [
        check(f"{tag}steady:discrete", "🎯", worst <= 1e-9, value=worst, tol=1e-9,
              note="どの槽も 入る速さ × 滞在時間（離散の定常: エンジンの 1 ステップをそのまま解いた式）、心臓は残り"),
        close(f"{tag}output", "🫀", pumped * 60.0, st["output"] * 60.0, 1e-9,
              note="送る速さ [L/min] = min(ポンプの速さ, 血液量 / (T_残り + dt))。心臓は戻ってきた分しか送れない"),
        conserved(f"{tag}conserve", res),
    ]


def e30_circulation(seed: int = 30) -> dict:
    """🫀 血液の循環。心臓が一定の速さで送り、器官を通って静脈・肺から心臓へ戻る。器官は血液を使わずに通すだけ
    （入る = 出る、±0）なので、総量 5 L がそのまま回り続ける。どの槽も 入る速さ × 滞在時間、一周は 血液量 / 心拍出量。
    血液が減るか、ポンプだけ速くすると、心臓は戻ってきた分しか送れない（拍動・圧力・自動調節は入れていない）。"""
    p, tau = CIRCULATION, circulation_taus()
    dt, names = p["dt"], circulation_columns()
    t_end = 40.0 * tau["vein"]
    lap = p["blood"] / p["output"]
    release = t_end - 4.0 * lap
    net = circulation_network()
    net.tracer("rbc", "heart", release=release, stop="heart")  # 🩸 赤血球 1 つ（先入れ先出し、いちばん多く流れる器官を通る）
    res = net.run(t_end)
    exact, smooth = circulation_steady(dt=dt), circulation_steady()
    got = {n: from_q(res["volumes_q"][n]) for n in names}
    rates = node_rates(res)
    balance = {n: abs(r["in"] - r["out"]) / r["in"] for n, r in rates.items()}
    links = res["links"]
    moved_in = sum(ln["volume_moved_q"] for ln in links if ln["dst"] == "brain")
    moved_out = sum(ln["volume_moved_q"] for ln in links if ln["src"] == "brain")
    worst_smooth = max(rel_err(got[n], smooth[n]) for n in names)
    k_max = max(1.0 / tau["organ"], 1.0 / tau["lung"], sum(1.0 / t for t in tau["artery"].values()))
    shares = {o: rates[o]["in"] / rates["heart"]["out"] for o in p["flows"]}
    worst_share = max(abs(shares[o] - s) for o, s in p["flows"].items())
    rbc = res["tracers"][0]
    lap_got = None if rbc["arrival"] is None else rbc["arrival"] - rbc["path"][0]["t"]
    checks = _circulation_checks("", res, exact) + [
        check("steady:central-volume", "📐", worst_smooth <= dt * k_max, value=worst_smooth, tol=dt * k_max,
              note="中心容積の原理 V = F·T（連続の式）との差（dt × いちばん速い槽の速さ まで）"),
        check("brain:±0", "🧠", balance["brain"] <= 1e-9, value=balance["brain"], tol=1e-9,
              note="脳は 入る速さ = 出る速さ（使わずに通すだけ）"),
        check("brain:ledger", "🧾", moved_in - moved_out == res["volumes_q"]["brain"], value=from_q(moved_in - moved_out),
              expected=got["brain"], tol=0, note="脳に入った量 − 出た量 = いま脳にある量（Q 整数、0 ulp）"),
        check("all:±0", "🔁", max(balance.values()) <= 1e-9, value=max(balance.values()), tol=1e-9,
              note="どの槽も 入る = 出る。心臓では 戻ってくる速さ（静脈還流）= 送る速さ（心拍出量）"),
        check("shares", "🥧", worst_share <= 1e-9, value=worst_share, tol=1e-9,
              note="器官へ行く割合 = 表の割合（脳 15 %・腎臓 22 % …）。動脈 → 器官の弁の速さの比で決まる"),
        check("lap", "🩸", lap_got is not None and abs(lap_got - lap) <= 2 * dt, value=lap_got, expected=lap,
              tol=2 * dt, note="赤血球が一周する時間 = 血液量 / 心拍出量（60 s）。どの器官も滞在時間が同じなので経路によらない"),
    ]
    for item in random_scales(seed, 3, n_min=1, n_max=2):  # 🩸 失う割合 x（1〜99.9 %）
        x = item["value"]
        blood = p["blood"] * (1.0 - x)
        r = circulation_network(blood=blood).run(t_end)
        checks += _circulation_checks(f"loss={item['label']}:", r, circulation_steady(blood=blood, dt=dt))
    for item in random_scales(seed + 1, 3, n_min=0, n_max=0):  # 🫀 ポンプの速さ c 倍（1〜9.99 倍）、血管は同じ
        c = item["value"]
        r = circulation_network(output=c * p["output"]).run(t_end)
        checks += _circulation_checks(f"rate={item['label']}:", r, circulation_steady(output=c * p["output"], dt=dt))
    reserve = exact["heart"] - exact["output"] * dt
    summary = {"steady_ml": {n: 1e3 * got[n] for n in names}, "output_l_min": 60.0 * rates["heart"]["out"],
               "transit_s": {"organ": tau["organ"], "vein": tau["vein"], "lung": tau["lung"],
                             "artery": exact["artery"] / exact["output"], "heart": exact["heart"] / exact["output"]},
               "lap": lap_got, "loss_before_output_falls": reserve / p["blood"],
               "max_output_l_min": 60.0 * p["blood"] / (exact["rest_time"] + dt)}
    return {"params": {"circulation": p, "seed": seed}, "summary": summary, "checks": checks,
            "runs": {"circulation": res["series"], "blood": flow_rows(*blood_links(), "L/min")}}


# --- E31 火力発電（🔥 燃料 → 熱 → 電気。エネルギーだけを流す） ------------------------------

# 時間は秒、エネルギーは MJ（速さは MW）。効率は LHV（低位発熱量）基準で持ち、燃料は HHV（高位発熱量）で入れる:
# HHV − LHV（燃えてできた水蒸気の潜熱）は、そのまま煙突へ出る。stages は熱が溜まる槽を順に (名前, τ [s], {行き先: 割合})。
# 割合 None は発電端の効率 efficiency に合うように解く。"power" は発電端の電気で、所内率 aux で 🏠 所内 と 🔌 送電 に分ける。
# 🪨 磯子新 1 号級の USC（60 万 kW、主蒸気 600 °C、発電端 45 %・送電端 42 % LHV）、ボイラ 94 % は仮置き。
# 🌀 西名古屋 7-1 号（118.82 万 kW、1600 °C 級、発電端 63.08 % LHV）、ガスタービン単体 42.6 %（7HA.02）、
#    排熱回収ボイラの煙突 10.7 %、所内率は 1650 °C 級の 発電端 57.5 → 送電端 56 % HHV から。τ はどれも仮置き。
# table は環境省 BAT 表の 発電端 (HHV, LHV)（🪨 90〜110 万 kW 級 USC、🌀 1650 °C 級 GTCC）。
# part_load は部分負荷の出典: 🪨 630 MW 級 USC の設計の熱消費率 (負荷, kJ/kWh) と試験の (MW, 効率)、🌀 は 50 % 負荷で全負荷の
# 何割か（UKERC の典型値）、aux は 🪨 660 MW 超臨界の所内率 (負荷, 発電端に対する割合)（🌀 はデータが無いので比例のまま）。
# min_load は最低負荷の目安（Agora: 最新の石炭 25〜40 %、GTCC 30〜40 %）。startup は Agora の最新の発電所の
# 起動時間（止まってから 8 時間未満 hot・48 時間超 cold、最低負荷に届くまで、分）と負荷上げの速さ（%/分）の範囲、
# metal は厚い金属（タービンのロータ・管寄せ・ドラム）が動いているあいだの温度 °C（🪨 主蒸気、🌀 高圧の蒸気）。
# THERMAL["metal"] は 🧱 起動を決める厚い金属。tau・rate は (下限, 使う値, 上限): tau は止まった後に外気へ冷える時定数 [h]
# （Gülen: 自然に冷える HP ドラム 60〜80 h）、rate は温め直してよい速さ [K/分]（Kosman: 超臨界の蒸気タービンで許される
# 2 K/min、最適化に使った 1 K/min / Gülen: 蒸気の温度上げ 5〜10 °F/min = 2.8〜5.6 K/min）、off は止まって
# いた時間 [h]（Agora の hot < 8 h・cold > 48 h の境目）、turning は 1 段目の金属が 1,000 °F（537.8 °C）から 500 °F（260 °C）
# まで冷えるのにかかる時間 [h] の範囲（米国特許 6,898,935: ターニングで 40〜80 h）。
# curzon_ahlborn は 1975 年の表（名前, 低温 °C, 高温 °C, 実測の効率）。
THERMAL = {
    "t_cold": 25.0, "dt": 4.0, "scene": (2048.0, 8),
    "plants": {
        "coal": {"icon": "🪨", "gross": 600.0, "lhv_hhv": 26.1 / 27.3, "efficiency": 0.45, "aux": 1.0 - 42.0 / 45.0,
                 "t_hot": 600.0, "table": (0.43, 0.45), "min_load": 0.25,
                 "part_load": {"design": ((1.0, 7485.0), (0.75, 7604.0), (0.5, 7920.0)),
                               "test": ((630.0, 0.4592), (300.0, 0.4281)), "aux": ((1.0, 0.0476), (0.65, 0.0595))},
                 "startup": {"hot": (80.0, 150.0), "cold": (180.0, 360.0), "ramp": (3.0, 6.0), "metal": 600.0},
                 "stages": (("boiler", 512.0, {"steam": 0.94, "stack": 0.06}),
                                            ("steam", 64.0, {"power": None, "sea": None}))},
        "lng": {"icon": "🌀", "gross": 1188.2, "lhv_hhv": 0.90, "efficiency": 0.6308, "aux": 1.0 - 56.0 / 57.5,
                "t_hot": 1600.0, "table": (0.575, 0.64), "min_load": 0.30, "part_load": {"at_half": 0.90},
                "startup": {"hot": (30.0, 40.0), "cold": (120.0, 180.0), "ramp": (4.0, 8.0), "metal": 565.0},
                "stages": (("gt", 8.0, {"power": 0.426, "hrsg": 0.574}),
                                            ("hrsg", 256.0, {"stack": 0.107, "steam": 0.893}),
                                            ("steam", 64.0, {"power": None, "sea": None}))},
    },
    "curzon_ahlborn": (("West Thurrock", 25.0, 565.0, 0.36), ("CANDU", 25.0, 300.0, 0.30),
                       ("Larderello", 80.0, 250.0, 0.16)),
    "metal": {"tau": (60.0, 70.0, 80.0), "rate": (1.0, 2.0, 5.6), "off": {"hot": 8.0, "cold": 48.0}, "turning": (537.8, 260.0, 40.0, 80.0),
              "dt": 60.0, "scene": (52.0 * 3600.0, 24)},
}


def willans_no_load(kind: str) -> float:
    """📉 Willans の直線: 燃料 = 無負荷の分 a + (1 − a)·x（全負荷の燃料を 1、x は出力の割合）。a を出典の 2 点から決める
    （🪨 設計の熱消費率の 100 % と 50 %、🌀 50 % 負荷の効率の割合）。効率は η(x) = η₁·x / (a + (1 − a)·x)。"""
    pl = THERMAL["plants"][kind]["part_load"]
    if "design" in pl:
        (x1, hr1), *_mid, (x2, hr2) = pl["design"]
        f1, f2 = x1 * hr1, x2 * hr2  # 燃料 ∝ 出力 × 熱消費率
        return (f1 - (f2 - f1) / (x2 - x1) * x1) / f1  # 直線を出力 0 まで伸ばした燃料 / 全負荷の燃料
    return 1.0 / pl["at_half"] - 1.0  # 0.5 / (a + 0.5·(1 − a)) = at_half → a


def part_load_efficiency(kind: str, x: float) -> float:
    """発電端の効率（LHV）を出力の割合 x で（Willans）。"""
    a = willans_no_load(kind)
    return THERMAL["plants"][kind]["efficiency"] * x / (a + (1.0 - a) * x)


def aux_fixed(kind: str) -> float:
    """🏠 所内電力も 固定 + 比例 の直線: 全負荷の所内電力のうち、出力によらない分の割合 b（ファン・ポンプは負荷を下げても
    減りきらない）。出典の 2 点 (x, 所内率) を通るように解く。データが無ければ 0（比例のまま）。"""
    pts = THERMAL["plants"][kind]["part_load"].get("aux")
    if not pts:
        return 0.0
    (x1, s1), (x2, s2) = pts
    r = s1 * x1 / (s2 * x2)  # 所内電力の比（全負荷 / 部分負荷）= (b + (1 − b)·x1) / (b + (1 − b)·x2)
    return (x1 - r * x2) / (r * (1.0 - x2) - (1.0 - x1))


def aux_share(kind: str, x: float = 1.0) -> float:
    """所内率（発電端に対する割合）を出力の割合 x で: aux₁·(b/x + 1 − b)。全負荷で表の値 aux₁、下げるほど上がる。"""
    b = aux_fixed(kind)
    return THERMAL["plants"][kind]["aux"] * (b / x + 1.0 - b)


def thermal_plan(kind: str, *, marginal: bool = False) -> dict:
    """段の割合を埋める（None は発電端の効率に合わせて解く）。reach = 燃料 LHV のうちその段に入る割合。
    marginal=True は Willans の比例の燃料（無負荷の分を除いた燃料）に対する効率 η₁ / (1 − a) に合わせて解く。"""
    p = THERMAL["plants"][kind]
    target = p["efficiency"] / (1.0 - willans_no_load(kind)) if marginal else p["efficiency"]
    stages, reach, known = [], {p["stages"][0][0]: 1.0}, 0.0
    for name, tau, outs in p["stages"]:
        if None in outs.values():
            x = (target - known) / reach[name]
            outs = {"power": x, "sea": 1.0 - x}
        known += reach[name] * outs.get("power", 0.0)
        for dst, s in outs.items():
            if dst not in ("power", "stack", "sea"):
                reach[dst] = reach.get(dst, 0.0) + reach[name] * s
        stages.append((name, tau, outs))
    return {**p, "stages": tuple(stages), "reach": reach}


def thermal_columns(kind: str) -> list[str]:
    """断面図の列: ⛽ 燃料 ｜ 熱が溜まる段 … ｜ 🔌 送電 ｜ 🏠 所内 ｜ 🏭 煙突 ｜ 🌊 海。"""
    return ["fuel", *(name for name, _t, _o in THERMAL["plants"][kind]["stages"]), "grid", "aux", "stack", "sea"]


def thermal_fuel(kind: str, load: float = 1.0) -> tuple[float, float]:
    """燃やす速さ (LHV, HHV) [MW]。load は出力の割合（効率は変えない）。"""
    p = THERMAL["plants"][kind]
    lhv = load * p["gross"] / p["efficiency"]
    return lhv, lhv / p["lhv_hhv"]


def thermal_network(kind: str, *, load=1.0, slow: float = 1.0, t_fuel: float | None = None,
                    dt: float | None = None, willans: bool = False) -> Network:
    """🔥 燃料 → 段（熱が溜まる槽）→ 🔌 送電・🏠 所内・🏭 煙突・🌊 海（復水器の冷却水）。燃料は t_fuel 秒ぶんちょうど。
    slow はどの τ も何倍にするか（割合は変わらない）。load は出力の割合、または ((始まり [s], 割合), …) の切り替え。
    willans=True なら 📉 燃料も所内電力も 固定 + 比例 の直線: 無負荷の燃料（全負荷の a 倍、負荷によらない）を海へ、残りの
    比例の燃料を段へ入れ、所内電力の固定の分（全負荷の所内電力の b）は 🔌 送電から一定の速さで 🏠 所内へ（出力が 0 より
    大きい最初の切り替えから。送電が足りないあいだは取れない = 系統から受ける分は数えない）。"""
    p = THERMAL["plants"][kind]
    plan, aux = thermal_plan(kind, marginal=willans), p["aux"]
    net = Network(dt=THERMAL["dt"] if dt is None else dt)
    for name in thermal_columns(kind):
        net.tank(name, 1.0)
    t_end = THERMAL["scene"][0] if t_fuel is None else t_fuel
    a = willans_no_load(kind) if willans else 0.0
    full = p["gross"] / p["efficiency"]  # 全負荷の燃料（LHV）
    latent = 1.0 / p["lhv_hhv"] - 1.0  # HHV − LHV（水蒸気の潜熱）は燃やした LHV に比例してそのまま煙突へ
    if a > 0:
        net.pump("fuel", "sea", rate=a * full)
        net.pump("fuel", "stack", rate=a * full * latent)
    plan_load = [(0.0, load)] if isinstance(load, (int, float)) else list(load)
    for i, (start, x) in enumerate(plan_load):
        until = plan_load[i + 1][0] if i + 1 < len(plan_load) else None
        if willans:
            lhv = (1.0 - a) * x * full
            stack = lhv * latent
        else:
            # Preserve the legacy arithmetic before fixed-Q rounding.
            lhv, hhv = thermal_fuel(kind, x)
            stack = hhv - lhv
        net.pump("fuel", plan["stages"][0][0], rate=lhv, start=start, until=until)
        net.pump("fuel", "stack", rate=stack, start=start, until=until)
    steps = int(round(t_end / net.dt))
    net.nodes["fuel"]["v"] = sum(  # ポンプの量（Q）を、動いているステップの数だけ（ちょうど燃やし切る）
        ln["amount"] * max(0, min(steps, ln["until_step"] or steps) - ln["start_step"]) for ln in net.links)
    b = aux_fixed(kind) if willans else 0.0
    on = next((start for start, x in plan_load if x > 0), None)
    if b > 0 and on is not None:
        net.pump("grid", "aux", rate=aux * b * p["gross"], start=on)
        aux *= 1.0 - b  # 弁で分けるのは比例の分
    for name, tau, outs in plan["stages"]:
        for dst, s in outs.items():
            if dst == "power":
                net.valve(name, "grid", tau=slow * tau / (s * (1.0 - aux)))
                net.valve(name, "aux", tau=slow * tau / (s * aux))
            else:
                net.valve(name, dst, tau=slow * tau / s)
    return net


def thermal_steady(kind: str, *, load: float = 1.0, slow: float = 1.0, dt: float | None = None,
                   willans: bool = False) -> dict:
    """行き渡った後（閉じた式）: 段に溜まる熱 = 入る速さ × τ（💴 🪙 🫀 と同じ steady_tank）、出口の速さは割合どおり。
    lag は 送電の遅れの平均 = 電気の経路ごとの τ の和を、その経路の電気の取り分で重み付けたもの。
    willans=True は thermal_network(willans=True) と同じ（段に入るのは比例の燃料、無負荷の燃料は海へ、所内率は負荷で変わる）。"""
    p = THERMAL["plants"][kind]
    plan = thermal_plan(kind, marginal=willans)
    a = willans_no_load(kind) if willans else 0.0
    lhv = (a + (1.0 - a) * load) * p["gross"] / p["efficiency"]
    hhv, burn = lhv / p["lhv_hhv"], (1.0 - a) * load * p["gross"] / p["efficiency"]
    flows = {plan["stages"][0][0]: (burn, burn)}  # 段: (入る速さ, 同じステップの中で数える速さ)
    delay, out, power, lag = {plan["stages"][0][0]: 0.0}, {}, 0.0, 0.0
    sinks = {"stack": hhv - lhv, "sea": lhv - burn}
    for name, tau, outs in plan["stages"]:
        actual, counted = flows.pop(name)
        out[name], base = steady_tank(1.0 / (slow * tau), actual, counted, dt)
        here = delay[name] + slow * tau
        for dst, s in outs.items():
            if dst == "power":
                power += actual * s
                lag += actual * s * here
            elif dst in sinks:
                sinks[dst] += actual * s
            else:
                a, c = flows.get(dst, (0.0, 0.0))
                flows[dst] = (a + actual * s, c + base * s)
                delay[dst] = here
    aux = aux_share(kind, load) if willans else plan["aux"]
    return {"fuel_lhv": lhv, "fuel_hhv": hhv, "gross": power, "grid": power * (1.0 - aux), "aux": power * aux,
            **sinks, "lag": lag / power, "stored": out}


def carnot(t_cold: float, t_hot: float) -> float:
    return 1.0 - (t_cold + 273.15) / (t_hot + 273.15)


def curzon_ahlborn(t_cold: float, t_hot: float) -> float:
    """最大出力のときの効率 1 − √(Tc/Th)（温度差に比例する熱の出入りで出力を最大にした熱機関）。"""
    return 1.0 - math.sqrt((t_cold + 273.15) / (t_hot + 273.15))


def _thermal_checks(tag: str, kind: str, res: dict, st: dict, t_end: float) -> list[dict]:
    rates = node_rates(res)  # 最後のステップの、槽ごとの 入る・出る 速さ
    lhv = st["fuel_lhv"]
    stored = {n: from_q(res["volumes_q"][n]) for n in st["stored"]}
    worst = max(rel_err(stored[n], st["stored"][n]) for n in stored)
    eff = rates["grid"]["in"] / lhv
    grid = from_q(res["volumes_q"]["grid"])
    lag_got = t_end - grid / st["grid"]  # 送電した量 = 送電の速さ × (t − 遅れ)
    dt = res["dt"]
    return [
        check(f"{tag}stored", "🎯", worst <= 1e-9, value=worst, tol=1e-9,
              note="段に溜まる熱 = 入る速さ × τ（離散の定常: エンジンの 1 ステップをそのまま解いた式）"),
        close(f"{tag}net", "🔌", eff, st["grid"] / lhv, 1e-9,
              note="送電端の効率（LHV）= 段ごとの割合の積 × (1 − 所内率)"),
        check(f"{tag}lag", "⏳", abs(lag_got - st["lag"]) <= dt, value=lag_got, expected=st["lag"], tol=dt,
              note="送電の遅れの平均 = 電気の経路ごとの τ の和（取り分で重み付け）。燃やし始めてから溜まった熱の分だけ遅れる"),
        conserved(f"{tag}conserve", res),
    ]


def e31_thermal_power(seed: int = 31) -> dict:
    """🔥 火力発電のエネルギーの行き先。🪨 石炭（USC）と 🌀 LNG のコンバインドサイクルを、燃料 → 熱が溜まる段 → 電気・煙突・海
    と流す。効率は段ごとの割合の積、遅れは τ の和。Carnot・Curzon–Ahlborn の式とも比べる。💧 2 段目は水・蒸気の循環（水は回り、熱は水に乗る。
    🪨 はデータに無い損、🌀 は 3 圧・再熱）。📉 部分負荷は燃料も所内電力も 固定 + 比例 の直線、🌅 起動は 暖機 → 並列 → 負荷上げ。"""
    t, checks, summary = THERMAL, [], {}
    t_end = 32.0 * max(tau for p in t["plants"].values() for _n, tau, _o in p["stages"])
    for i, (kind, p) in enumerate(t["plants"].items()):
        plan = thermal_plan(kind)
        st = thermal_steady(kind, dt=t["dt"])
        res = thermal_network(kind, t_fuel=t_end).run(t_end)
        checks += _thermal_checks(f"{kind}:", kind, res, st, t_end)
        gross = st["gross"] / st["fuel_lhv"]
        hhv_net = st["grid"] / st["fuel_hhv"]
        ca, cn = curzon_ahlborn(t["t_cold"], p["t_hot"]), carnot(t["t_cold"], p["t_hot"])
        checks += [
            close(f"{kind}:gross", p["icon"], gross, p["efficiency"], 1e-12,
                  note="発電端の効率（LHV）= 出典の値（割合 None の段をこれに合うように解いた）"),
            close(f"{kind}:hhv", "🔥", p["table"][0] / p["table"][1], p["lhv_hhv"], 0.012,
                  note="表の HHV 基準 / LHV 基準 = 燃料の LHV/HHV（潜熱はそのまま煙突へ出るので、効率の比は燃料の比。表の丸め 0.5 点まで）"),
            check(f"{kind}:carnot", "📐", gross < cn, value=gross, expected=cn,
                  note="Carnot の上限 1 − Tc/Th より低い（Tc 25 °C、Th は主蒸気・タービン入口の温度）"),
            check(f"{kind}:beyond-ca", "📈", gross > ca, value=gross, expected=ca,
                  note="いまの発電所は Curzon–Ahlborn の 1 − √(Tc/Th) を超える（再熱・再生・複合で熱の出入りの損を減らした）"),
        ]
        sinks = {k: st[k] / st["fuel_hhv"] for k in ("grid", "aux", "stack", "sea")}
        summary[kind] = {"efficiency_lhv": {"gross": gross, "net": st["grid"] / st["fuel_lhv"]},
                         "efficiency_hhv": {"gross": st["gross"] / st["fuel_hhv"], "net": hhv_net},
                         "share_of_hhv": sinks, "lag_s": st["lag"], "carnot": cn, "curzon_ahlborn": ca,
                         "plan": {name: outs for name, _tau, outs in plan["stages"]}}
        for item in random_scales(seed + 2 * i, 2, n_min=1, n_max=1):  # 🎲 どの τ も c 倍速く（0.1〜0.999 倍）
            c = item["value"]
            r = thermal_network(kind, slow=c, t_fuel=t_end).run(t_end)
            checks += _thermal_checks(f"{kind}:slow={item['label']}:", kind, r,
                                      thermal_steady(kind, slow=c, dt=t["dt"]), t_end)
        for item in random_scales(seed + 2 * i + 1, 2, n_min=1, n_max=1):  # 🎲 出力を x 倍（0.1〜0.999、効率は変えない）
            x = item["value"]
            r = thermal_network(kind, load=x, t_fuel=t_end).run(t_end)
            more = _thermal_checks(f"{kind}:load={item['label']}:", kind, r, thermal_steady(kind, load=x, dt=t["dt"]),
                                   t_end)
            ratio = from_q(r["volumes_q"]["grid"]) / from_q(res["volumes_q"]["grid"])
            checks += more + [close(f"{kind}:load={item['label']}:linear", "🎲📈", ratio, x, 1e-9,
                                    note="燃料を x 倍にすると送電も x 倍（部分負荷で効率が落ちることは入れていない）")]
    # 📚 Curzon–Ahlborn（1975）の表: 実測は Carnot より低く、CA の方が近い
    for name, tc, th, seen in t["curzon_ahlborn"]:
        cn, ca = carnot(tc, th), curzon_ahlborn(tc, th)
        checks.append(check(f"ca:{name}", "📚", seen < cn and abs(seen - ca) < abs(seen - cn), value=seen, expected=ca,
                            note=f"Carnot {cn:.2f}・CA {ca:.3f}。実測は Carnot より低く、CA に近い"))
        summary.setdefault("curzon_ahlborn_table", {})[name] = {"carnot": cn, "curzon_ahlborn": ca, "observed": seen}
    part_checks, summary["part_load"] = _part_load_checks(seed, t_end)
    checks += part_checks
    start_checks, summary["startup"] = _startup_checks(seed)
    metal_checks, summary["metal"] = _metal_checks(seed)
    checks += start_checks + metal_checks
    loop_checks, summary["steam_loop"] = _steam_loop_checks(seed)
    checks += loop_checks
    return {"params": {"thermal": t, "steam_loop": STEAM_LOOP, "seed": seed}, "summary": summary, "checks": checks,
            "runs": {**{kind: thermal_network(kind).run(t["scene"][0], sample_every=t["scene"][1])["series"]
                        for kind in t["plants"]},
                     **{f"energy:{kind}": flow_rows(*energy_links(kind), "%") for kind in t["plants"]},
                     "tq": tq_rows(), "startup": startup_rows(), "sensitivity": sensitivity_rows(),
                     "rate-window": rate_window_rows()}}


def _part_load_checks(seed: int, t_end: float) -> tuple[list[dict], dict]:
    """📉 部分負荷（Willans）: 無負荷の燃料は出力によらないので、出力を下げるほど効率が下がる。"""
    checks, summary = [], {}
    for i, (kind, p) in enumerate(THERMAL["plants"].items()):
        a = willans_no_load(kind)
        full = p["gross"] / p["efficiency"]
        for item in random_scales(seed + 20 + i, 2, n_min=1, n_max=1):  # 🎲 出力の割合 0.1〜0.999
            x = item["value"]
            res = thermal_network(kind, load=x, willans=True, t_fuel=t_end).run(t_end)
            grid = node_rates(res)["grid"]
            eta = (grid["in"] - grid["out"]) / ((a + (1.0 - a) * x) * full)
            want = part_load_efficiency(kind, x) * (1.0 - aux_share(kind, x))
            checks += [close(f"{kind}:part={item['label']}:net", "📉", eta, want, 1e-9,
                             note=f"送電端の効率 = η₁·x / (a + (1 − a)·x) × (1 − 所内率 aux₁·(b/x + 1 − b))、a = {a:.4f}"
                                  f"（無負荷の燃料）、b = {aux_fixed(kind):.4f}（所内電力の固定の分）"
                                  + ("。最低負荷の目安より下は式の上だけ" if x < p["min_load"] else "")),
                       conserved(f"{kind}:part={item['label']}:conserve", res)]
        pl = p["part_load"]
        loads = (1.0, 0.75, 0.5, p["min_load"])
        curve = {x: part_load_efficiency(kind, x) for x in loads}
        summary[kind] = {"no_load": a, "efficiency": {f"{x:g}": e for x, e in curve.items()}, "aux_fixed": aux_fixed(kind),
                         "aux": {f"{x:g}": aux_share(kind, x) for x in loads},
                         "net": {f"{x:g}": curve[x] * (1.0 - aux_share(kind, x)) for x in loads}}
        if "design" in pl:
            (x1, hr1), (xm, hrm), _last = pl["design"]
            predicted = hr1 * (a + (1.0 - a) * xm) / xm  # 熱消費率 ∝ 1/効率
            (p1, e1), (p2, e2) = pl["test"]
            ratio = part_load_efficiency(kind, p2 / p1) / p["efficiency"]
            checks += [
                close(f"{kind}:part:willans", "📏", predicted, hrm, 0.01,
                      note=f"100 % と 50 % の直線から 75 % の熱消費率を出すと {predicted:.0f} kJ/kWh（設計 {hrm:.0f}）。"
                           "3 点がほぼ 1 本の直線に乗る（Willans の直線）"),
                check(f"{kind}:part:test", "🔮", abs(ratio - e2 / e1) <= 0.01, value=ratio, expected=e2 / e1, tol=0.01,
                      note=f"設計の直線から、試験の {p2:.0f} MW / {p1:.0f} MW の効率の比を出すと {ratio:.3f}（試験 {e2 / e1:.3f}）"),
            ]
            summary[kind]["predicted_75_kj_kwh"] = predicted
        else:
            checks.append(close(f"{kind}:part:half", "📚", curve[0.5] / p["efficiency"], pl["at_half"], 1e-12,
                                note="50 % 負荷で全負荷の 9 割（UKERC の典型値）。Willans の直線を当てはめた"))
    return checks, summary


def metal_warmup(kind: str, off: float, *, tau: float | None = None, rate: float | None = None) -> float:
    """🧱 暖機の時間 [s]: 止まっていた off 時間 [h] で、厚い金属は外気へ冷える（θ = θ₀·e^(−t/τ)、θ は外気より何 K 熱いか）。
    許される速さ r で温め直すと、温めているあいだも冷え続けるので（dθ/dt = r − θ/τ）θ₀ に戻るまで τ·ln((rτ − θ_off)/(rτ − θ₀))。"""
    m = THERMAL["metal"]
    tau_h = m["tau"][1] if tau is None else tau  # tau [h]・rate [K/分] を渡すと 🎚️ その値で
    t, r = 3600.0 * tau_h, (m["rate"][1] if rate is None else rate) / 60.0
    hot = THERMAL["plants"][kind]["startup"]["metal"] - THERMAL["t_cold"]
    return t * math.log((r * t - hot * math.exp(-off / tau_h)) / (r * t - hot))


def metal_network(kind: str = "coal", offs=(8.0, 48.0), *, t_end: float | None = None, dt: float | None = None) -> Network:
    """🧱 厚い金属の温度（槽の高さ = 外気より何 K 熱いか、槽ごとに off-8h など）: 動いているあいだの温度で止まり、🌬️ 外気へ
    τ で冷える。止まっていた時間 off [h] で 🔥 許される速さで温め直し（metal_warmup のあいだ）、戻ったらその温度を保つ。"""
    m = THERMAL["metal"]
    net = Network(dt=m["dt"] if dt is None else dt)
    tau, r = 3600.0 * m["tau"][1], m["rate"][1] / 60.0
    hot = THERMAL["plants"][kind]["startup"]["metal"] - THERMAL["t_cold"]
    t_end = m["scene"][0] if t_end is None else t_end
    net.tank("burner", 1.0).tank("air", 1.0)
    for off in offs:
        name, start = f"off-{off:g}h", 3600.0 * off
        back = start + metal_warmup(kind, off)
        net.tank(name, 1.0, volume=hot)
        net.valve(name, "air", tau=tau)
        net.pump("burner", name, rate=r, start=start, until=back)
        net.pump("burner", name, rate=hot / tau, start=back)
    steps = int(round(t_end / net.dt))
    net.nodes["burner"]["v"] = sum(  # 温めるのに使う分だけ（ちょうど使い切る）
        ln["amount"] * max(0, min(steps, ln["until_step"] or steps) - ln["start_step"]) for ln in net.links if ln["src"] == "burner")
    return net


def thermal_startup(kind: str, start="hot", *, at: float = 0.0, tau: float | None = None, rate: float | None = None) -> tuple:
    """🌅 起動の手順（thermal_network の load に渡す切り替え）: 暖機（無負荷の燃料を燃やし、電気は出さない）→ 最低負荷に並列 →
    1 分ごとに負荷上げの速さだけ上げて全負荷へ。起動時間 = 🧱 金属を温め直す時間（start は止まっていた時間: "hot" = 8 h・
    "cold" = 48 h、または h の数）+ 金属によらない手順の分（Agora の hot の起動時間 − 8 h の暖機、負にはしない）。
    at は出典の範囲のどこか（0 = 速い側、1 = 遅い側）。起動時間は分に丸める（切り替えが刻み dt に揃う）。
    tau・rate は 🎚️ 金属の仮定を振るとき（metal_warmup と同じ）。"""
    p, m = THERMAL["plants"][kind], THERMAL["metal"]
    lo, hi = p["startup"]["hot"]
    fast, slow = p["startup"]["ramp"][::-1]
    off = m["off"][start] if isinstance(start, str) else float(start)
    warm = {o: metal_warmup(kind, o, tau=tau, rate=rate) for o in (m["off"]["hot"], off)}
    fixed = max(0.0, 60.0 * (lo + at * (hi - lo)) - warm[m["off"]["hot"]])
    t, x, ramp = 60.0 * round((fixed + warm[off]) / 60.0), p["min_load"], (fast + at * (slow - fast)) / 100.0
    plan = [(0.0, 0.0)]
    while x < 1.0:
        plan.append((t, x))
        t, x = t + 60.0, x + ramp
    return (*plan, (t, 1.0))


def _metal_checks(seed: int) -> tuple[list[dict], dict]:
    """🧱 止まっていた時間で金属が冷えた分 ÷ 温め直してよい速さ = hot と cold の起動時間の差。"""
    m = THERMAL["metal"]
    t0, t1, lo_h, hi_h = m["turning"]
    hours = m["tau"][1] * math.log((t0 - THERMAL["t_cold"]) / (t1 - THERMAL["t_cold"]))
    checks = [check("metal:turning", "📚", lo_h <= hours <= hi_h, value=hours, expected=lo_h,
                    note=f"τ = {m['tau'][1]:g} h なら 1 段目の金属が 1,000 °F から 500 °F まで冷えるのに {hours:.0f} h。"
                         f"ターニングを回し続ける時間 {lo_h:.0f}〜{hi_h:.0f} h（米国特許 6,898,935）の中")]
    summary = {"turning_h": hours}
    for i, (kind, p) in enumerate(THERMAL["plants"].items()):
        hot = p["startup"]["metal"] - THERMAL["t_cold"]
        temp = {s: THERMAL["t_cold"] + hot * math.exp(-off / m["tau"][1]) for s, off in m["off"].items()}
        warm = {s: metal_warmup(kind, off) / 60.0 for s, off in m["off"].items()}
        lo, hi = p["startup"]["cold"]
        for item in random_scales(seed + 40 + i, 2, n_min=1, n_max=1):  # 🎲 hot の手順の分を出典の範囲のどこかから
            cold = thermal_startup(kind, "cold", at=item["value"])[1][0] / 60.0
            checks.append(check(f"{kind}:metal:cold={item['label']}", "🧱", lo <= cold <= hi, value=cold, expected=lo,
                                note=f"止まって 8 h で金属は {temp['hot']:.0f} °C、48 h で {temp['cold']:.0f} °C まで冷える。"
                                     f"{m['rate'][1]:g} K/min で温め直すと {warm['hot']:.0f} 分 / {warm['cold']:.0f} 分。hot の起動時間から"
                                     f"決めた手順の分を足した cold の起動時間 {cold:.0f} 分は Agora の {lo:.0f}〜{hi:.0f} 分の中"))
        off = m["off"]["cold"]
        dt = m["dt"]
        t_end = dt * round((3600.0 * off + metal_warmup(kind, off)) / dt)  # 温め直し終わるステップ
        net = metal_network(kind, (off,), t_end=t_end)
        res = net.run(t_end)
        tau, r = 3600.0 * m["tau"][1], m["rate"][1] / 60.0
        t_w = t_end - 3600.0 * off
        want = r * tau + (hot * math.exp(-off / m["tau"][1]) - r * tau) * math.exp(-t_w / tau)
        checks += [close(f"{kind}:metal:engine", "🧱", from_q(res["volumes_q"][f"off-{off:g}h"]), want, 1e-9,
                         note=f"金属の槽（高さ = 外気より何 K 熱いか）: {off:g} h 冷えて {t_w / 60:.0f} 分温め直した温度 = "
                              "閉じた式 rτ + (θ_止 − rτ)·e^(−t/τ)"),
                   conserved(f"{kind}:metal:conserve", res)]
        summary[kind] = {"metal_c": temp, "warm_min": warm, "cold_minus_hot_min": warm["cold"] - warm["hot"]}
    # 🎚️ τ を出典の幅で、温め直しの速さを出典の幅の中で振る: cold の起動時間がどの発電所でも（手順の分は範囲の両端 at = 0・1）
    # Agora の範囲に入る速さの窓。速いほど起動は短い（単調）ので、窓の両端は二分法で
    lo_r, rate, hi_r = m["rate"]

    def colds(tau: float, r: float) -> list[tuple[float, tuple]]:
        return [(thermal_startup(k, "cold", at=at, tau=tau, rate=r)[1][0] / 60.0, p["startup"]["cold"])
                for k, p in THERMAL["plants"].items() for at in (0.0, 1.0)]

    windows = {tau: (bisection(lambda r: any(c > hi for c, (_lo, hi) in colds(tau, r)), lo_r, hi_r, steps=48),
                     bisection(lambda r: all(c >= lo for c, (lo, _hi) in colds(tau, r)), lo_r, hi_r, steps=48))
               for tau in sweep_points(seed + 50, m["tau"])}
    checks.append(check("metal:rate-window", "🎚️", all(a <= rate <= b for a, b in windows.values()), value=rate,
                        expected=min(b for _a, b in windows.values()),
                        note="cold の起動時間が Agora の範囲に入る温め直しの速さは τ " + "・".join(
                            f"{t:.0f} h で {a:.2f}〜{b:.2f}" for t, (a, b) in windows.items())
                             + f" K/min。出典の {lo_r:g}〜{hi_r:g} K/min のうちこの狭い窓だけで、使った {rate:g} K/min（Kosman の許容）"
                               "はどの τ でも中。Kosman の最適化の 1 K/min、Gülen の蒸気の温度上げ 2.8 K/min 以上だと外れる"))
    summary["rate_window"] = {f"{t:.1f}": w for t, w in windows.items()}
    return checks, summary


def _startup_checks(seed: int) -> tuple[list[dict], dict]:
    """🌅 起動: 止まった状態から 暖機 → 並列 → 負荷上げ。燃料の比例の分は段を通って電気になり、残りは段に溜まる。"""
    checks, summary = [], {}
    for i, (kind, p) in enumerate(THERMAL["plants"].items()):
        a, full = willans_no_load(kind), p["gross"] / p["efficiency"]
        settle = 32.0 * max(tau for _n, tau, _o in p["stages"])
        st = thermal_steady(kind, willans=True, dt=THERMAL["dt"])
        plan = thermal_plan(kind, marginal=True)
        for j, start in enumerate(("hot", "cold")):
            for item in random_scales(seed + 30 + 2 * i + j, 1, n_min=1, n_max=1):  # 🎲 出典の範囲のどこか（0.1〜0.999）
                load = thermal_startup(kind, start, at=item["value"])
                t_full = load[-1][0]
                t_end = THERMAL["dt"] * math.ceil((t_full + settle) / THERMAL["dt"])
                net = thermal_network(kind, load=load, willans=True, t_fuel=t_end)
                res = net.run(t_end)
                tag = f"{kind}:startup:{start}={item['label']}:"
                ends = [*(s for s, _x in load[1:]), t_end]
                into = {plan["stages"][0][0]: sum((1.0 - a) * x * full * (e - s) for (s, x), e in zip(load, ends))}
                gross = 0.0  # 起動からの電気の合計 = 段ごとに (入った − 溜まっている) × 電気の割合
                for name, _tau, outs in plan["stages"]:
                    gone = into.pop(name) - st["stored"][name]
                    for dst, share in outs.items():
                        if dst == "power":
                            gross += gone * share
                        elif dst not in ("stack", "sea"):
                            into[dst] = into.get(dst, 0.0) + gone * share
                got = from_q(res["volumes_q"]["grid"]) + from_q(res["volumes_q"]["aux"])
                grid = node_rates(res)["grid"]
                fixed = [ln for ln in net.links if ln["src"] == "grid"]  # 所内電力の固定の分（データの無い 🌀 には無い）
                owed = sum(from_q(ln["amount"] * (res["steps"] - ln["start_step"]) - ln["flow_sum"]) for ln in fixed)
                checks += [
                    close(f"{tag}gross", "⚡", got, gross, 1e-9,
                          note=f"起動から {t_end / 3600:.2f} h の電気（送電 + 所内）= 燃やした比例の燃料が段を通った分 − 段に溜まった熱"
                               f"（暖機 {load[1][0] / 60:.0f} 分 → {p['min_load']:.0%} で並列 → {len(load) - 2} 分で全負荷）"),
                    close(f"{tag}full", "🔌", grid["in"] - grid["out"], st["grid"], 1e-9,
                          note="全負荷に上げきって落ち着いた後の送電 = 全負荷の定常（所内率は表の値に戻る）"),
                    conserved(f"{tag}conserve", res),
                ]
                if fixed:
                    checks.append(check(f"{tag}import", "🏠", owed > 0, value=owed / 3600.0, expected=0.0,
                                        note=f"並列の直後は送電が所内電力の固定の分に足りず、{owed / 3600:.2f} MWh は系統から受ける"
                                             "（ここでは数えない）"))
        best = {start: thermal_startup(kind, start) for start in ("hot", "cold")}
        summary[kind] = {start: {"sync_min": load[1][0] / 60.0, "full_min": load[-1][0] / 60.0,
                                 "warm_fuel_full_load_min": a * load[1][0] / 60.0} for start, load in best.items()}
    return checks, summary


# --- E31 の 2 段目: 水・蒸気の循環（💧 水は回り、熱は水に乗って運ばれる） ----------------------------

# 時間は秒、水は kg・熱は MJ（速さは kg/s・MW）、比エンタルピー h は kJ/kg（IAPWS、Wolfram の ThermodynamicData）。
# 🪨 電源開発 橘湾 2 号（1,050 MW、USC）: 主蒸気 3,000 t/h・再熱蒸気 2,490 t/h、主蒸気 25.0 MPa(g)/600 °C、再熱 4.28 MPa(g)/610 °C、
#    再熱器入口 4.59 MPa(g)/345 °C、節炭器入口の給水 288 °C（日立評論 2001）。給水は主蒸気の圧力、復水は 5 kPa・32 °C（仮置き）。
#    高圧の抽気 = 主蒸気 − 再熱蒸気（給水へ）。低圧側は stages 段に分け、膨張線を h で・給水の上がり幅を等分（仮置き）。
#    排気の h は出力 1,050 MW に合うように解く。ボイラ 94 % と τ [s] は仮置き（heaters・lp は段の合計）。
# 🌀 西名古屋 7-1 号の 1 段目の割合（排熱のうち蒸気へ・蒸気タービンの出力）に、H 級の 3 圧・再熱の蒸気（GE GER-3936A の表 2）を当てた
#    排熱回収ボイラ: 高圧 165 bar(g)・565 °C、再熱 23.8 bar(g)・565 °C、低圧 2.2 bar(g)・277 °C、排気 0.04064 bar(a)。h・s は蒸気表、
#    hp_s・ip_s は高圧・中圧の等エントロピー膨張の行き先（再熱・低圧の圧力。再熱器の圧損は入れない）。saturation は排気の圧力の飽和。
#    t_steam は高圧・再熱の出口 °C（📒 T–Q の記録の熱い端）。
#    drums は高圧・中圧（= 再熱の圧力）・低圧のドラムの飽和 (°C, h_f, h_g)、pinch は蒸発器の出口のガスと飽和の温度差
#    (下限, 使う値, 上限) [K]（ふつう 5〜15 K）。中圧の出口 = 低圧の入口（GE は低圧の入口を中圧の排気に揃える）、中圧の過熱器を出た蒸気 = 高圧の排気
#    （再熱器の前で混ぜても冷えない）と置く。plant は西名古屋 7 号系列（環境アセスの審査書、1 系列・大気 5 °C）: 煙突の出口 °C、
#    排熱回収ボイラ 3 基の 高圧・中圧・低圧の蒸発量 [t/h/基]、GT・ST の出力 [MW]、復水器の冷却水 [m³/s]。
# losses はデータに無い損（Thermoflow の 900 MW 級の熱勘定）: (蒸気タービンの軸の仕事 kW, 発電機 kWe, 給水ポンプを回すタービン kW,
#    排気損 kJ/kg = 8.377 BTU/lb)。発電機の効率 = 発電機 / 軸、🪨 給水ポンプ = 出力 × ポンプ / 発電機（仕事は給水に入る）。
#    排気損（羽根を出る蒸気の速さの分）は、排気の乾き度を出すときにエネルギーの釣り合いの h から引く。
# saturation は 5 kPa の飽和 (h_f, h_g, s_f, s_g)、wet は火力の低圧排気の乾き度の目安（湿り 10 % 前後）。
# cooling は (取放水の温度差 K, 海水の比熱 kJ/(kg·K), 密度 kg/m³)、nanko は関西電力 南港の既設（汽力 60 万 kW）の冷却水
# (m³/s, MW)。
STEAM_LOOP = {
    "power": 1050.0, "main": 3000.0 / 3.6, "reheat": 2490.0 / 3.6, "boiler_efficiency": 0.94, "t_hot": 600.0,
    "h": {"feed": 1270.11, "boiler": 3492.58, "hp": 3063.79, "reheater": 3694.89, "condenser": 134.09},
    "taus": {"feed": 60.0, "heaters": 60.0, "boiler": 60.0, "hp": 4.0, "reheater": 16.0, "lp": 8.0, "condenser": 120.0,
             "furnace": 4.0},
    "stages": 6, "saturation": (137.73, 2560.76, 0.47614, 8.39389), "wet": 0.90,
    "losses": (1088803.0, 1077400.0, 42894.0, 8.377 * 2.326),
    "lng": {"h": {"hp": 3474.73, "hp_s": 2932.05, "reheat": 3608.14, "ip_s": 2979.81, "lp": 3022.14}, "t_steam": 565.0,
            "s": {"hp": 6.50996, "lp": 7.58779}, "saturation": (122.54, 2554.16, 0.42619, 8.46769),
            "drums": ((350.35, 1673.79, 2561.23), (223.55, 960.05, 2801.85), (135.88, 571.51, 2728.04)),
            "pinch": (5.0, 10.0, 15.0),
            "plant": {"stack": 88.0, "steam": (246.0, 59.0, 36.0), "gt_st": (769.5, 388.5), "cooling": 25.0},
            "taus": {"gt": 4.0, "exhaust": 16.0, "hrsg": 120.0, "st": 8.0, "condenser": 120.0}},
    "cooling": (7.0, 3.99, 1025.0), "nanko": {"cooling": (26.4, 600.0)},
    "dt": 1.0, "scene": (2048.0, 16),
}


def _coal_lp(stages: int, h_ex: float, m: float, r: float, pump: float) -> dict:
    """🪨 低圧側を stages 段に: 膨張線（再熱 → 排気）を h で等分、給水の上がり幅（復水 → 高圧抽気と混ざる前）を等分。
    段 k の抽気は上から k 番目の給水加熱器へ。熱の釣り合いで抽気の量を上の加熱器から順に解き、残りが復水器へ。
    m・r は主蒸気・再熱蒸気 [kg/s]、pump は給水ポンプの仕事 [MW]（給水に入るので、加熱器で上げる分はその分少ない）。"""
    h = STEAM_LOOP["h"]
    h_top = (m * h["feed"] - (m - r) * h["hp"] - pump * 1e3) / r  # 低圧側の加熱器を出た給水（高圧抽気と混ざる前）
    steam = [h["reheater"] - k * (h["reheater"] - h_ex) / stages for k in range(1, stages + 1)]
    water = [h["condenser"] + j * (h_top - h["condenser"]) / stages for j in range(stages + 1)]
    flow, ext = r, [0.0] * stages  # ext[k-1]: 段 k の抽気（加熱器 stages+1−k へ）
    for k in range(1, stages + 1):
        j = stages + 1 - k
        ext[k - 1] = flow * (water[j] - water[j - 1]) / (steam[k - 1] - water[j - 1])
        flow -= ext[k - 1]
    through, work, h_prev = r, 0.0, h["reheater"]
    passing = []
    for k in range(stages):
        passing.append(through)
        work += through * (h_prev - steam[k]) / 1e3
        through -= ext[k]
        h_prev = steam[k]
    return {"steam": steam, "water": water, "ext": ext, "passing": passing, "condenser": flow, "work": work}


def _hrsg_plan(eta_g: float, ex_loss: float, pinch: float | None = None) -> dict:
    """🌀 3 圧・再熱: 復水器 → 高圧・中圧・低圧の節 → 高圧タービン → 再熱器（中圧の蒸気と合流）→ 中圧タービン → 低圧タービン
    （低圧の蒸気と合流）→ 復水器。📍 節ごとの蒸気の量は pinch で解く: ガスの温度は熱を渡すほど一直線に下がり（比熱一定、
    煙突の損と煙突の温度から）、高圧・中圧の蒸発器の出口でガスは飽和 + pinch。そこより熱い側で渡せる熱で、高圧（蒸発・過熱・
    再熱）、中圧（高圧の節炭・中圧の蒸発と過熱・低圧の過熱）を、残りを低圧と節炭器へ。高圧・低圧タービンの内部効率 η（共通）を、
    1 段目で解いた蒸気タービンの出力 ÷ 発電機の効率 に合うように解く（中圧は出口が低圧の入口なので、データから決まる）。"""
    q = STEAM_LOOP["lng"]
    h, (h_f, h_g, s_f, s_g), t = q["h"], q["saturation"], q["taus"]
    (t_h, f_h, _g_h), (t_i, f_i, _g_i), (t_l, f_l, g_l) = q["drums"]
    plant = THERMAL["plants"]["lng"]
    stage = {name: outs for name, _t, outs in thermal_plan("lng")["stages"]}
    fuel = plant["gross"] / plant["efficiency"]
    q_exhaust = fuel * stage["gt"]["hrsg"]
    q_steam = q_exhaust * stage["hrsg"]["steam"]
    power = q_steam * stage["steam"]["power"]  # 蒸気タービンの発電端
    c_g = q_exhaust * stage["hrsg"]["stack"] / (q["plant"]["stack"] - THERMAL["t_cold"])  # ガスの熱容量流量 MW/K
    t_exh = q["plant"]["stack"] + q_steam / c_g  # ガスタービンの出口
    pinch = q["pinch"][1] if pinch is None else pinch  # 🎚️ 渡すとその値で
    hot = [c_g * (t_exh - t_s - pinch) * 1e3 for t_s in (t_h, t_i)]  # 高圧・中圧の pinch より熱い側で渡せる熱 kW
    h_lps = h_f + (q["s"]["lp"] - s_f) / (s_g - s_f) * (h_g - h_f)  # 低圧の入口から等エントロピーで膨張した排気

    def det(m: list) -> float:
        return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))

    def cycle(eta: float) -> tuple:
        h_crh, h_ex = h["hp"] - eta * (h["hp"] - h["hp_s"]), h["lp"] - eta * (h["lp"] - h_lps)
        rows = [[h["hp"] - f_h + h["reheat"] - h_crh, h["reheat"] - h_crh, 0.0, hot[0]],  # 高圧の pinch より熱い側
                [f_h - f_i, h_crh - f_i, h["lp"] - g_l, hot[1] - hot[0]],  # 高圧と中圧の pinch のあいだ
                [h["hp"] - h_f + h["reheat"] - h_crh, h["reheat"] - h_f, h["lp"] - h_f, q_steam * 1e3]]  # 全部
        d = det([r[:3] for r in rows])
        m_hp, m_ip, m_lp = (det([[r[3] if k == j else r[k] for k in range(3)] for r in rows]) / d for j in range(3))
        m = m_hp + m_ip + m_lp
        work = (m_hp * (h["hp"] - h_crh) + (m_hp + m_ip) * (h["reheat"] - h["lp"]) + m * (h["lp"] - h_ex)) / 1e3
        return (m_hp, m_ip, m_lp), h_crh, h_ex, work

    eta = bisection(lambda e: cycle(e)[3] < power / eta_g, 0.5, 1.0)
    (m_hp, m_ip, m_lp), h_crh, h_ex, _shaft = cycle(eta)
    m = m_hp + m_ip + m_lp
    t_lp = t_exh - (hot[1] + (m_hp + m_ip) * (f_i - f_l) + m_lp * (g_l - f_l)) / c_g / 1e3  # 低圧の蒸発器の出口のガス
    gen = {k: v for k, v in (("grid", eta_g), ("sea", 1.0 - eta_g)) if v > 0}
    water = {"hp-hrsg": {"flow": m_hp, "h": h["hp"], "tau": t["hrsg"], "outs": {"hp-st": 1.0}},
             "hp-st": {"flow": m_hp, "h": h_crh, "tau": t["st"], "outs": {"reheater": 1.0}, "pull": (h["hp"], gen)},
             "ip-hrsg": {"flow": m_ip, "h": h_crh, "tau": t["hrsg"], "outs": {"reheater": 1.0}},
             "reheater": {"flow": m_hp + m_ip, "h": h["reheat"], "tau": t["hrsg"], "outs": {"ip-st": 1.0}},
             "ip-st": {"flow": m_hp + m_ip, "h": h["lp"], "tau": t["st"], "outs": {"lp-st": 1.0}, "pull": (h["reheat"], gen)},
             "lp-hrsg": {"flow": m_lp, "h": h["lp"], "tau": t["hrsg"], "outs": {"lp-st": 1.0}},
             "lp-st": {"flow": m, "h": h_ex, "tau": t["st"], "outs": {"condenser": 1.0}, "pull": (h["lp"], gen)},
             "condenser": {"flow": m, "h": h_f, "tau": t["condenser"],
                           "outs": {"hp-hrsg": m_hp / m, "ip-hrsg": m_ip / m, "lp-hrsg": m_lp / m}, "pull": (h_ex, {"sea": 1.0})}}
    duty = {n: water[n]["flow"] * (water[n]["h"] - h_in) / 1e3
            for n, h_in in (("hp-hrsg", h_f), ("ip-hrsg", h_f), ("reheater", h_crh), ("lp-hrsg", h_f))}
    steam = stage["hrsg"]["steam"]
    return {"fuel": ("gt", fuel), "cold": h_f, "power": plant["gross"], "eta": eta, "eta_g": eta_g, "h_ex": h_ex,
            "x": (h_ex - ex_loss - h_f) / (h_g - h_f), "t_exh": t_exh, "lp_pinch": t_lp - t_l, "pinch": pinch,
            "heat": {"gt": {"tau": t["gt"], "outs": {"grid": stage["gt"]["power"], "exhaust": stage["gt"]["hrsg"]}},
                     "exhaust": {"tau": t["exhaust"], "outs": {**{n: steam * d / q_steam for n, d in duty.items()},
                                                               "stack": stage["hrsg"]["stack"]}}},
            "water": water, "q_in": q_steam, "sea": q_steam - power,
            "work": {"gt": fuel * stage["gt"]["power"], "st": power}}


def r32_h(*, low: bool = False, eta: float | None = None) -> dict:
    """🫧 R32 の h [kJ/kg]: 外気 7 °C（low=True で 2 °C）。圧縮機の出口 = 入口 + 等エントロピーの上がり / 効率
    （eta を渡すと 🎚️ その値で）。凝縮器の出口は凝縮温度で決まるので外気によらない。"""
    q = HEAT_PUMP["r32"]
    st = q["low"] if low else q["h"]
    e = q["eta"][1] if eta is None else eta
    return {"evaporator": st["evaporator"], "compressor": st["evaporator"] + (st["compressor_s"] - st["evaporator"]) / e,
            "condenser": q["h"]["condenser"]}


def _refrigerant_plan(eta: float | None = None) -> dict:
    """🫧 冷媒（R32）の輪: 室外の熱交換器（蒸発器）→ 圧縮機 → 室内の熱交換器（凝縮器）→ 膨張弁（h はそのまま）→ 蒸発器。
    🔌 電気は圧縮機で、🌬️ 外気の熱は蒸発器で冷媒に入り、凝縮器で 🏠 部屋へ出る。冷媒の量 = 出す熱 / (h_圧縮機 − h_凝縮器)、
    τ の比は 1 段目のまま、合計は 封入量 / 冷媒の量（冷媒は g、熱は kJ、速さは kW）。eta は圧縮機の効率（r32_h）。"""
    p = HEAT_PUMP
    q = p["r32"]
    h = r32_h(eta=eta)
    m = p["heat"] * 1e3 / (h["compressor"] - h["condenser"])  # g/s
    stages = dict(p["stages"])
    tau = {k: v * q["charge"] * 1e3 / m / sum(stages.values()) for k, v in stages.items()}
    work, outdoor = m * (h["compressor"] - h["evaporator"]) / 1e3, m * (h["evaporator"] - h["condenser"]) / 1e3
    water = {"evaporator": {"flow": m, "h": h["evaporator"], "tau": tau["evaporator"], "outs": {"compressor": 1.0}},
             "compressor": {"flow": m, "h": h["compressor"], "tau": tau["compressor"], "outs": {"condenser": 1.0}},
             "condenser": {"flow": m, "h": h["condenser"], "tau": tau["condenser"], "outs": {"evaporator": 1.0},
                           "pull": (h["compressor"], {"room": 1.0})}}
    return {"fuel": {"power": ("compressor", work), "outdoor": ("evaporator", outdoor)}, "cold": h["condenser"],
            "power": p["heat"], "out": "room", "heat": {}, "water": water, "q_in": work + outdoor, "sea": 0.0,
            "work": {"compressor": work}, "cop": (h["compressor"] - h["condenser"]) / (h["compressor"] - h["evaporator"])}


def steam_loop_plan(kind: str = "coal", stages: int | None = None, *, losses: bool = True, scale: float = 1.0,
                    pinch: float | None = None) -> dict:
    """水の槽（通る水 flow・行き渡った後の h・τ・行き先の割合 outs・熱だけを抜く pull = (入ってくる水の h, {行き先: 割合})）と、
    熱だけの槽（τ・行き先の割合）を、データと蒸気表から閉じた式で決める。losses=True でデータに無い損（発電機・🪨 給水ポンプ・
    排気損）を入れる。scale は 🪨 の主蒸気・再熱蒸気の量の倍率、pinch は 🌀 の pinch [K]。x は排気の乾き度（排気損を引いた h で）。
    kind="r32" は ❄️ ヒートポンプの 🫧 冷媒の輪（同じ組み立てで、熱の入口が 🔌 電気と 🌬️ 外気、出口が 🏠 部屋）。"""
    if kind == "r32":
        return _refrigerant_plan()
    p = STEAM_LOOP
    st_kw, gen_kw, bfp_kw, ex_loss = p["losses"] if losses else (1.0, 1.0, 0.0, 0.0)
    eta_g = gen_kw / st_kw
    if kind == "lng":
        return _hrsg_plan(eta_g, ex_loss, pinch)
    n = p["stages"] if stages is None else stages
    h, t = p["h"], p["taus"]
    m, r = scale * p["main"], scale * p["reheat"]
    pump = p["power"] * bfp_kw / gen_kw
    shaft = p["power"] / eta_g + pump  # タービンの軸 = 発電機に入る仕事 + 給水ポンプを回す仕事
    w_hp = m * (h["boiler"] - h["hp"]) / 1e3
    h_ex = bisection(lambda x: _coal_lp(n, x, m, r, pump)["work"] > shaft - w_hp, p["saturation"][0], h["reheater"])
    lp = _coal_lp(n, h_ex, m, r, pump)
    gen = {k: v / shaft for k, v in (("grid", p["power"]), ("sea", p["power"] / eta_g - p["power"]), ("feed", pump)) if v > 0}
    q_main, q_re = m * (h["boiler"] - h["feed"]) / 1e3, r * (h["reheater"] - h["hp"]) / 1e3
    fuel = (q_main + q_re) / p["boiler_efficiency"]
    heaters = [f"heater-{j}" for j in range(1, n + 1)]
    stages_ = [f"lp-{k}" for k in range(1, n + 1)]
    water = {"feed": {"flow": m, "h": h["feed"], "tau": t["feed"], "outs": {"boiler": 1.0}},
             "boiler": {"flow": m, "h": h["boiler"], "tau": t["boiler"], "outs": {"hp": 1.0}},
             "hp": {"flow": m, "h": h["hp"], "tau": t["hp"], "outs": {"reheater": r / m, "feed": (m - r) / m},
                    "pull": (h["boiler"], gen)},
             "reheater": {"flow": r, "h": h["reheater"], "tau": t["reheater"], "outs": {stages_[0]: 1.0}}}
    for k, name in enumerate(stages_):
        through, e = lp["passing"][k], lp["ext"][k]
        down = stages_[k + 1] if k + 1 < n else "condenser"
        water[name] = {"flow": through, "h": lp["steam"][k], "tau": t["lp"] / n,
                       "outs": {down: (through - e) / through, heaters[n - 1 - k]: e / through},
                       "pull": (lp["steam"][k - 1] if k else h["reheater"], gen)}
    water["condenser"] = {"flow": lp["condenser"], "h": h["condenser"], "tau": t["condenser"], "outs": {heaters[0]: 1.0},
                          "pull": (h_ex, {"sea": 1.0})}
    feed_flow = lp["condenser"]
    for j, name in enumerate(heaters):
        feed_flow += lp["ext"][n - 1 - j]
        water[name] = {"flow": feed_flow, "h": lp["water"][j + 1], "tau": t["heaters"] / n,
                       "outs": {heaters[j + 1] if j + 1 < n else "feed": 1.0}}
    h_f, h_g = p["saturation"][:2]
    return {"fuel": ("furnace", fuel), "cold": h["condenser"], "power": p["power"], "stages": n, "h_ex": h_ex,
            "x": (h_ex - ex_loss - h_f) / (h_g - h_f), "eta_g": eta_g, "pump": pump,
            "heat": {"furnace": {"tau": t["furnace"], "outs": {"boiler": q_main / fuel, "reheater": q_re / fuel,
                                                               "stack": 1.0 - (q_main + q_re) / fuel}}},
            "water": water, "q_in": q_main + q_re, "sea": q_main + q_re - p["power"],
            "work": {"hp": w_hp, "lp": lp["work"]}}


def steam_loop_columns(kind: str = "coal") -> list[str]:
    """断面図の列（比エンタルピー）。🪨: ボイラ ｜ 高圧 ｜ 再熱器 ｜ 低圧の段 … ｜ 復水器 ｜ 給水加熱器 … ｜ 給水、
    🌀: 高圧の節 ｜ 高圧タービン ｜ 中圧の節 ｜ 再熱器 ｜ 中圧タービン ｜ 低圧の節 ｜ 低圧タービン ｜ 復水器、
    🫧: 室外の熱交換器 ｜ 圧縮機 ｜ 室内の熱交換器。"""
    water = list(steam_loop_plan(kind)["water"])
    if kind in ("lng", "r32"):
        return water
    n = STEAM_LOOP["stages"]
    return ["boiler", "hp", "reheater", *(f"lp-{k}" for k in range(1, n + 1)), "condenser",
            *(f"heater-{j}" for j in range(1, n + 1)), "feed"]


def _specific_enthalpy(parts: dict) -> float:
    """断面図で読む高さ = 比エンタルピー [MJ/kg]（水が無い槽は 0）。"""
    return parts["heat"] / parts["water"] if parts.get("water", 0.0) > 0 else 0.0


def steam_loop_network(kind: str = "coal", *, stages: int | None = None, slow: float = 1.0, t_fuel: float | None = None,
                       dt: float | None = None) -> Network:
    """💧 水（kg）と熱（MJ）を成分に持つ水の輪（🫧 kind="r32" は冷媒 g と熱 kJ）。熱だけの槽（🪨 火炉、🌀 ガスタービン・排ガス）
    や熱の入口（🫧 電気・外気）が水の槽に熱を入れ、タービンは熱だけを ⚡ に、復水器は熱だけを 🌊 海（🫧 は 🏠 部屋）に出す。水は行き渡った量から始め、どの槽も復水の h（冷えた状態）。
    熱の抜き方は、1 ステップを厳密に解いたときに行き渡った後の h が目標になるように決める（出ていく水の割合 Φ、
    熱の割合 g = Φ·(h_入 / h_出 − 1)。dt → 0 で k = k_弁·(h_入 / h_出 − 1)）。弁は 1 ステップに水の Φ·割合 を取るように
    τ を置く（熱だけを抜く槽では弁どうしを合わせて解かないため。1 本なら τ のまま）。熱の行き先が何本か（⚡・発電機の損・
    給水ポンプ）あるときも同じく、1 ステップに熱の g·割合 を取る。"""
    plan = steam_loop_plan(kind, stages)
    net = Network(dt=STEAM_LOOP["dt"] if dt is None else dt, species=("water", "heat"), volume_of=_specific_enthalpy)
    phi = {k: -math.expm1(-net.dt / (slow * w["tau"])) for k, w in plan["water"].items()}
    for name, w in plan["water"].items():
        water = w["flow"] * net.dt / phi[name]  # 行き渡った水の量（1 ステップを厳密に解いた値）
        net.tank(name, 1.0, parts={"water": water, "heat": water * plan["cold"] / 1e3})
    feeds = plan["fuel"] if isinstance(plan["fuel"], dict) else {"fuel": plan["fuel"]}  # 熱の入口: {元: (行き先, 速さ)}
    ends = {d for spec in plan["heat"].values() for d in spec["outs"]} | {
        d for w in plan["water"].values() if "pull" in w for d in w["pull"][1]}
    ends -= {*plan["water"], *plan["heat"]}
    sinks = [n for n in ("grid", "sea", "stack") if n in ends] + sorted(ends - {"grid", "sea", "stack"})
    for name in (*plan["heat"], *feeds, *sinks):
        net.tank(name, 1.0, parts={})
    steps = int(round((STEAM_LOOP["scene"][0] if t_fuel is None else t_fuel) / net.dt))
    for src, (dst, rate) in feeds.items():
        feed = net.pump(src, dst, rate=rate)
        net.nodes[src]["parts"]["heat"] = net.nodes[src]["v"] = steps * feed["amount"]
    for src, spec in plan["heat"].items():  # 熱だけの槽 → 水の槽・出口（弁だけなので合わせて厳密に解かれる）
        for dst, share in spec["outs"].items():
            net.valve(src, dst, tau=slow * spec["tau"] / share)
    for src, w in plan["water"].items():
        for dst, share in w["outs"].items():
            net.valve(src, dst, tau=-net.dt / math.log1p(-phi[src] * share))
        if "pull" in w:
            h_in, outs = w["pull"]
            g = phi[src] * (h_in / w["h"] - 1.0)
            for dst, share in outs.items():
                net.evaporate(src, dst, rates={"heat": -math.log1p(-g * share) / net.dt})
    return net


def _loop_run_checks(tag: str, plan: dict, net: Network, res: dict) -> list[dict]:
    names = list(plan["water"])
    parts = {n: {sp: from_q(q) for sp, q in res["parts_q"][n].items()} for n in names}
    last = [(ln, {sp: from_q(q) / res["dt"] for sp, q in moved.items()}) for ln, moved in zip(net.links, net.last_parts)]
    net_in = {}  # 最後のステップの 入る − 出る（槽ごと・成分ごと、速さ）
    for ln, moved in last:
        for end, sign in ((ln["src"], -1.0), (ln["dst"], 1.0)):
            for sp, q in moved.items():
                net_in.setdefault(end, {}).setdefault(sp, 0.0)
                net_in[end][sp] += sign * q
    water_balance = max(abs(net_in[n]["water"]) / plan["water"][n]["flow"] for n in names)
    start = {n: net.nodes[n]["parts"]["water"] for n in names}
    mass = max(abs(res["parts_q"][n]["water"] - start[n]) / start[n] for n in names)
    worst_h = max(rel_err(1e3 * _specific_enthalpy(parts[n]), plan["water"][n]["h"]) for n in names)
    q_in = sum(m["heat"] for ln, m in last if ln["src"] not in plan["water"] and ln["dst"] in plan["water"])
    pulled = sum(m["heat"] for ln, m in last if ln["src"] in plan["water"] and ln["dst"] not in plan["water"])
    out = plan.get("out", "grid")
    fluid = "水" if out == "grid" else "冷媒"
    return [
        check(f"{tag}water:±0", "💧", water_balance <= 1e-9, value=water_balance, tol=1e-9,
              note=f"どの槽も 入る{fluid} = 出る{fluid}（{fluid}は回るだけ）"),
        check(f"{tag}water:steady", "🎯", mass <= 1e-9, value=mass, tol=1e-9,
              note=f"槽の{fluid} = 通る{fluid} × 滞在時間（行き渡った量から始めて、熱が入っても動かない）"),
        check(f"{tag}enthalpy", "🌡️", worst_h <= 1e-9, value=worst_h, tol=1e-9,
              note="行き渡った後の比エンタルピー = 物性表・解いた値（混ざる槽は混ざった結果）"),
        close(f"{tag}power", "⚡" if out == "grid" else "🏠", net_in[out]["heat"], plan["power"], 1e-9,
              note="発電端の出力" if out == "grid" else "部屋に出す熱"),
        close(f"{tag}first-law", "🧾", q_in, pulled, 1e-9,
              note="水に入る熱 = タービンの仕事 + 復水器から海へ（第 1 法則。水が運ぶ熱は回るだけ）" if out == "grid" else
              "冷媒に入る熱（🔌 電気 + 🌬️ 外気）= 🏠 部屋へ出す熱（第 1 法則。冷媒が運ぶ熱は回るだけ）"),
        check(f"{tag}conserve", "⚖️", res["parts_error_q"]["water"] == 0 and res["parts_error_q"]["heat"] == 0,
              value=max(res["parts_error_q"].values()), expected=0, tol=0, note=f"{fluid}と熱を別々に 0 ulp"),
    ]


def _steam_loop_checks(seed: int) -> tuple[list[dict], dict]:
    p = STEAM_LOOP
    t_end = 64.0 * max(p["taus"].values())
    checks, summary = [], {}
    for kind in ("coal", "lng"):
        plan = steam_loop_plan(kind)
        net = steam_loop_network(kind, t_fuel=t_end)
        checks += _loop_run_checks(f"{kind}:loop:", plan, net, net.run(t_end))
    plan = steam_loop_plan("coal")
    ex_loss = p["losses"][3]
    # 🪨 段に分けた抽気: 段が多いほど排気は湿った側へ。まとめた 1 本だと過熱に見える
    xs = {n: steam_loop_plan("coal", n)["x"] for n in (1, 2, 4, 6, 8)}
    x0 = steam_loop_plan("coal", losses=False)["x"]
    flow = bisection(lambda c: steam_loop_plan("coal", scale=c)["x"] < p["wet"], 0.5, 1.0, steps=64)
    share = p["power"] / plan["q_in"]
    coal = thermal_plan("coal")["stages"][1][2]["power"]
    cn, ca = carnot(THERMAL["t_cold"], p["t_hot"]), curzon_ahlborn(THERMAL["t_cold"], p["t_hot"])
    d_t, c_p, rho = p["cooling"]
    cooling = plan["sea"] / (c_p * d_t) * 1e3 / rho  # m³/s
    per_100 = 100.0 * cooling / p["power"]
    nanko = 100.0 * p["nanko"]["cooling"][0] / p["nanko"]["cooling"][1]
    checks += [
        check("coal:loop:stages", "🪜", all(xs[a] > xs[b] for a, b in zip(list(xs)[:-1], list(xs)[1:])) and xs[1] > 1.0,
              value=xs[p["stages"]], expected=1.0,
              note="抽気を段に分けるほど排気は湿った側へ下がる（1 本 " + " → ".join(f"{n} 段 {x:.3f}" for n, x in xs.items())
                   + "）。まとめた 1 本だと乾き度が 1 を超え、過熱に見える"),
        check("coal:loop:wet", "💦", p["wet"] < plan["x"] < x0 and flow < 1.0, value=plan["x"], expected=p["wet"],
              note=f"データに無い損（発電機 {plan['eta_g']:.2%}・給水ポンプを回す {plan['pump']:.1f} MW・排気損 {ex_loss:.1f} kJ/kg、"
                   f"Thermoflow の熱勘定）を入れると、6 段の排気の乾き度は {x0:.3f} → {plan['x']:.3f}（湿り {1 - plan['x']:.1%}）。"
                   f"湿り 10 % までの残りは蒸気の量: 3,000 t/h は最大連続負荷のときの値で、定格の蒸発量がその {flow:.1%}"
                   f"（{p['main'] * 3.6 * flow:,.0f} t/h）なら湿り 10 %"),
        check("coal:loop:carnot", "📐", ca < share < cn, value=share, expected=cn,
              note=f"蒸気 → 電気 {share:.3f} は Carnot {cn:.3f} より低く、CA {ca:.3f} より高い"),
        check("coal:loop:coal", "📚", abs(share - coal) <= 0.03, value=share, expected=coal, tol=0.03,
              note="橘湾 2 号の流量と蒸気表から出した 蒸気 → 電気 の割合は、🪨（磯子級、再熱 620 °C）で解いた割合と 3 点以内"),
        check("coal:loop:cooling", "🌊", abs(per_100 / nanko - 1.0) <= 0.10, value=per_100, expected=nanko, tol=0.10,
              note=f"海へ出す {plan['sea']:.0f} MW を 7 °C 上げる海水 = {cooling:.1f} m³/s（100 MW あたり {per_100:.2f}）。"
                   "南港の既設 26.4 m³/s / 600 MW（補機の冷却水込み・古い汽力）と 1 割以内"),
    ]
    # 🌀 3 圧・再熱: pinch で解いた節の割合、解いた内部効率、排気の湿り（再熱なしと比べる）、蒸気の量、GT : ST、冷却水
    lng, q = steam_loop_plan("lng"), p["lng"]
    h, (h_f, h_g, s_f, s_g), nn = q["h"], q["saturation"], q["plant"]
    eta_ip = (h["reheat"] - h["lp"]) / (h["reheat"] - h["ip_s"])
    h_once = h["hp"] - lng["eta"] * (h["hp"] - (h_f + (q["s"]["hp"] - s_f) / (s_g - s_f) * (h_g - h_f)))
    wet, wet_once = 1.0 - lng["x"], 1.0 - (h_once - ex_loss - h_f) / (h_g - h_f)
    st = lng["work"]["st"]
    flows = [lng["water"][n]["flow"] for n in ("hp-hrsg", "ip-hrsg", "lp-hrsg")]
    split, nn_split = [f / sum(flows) for f in flows], [f / sum(nn["steam"]) for f in nn["steam"]]
    gt_st, nn_ratio = lng["work"]["gt"] / st, nn["gt_st"][0] / nn["gt_st"][1]
    per_mw, nn_flow = sum(flows) * 3.6 / st, 3.0 * sum(nn["steam"]) / nn["gt_st"][1]  # t/h / MW
    rise = lng["sea"] / (c_p * nn["cooling"] * rho / 1e3)  # 西名古屋の冷却水で何 °C 上がるか
    checks += [
        check("lng:loop:pinch", "📍", max(abs(a - b) for a, b in zip(split, nn_split)) <= 0.03, value=split[0],
              expected=nn_split[0], tol=0.03,
              note=f"高圧・中圧の蒸発器の出口でガス = 飽和 + {lng['pinch']:.0f} K と置いて解いた 高圧 : 中圧 : 低圧 = "
                   + " : ".join(f"{x:.1%}" for x in split) + "。西名古屋 7 号系列の蒸発量（246 : 59 : 36 t/h）は "
                   + " : ".join(f"{x:.1%}" for x in nn_split) + f"、どれも 3 点以内（ガスタービンの出口 {lng['t_exh']:.0f} °C）"),
        check("lng:loop:lp-pinch", "📍", lng["lp_pinch"] >= lng["pinch"], value=lng["lp_pinch"], expected=lng["pinch"],
              note=f"残りで作る低圧の蒸発器の出口でも、ガスは飽和より {lng['lp_pinch']:.1f} K 熱い（{lng['pinch']:.0f} K 以上なので作れる）"),
        check("lng:loop:isentropic", "📐", 0.85 <= lng["eta"] <= 0.95 and 0.85 <= eta_ip <= 0.95, value=lng["eta"],
              expected=0.9,
              note=f"蒸気タービンの出力に合わせて解いた 高圧・低圧の内部効率 {lng['eta']:.3f}、データから決まる中圧 {eta_ip:.3f}"
                   "（565 °C・23.8 bar(g) → 277 °C・2.2 bar(g)）。どちらも大型の蒸気タービンの範囲（摩擦で 1 割ほど失う）"),
        check("lng:loop:reheat", "💦", wet < 0.12 < wet_once, value=wet, expected=0.12,
              note=f"3 圧・再熱の排気の湿りは {wet:.1%}。同じ η で再熱しないと {wet_once:.1%} で、翼の寿命が縮む目安 12 % を超える"),
        check("lng:loop:flow", "🚿", abs(per_mw / nn_flow - 1.0) <= 0.10, value=per_mw, expected=nn_flow, tol=0.10,
              note=f"蒸気 {sum(flows) * 3.6:,.0f} t/h で 蒸気タービン {st:.0f} MW（1 MW あたり {per_mw:.2f} t/h）。"
                   f"西名古屋 7 号系列の {3.0 * sum(nn['steam']):,.0f} t/h / {nn['gt_st'][1]:.1f} MW（{nn_flow:.2f}）と 1 割以内"),
        check("lng:loop:gt-st", "🌀", abs(gt_st / nn_ratio - 1.0) <= 0.10, value=gt_st, expected=nn_ratio, tol=0.10,
              note=f"ガスタービン : 蒸気タービン の出力 = {gt_st:.2f}。西名古屋 7 号系列の 769.5 : 388.5 MW（{nn_ratio:.2f}）と 1 割以内"),
        check("lng:loop:cooling", "🌊", rise <= d_t, value=rise, expected=d_t,
              note=f"海へ出す {lng['sea']:.0f} MW を西名古屋の冷却水 {nn['cooling']:.0f} m³/s で受けると {rise:.1f} °C 上がる"
                   f"（取放水の温度差 {d_t:.0f} °C 以下に収まる）"),
    ]
    for item in random_scales(seed + 10, 2, n_min=1, n_max=1):  # 🎲 どの τ も c 倍（0.1〜0.999）: 状態点は変わらない
        c = item["value"]
        net = steam_loop_network("coal", slow=c, t_fuel=t_end)
        checks += _loop_run_checks(f"coal:loop:slow={item['label']}:", plan, net, net.run(t_end))
    # 🎚️ pinch を出典の幅で振る: 西名古屋と 3 点以内に収まる範囲と、どの値でも低圧が作れるか
    def nn_gap(pin: float) -> tuple[float, float]:
        pl = steam_loop_plan("lng", pinch=pin)
        fl = [pl["water"][n]["flow"] for n in ("hp-hrsg", "ip-hrsg", "lp-hrsg")]
        return max(abs(x / sum(fl) - b) for x, b in zip(fl, nn_split)), pl["lp_pinch"] - pin

    lo_p, use_p, hi_p = q["pinch"]
    gaps = {x: nn_gap(x) for x in sweep_points(seed + 11, q["pinch"])}
    edge = bisection(lambda x: nn_gap(x)[0] <= 0.03, lo_p, hi_p, steps=48)
    checks.append(check("lng:loop:pinch:sweep", "🎚️", use_p < edge and all(g[1] >= 0.0 for g in gaps.values()), value=edge,
                        expected=use_p, note=f"pinch を出典の {lo_p:g}〜{hi_p:g} K で振ると、西名古屋からのずれは "
                        + "・".join(f"{x:.1f} K で {100 * g[0]:.1f} 点" for x, g in gaps.items())
                        + f"。3 点以内に収まるのは pinch ≤ {edge:.1f} K（使った {use_p:g} K はその中）。低圧の蒸発器はどの値でも作れる"))
    # 🎚️ 🪨 ボイラ 94 % は仮置き: 1 段目の 蒸気 → 電気 = 発電端の効率 / ボイラ なので、3 点以内に収まるボイラの効率の下限を出す
    eff = thermal_plan("coal")["efficiency"]
    b_edge = bisection(lambda b: abs(share - eff / b) > 0.03, 0.5, 1.0, steps=48)
    b_use = thermal_plan("coal")["reach"]["steam"]
    checks.append(check("coal:loop:coal:boiler", "🎚️", b_edge < b_use, value=b_edge, expected=b_use,
                        note=f"1 段目の 蒸気 → 電気 は 発電端 {eff:.0%} / ボイラ。橘湾 2 号から出した {share:.3f} と 3 点以内なのは"
                             f"ボイラ ≥ {b_edge:.1%}（仮置きの {b_use:.0%} はその上。出典の幅は見つからず、境目だけ）"))
    summary["coal"] = {"h_ex": plan["h_ex"], "x_by_stages": {str(n): x for n, x in xs.items()}, "x_lossless": x0,
                       "generator": plan["eta_g"], "feed_pump_mw": plan["pump"], "flow_for_wet": flow,
                       "steam_to_power": share, "q_in_mw": plan["q_in"], "sea_mw": plan["sea"],
                       "plant_lhv": p["power"] / plan["fuel"][1], "cooling_m3_s": cooling, "cooling_per_100mw": per_100,
                       "condenser_share_of_main": plan["water"]["condenser"]["flow"] / p["main"]}
    summary["lng"] = {"steam_t_h": sum(flows) * 3.6, "split": split, "t_exh": lng["t_exh"], "lp_pinch": lng["lp_pinch"],
                      "h_ex": lng["h_ex"], "eta": lng["eta"], "eta_ip": eta_ip, "wetness": wet, "wetness_no_reheat": wet_once,
                      "t_h_per_mw": per_mw, "gt_st": gt_st, "sea_mw": lng["sea"], "cooling_rise": rise,
                      "pinch_edge_k": edge, "pinch_gap": {f"{x:.2f}": g[0] for x, g in gaps.items()}}
    summary["coal"]["boiler_edge"] = b_edge
    return checks, summary


# --- E32 ヒートポンプ（❄️ 🔥 の裏返し: 電気 + 外気の熱 → 部屋） --------------------------------

# 時間は秒、エネルギーは kJ（速さは kW）。2.8 kW 級のエアコン（Panasonic CS-TX286D2）の暖房: JIS C 9612 の条件
# （外気 7 °C・室内 20 °C）で 3.6 kW を 705 W で出す。stages は熱が溜まる槽を順に (名前, τ [s])、τ は仮置き。
# 電気は圧縮機へ、外気の熱は室外の熱交換器へ入り、どちらも冷媒に乗って室内の熱交換器から部屋へ出る。
# r32 は 🫧 冷媒の輪（同じ機種: R32 を 1.06 kg 封入）。h [kJ/kg] は R32 の物性（Wolfram）で、蒸発温度・過熱度・過冷却度は
# 建築研究所の式（蒸発 = 外気 − (0.100·凝縮 + 2.95)、過冷却 = 0.245·凝縮 − 1.72、過熱 = 4.49 − 0.036·凝縮）。圧縮機は
# 等エントロピーの出口 compressor_s と効率 eta = (下限, 使う値, 上限)（0.70 は Daikin の R32 の検討、モータの損込み。密閉形なので
# 損も冷媒に入る。幅は ±0.1: R32 のロータリーを測った例は 160 Hz で約 0.68 の山、arXiv 2407.11080）。凝縮温度 t_cnd はカタログの COP 3.6/0.705 に合う
# 39 °C（送風機の分は入れない）。outlet は 2.8 kW 級の暖房定格で測った吹き出しの空気 °C（JSRAE）。low は外気 2 °C（JIS の
# 暖房低温）の h、defrost は 5 °C 未満・湿度 8 割以上で掛ける除霜の補正（建築研究所）、jsrae は同じ 2.8 kW 級（M社）の
# カタログの暖房 (kW, W): 標準（外気 7 °C）・低温（2 °C）。
# max は 🏠 外気で落ちる最大出力: 建築研究所の式（最大出力比 = a₂(θ−7)² + a₁(θ−7) + a₀、a = b·最大能力比 + c、
# b・c = 傾き × 定格冷房能力 [kW] + 切片）と、この機種の定格冷房 2.8 kW・最大暖房 8.8 kW、カタログの低温の最大能力 (°C, kW)。
# frost °C 未満は着霜として除霜の補正を掛ける。house は JIS の APF の家（同じ部屋の冷暖房）: 暖房負荷は外気 t0 °C で 0、
# t100 °C で 0.82 × (1.25 × 定格冷房能力)（ISO 16358-2 の表 2 は 0.82 × 7 °C の暖房の全負荷能力。0.82 は 0 °C で出せる
# 能力の 7 °C に対する比。JIS は同じ部屋にするため暖房の能力を 冷房 × 1.25 に置く: 経産省 WG）。load = (0.82, 1.25)。
# tau は部屋の熱の滞在時間 [s]（仮置き、
# 行き渡った後の温度には効かない）、outs は並べる 2 部屋の外気 °C。
HEAT_PUMP = {"heat": 3.6, "power": 0.705, "t_room": 20.0, "t_out": 7.0, "dt": 1.0, "scene": (256.0, 2),
             "stages": (("evaporator", 16.0), ("compressor", 8.0), ("condenser", 32.0)),
             "r32": {"h": {"evaporator": 519.12, "compressor_s": 563.70, "condenser": 257.36}, "eta": (0.60, 0.70, 0.80),
                     "t_cnd": 39.0, "charge": 1.06, "outlet": 32.4,
                     "low": {"t_out": 2.0, "evaporator": 517.94, "compressor_s": 569.99},
                     "defrost": 0.77, "jsrae": ((3.6, 910.0), (3.5, 1320.0)), "scene": (512.0, 4)},
             "max": {"cooling": 2.8, "heat_max": 8.8, "frost": 5.0, "catalog": ((2.0, 6.4), (-7.0, 5.9), (-15.0, 5.2)),
                     "bri": ((0.000181, 0.002322, 0.003556), (-0.000184, 0.013904, 0.993431),
                             (-0.000173, -0.003980, -0.002870), (0.000367, 0.003983, 0.006376))},
             "house": {"t0": 17.0, "t100": 0.0, "load": (0.82, 1.25), "tau": 1800.0, "outs": (2.0, -15.0),
                       "scene": (4.0 * 3600.0, 30)}}


def heat_pump_max(t_out: float) -> float:
    """🏠 外気 t_out °C での最大暖房出力 [kW]: 建築研究所の式（定格暖房 × 最大出力比、外気の 2 次式）× 除霜の補正
    （frost °C 未満）。外気が冷たいほど蒸発温度が下がり、吸い込む冷媒が薄くなって 1 回転で運べる熱が減る。"""
    p = HEAT_PUMP
    m = p["max"]
    b1, b0, c1, c0 = m["bri"]
    ratio = m["heat_max"] / p["heat"]
    a2, a1, a0 = ((b1[i] * m["cooling"] + b0[i]) * ratio + c1[i] * m["cooling"] + c0[i] for i in range(3))
    x = t_out - p["t_out"]
    return (a2 * x * x + a1 * x + a0) * p["heat"] * (p["r32"]["defrost"] if t_out < m["frost"] else 1.0)


def house_load(t_out: float, *, iso: bool = False) -> float:
    """🏠 JIS の家の暖房負荷 [kW]: 外気 t0 で 0、t100 で 0.82 × 1.25 × 定格冷房能力 の直線（室温 20 °C を保つのに要る熱。
    t0 が 20 °C より低いのは、人・家電・日射の熱で 3 K ぶんまかなえるから）。iso=True は ISO 16358-2 のまま
    （t100 で 0.82 × 定格暖房能力）。"""
    p = HEAT_PUMP
    h = p["house"]
    ratio, heat_per_cool = h["load"]
    full = p["heat"] if iso else heat_per_cool * p["max"]["cooling"]
    return ratio * full * (h["t0"] - t_out) / (h["t0"] - h["t100"])


def house_balance(*, iso: bool = False) -> float:
    """🏠 バランス点 °C: 最大暖房出力 = 家の暖房負荷 になる外気。これより冷たいと 1 台では 20 °C を保てない。"""
    return bisection(lambda t: heat_pump_max(t) < house_load(t, iso=iso), -40.0, HEAT_PUMP["max"]["frost"])


def house_network(t_outs=None, *, slow: float = 1.0, t_end: float | None = None, dt: float | None = None) -> Network:
    """🏠 部屋ごと（room+2 など、外気 °C）の熱: ヒートポンプは要る分（家の暖房負荷）を出し、最大出力が足りなければ最大のまま。
    人・家電・日射の分 UA·(室温 − t0) も入る。部屋の熱は壁と換気から外気へ UA·(室温 − 外気) で逃げる（τ の弁）。
    槽の面積 = 熱容量 UA·τ [kJ/K] なので、高さ = 室温 − 外気 [K]。外気の温度の部屋から始める。slow は τ の倍率。"""
    p = HEAT_PUMP
    h = p["house"]
    net = Network(dt=p["dt"] if dt is None else dt)
    ua = house_load(h["t100"]) / (h["t0"] - h["t100"])  # kW/K
    t_end = h["scene"][0] if t_end is None else t_end
    net.tank("heat", 1.0).tank("outdoor", 1.0)
    for t in h["outs"] if t_outs is None else t_outs:
        name = f"room{t:+g}"
        net.tank(name, ua * slow * h["tau"])
        net.pump("heat", name, rate=min(house_load(t), heat_pump_max(t)) + ua * (p["t_room"] - h["t0"]))
        net.valve(name, "outdoor", tau=slow * h["tau"])
    steps = int(round(t_end / net.dt))
    net.nodes["heat"]["v"] = sum(ln["amount"] for ln in net.links if ln["src"] == "heat") * steps
    return net


def heat_pump_columns() -> list[str]:
    """断面図の列: 🔌 電気 ｜ 🌬️ 外気 ｜ 室外の熱交換器 ｜ 圧縮機 ｜ 室内の熱交換器 ｜ 🏠 部屋。"""
    return ["power", "outdoor", *(name for name, _t in HEAT_PUMP["stages"]), "room"]


def heat_pump_network(*, slow: float = 1.0, t_run: float | None = None, dt: float | None = None) -> Network:
    """🔌 電気 → 圧縮機、🌬️ 外気 → 室外の熱交換器 → 圧縮機 → 室内の熱交換器 → 🏠 部屋。電気と外気は t_run 秒ぶんちょうど。"""
    p = HEAT_PUMP
    net = Network(dt=p["dt"] if dt is None else dt)
    for name in heat_pump_columns():
        net.tank(name, 1.0)
    (evap, t_e), (comp, t_c), (cond, t_h) = p["stages"]
    steps = int(round((p["scene"][0] if t_run is None else t_run) / net.dt))
    for src, dst, rate in (("power", comp, p["power"]), ("outdoor", evap, p["heat"] - p["power"])):
        ln = net.pump(src, dst, rate=rate)
        net.nodes[src]["v"] = steps * ln["amount"]  # ちょうど steps 回分
    net.valve(evap, comp, tau=slow * t_e)
    net.valve(comp, cond, tau=slow * t_c)
    net.valve(cond, "room", tau=slow * t_h)
    return net


def heat_pump_steady(*, slow: float = 1.0, dt: float | None = None) -> dict:
    """行き渡った後: 段の熱 = 入る速さ × τ（steady_tank）、部屋に入る速さ = 電気 + 外気。lag は部屋に届く熱の遅れの平均。"""
    p = HEAT_PUMP
    (_e, t_e), (_c, t_c), (_h, t_h) = p["stages"]
    w, q_c = p["power"], p["heat"] - p["power"]
    evap, e_base = steady_tank(1.0 / (slow * t_e), q_c, q_c, dt)
    comp, c_base = steady_tank(1.0 / (slow * t_c), w + q_c, w + e_base, dt)
    cond, _base = steady_tank(1.0 / (slow * t_h), w + q_c, c_base, dt)
    lag = slow * (t_c + t_h + t_e * q_c / (w + q_c))
    return {"room": w + q_c, "cop": (w + q_c) / w, "lag": lag,
            "stored": {"evaporator": evap, "compressor": comp, "condenser": cond}}


def carnot_cop(t_cold: float, t_hot: float) -> float:
    """暖房の COP の上限 Th/(Th − Tc) = 1 / (Carnot の効率)。"""
    return 1.0 / carnot(t_cold, t_hot)


def e32_heat_pump(seed: int = 32) -> dict:
    """❄️ ヒートポンプは 🔥 の裏返し: 電気 W で外気から Q_c を汲み、部屋に W + Q_c を出す（COP = (W + Q_c)/W > 1）。
    上限は Carnot の COP = 1/η_Carnot。🔥 の発電所の電気で回すと、燃料 1 から部屋に η·COP の熱が届く。
    🫧 2 段目は冷媒（R32）の輪（冷媒は回り、熱は冷媒に乗る）。外気が冷たいと蒸発温度が下がり、COP も下がる。
    🏠 3 段目は外気で落ちる最大出力と部屋から逃げる熱: 交わる外気（バランス点）より冷たいと 1 台では 20 °C を保てない。"""
    p, checks = HEAT_PUMP, []
    t_end = 32.0 * max(t for _n, t in p["stages"])
    st = heat_pump_steady(dt=p["dt"])

    def run_checks(tag: str, res: dict, st: dict) -> list[dict]:
        rates = node_rates(res)
        stored = {n: from_q(res["volumes_q"][n]) for n in st["stored"]}
        worst = max(rel_err(stored[n], st["stored"][n]) for n in stored)
        room_in, w = rates["room"]["in"], rates["power"]["out"]
        lag_got = t_end - from_q(res["volumes_q"]["room"]) / st["room"]
        return [
            check(f"{tag}stored", "🎯", worst <= 1e-9, value=worst, tol=1e-9,
                  note="段の熱 = 入る速さ × τ（離散の定常: エンジンの 1 ステップをそのまま解いた式）"),
            close(f"{tag}cop", "❄️", room_in / w, st["cop"], 1e-9, note="COP = 部屋に出る熱 / 電気（3.6 kW / 705 W）"),
            close(f"{tag}balance", "🧾", room_in, w + rates["outdoor"]["out"], 1e-12,
                  note="部屋に出る熱 = 電気 + 外気から汲んだ熱（エネルギーは増えない。外から運んでくるだけ）"),
            check(f"{tag}lag", "⏳", abs(lag_got - st["lag"]) <= res["dt"], value=lag_got, expected=st["lag"],
                  tol=res["dt"], note="部屋に届く熱の遅れの平均 = 経路ごとの τ の和（電気と外気の取り分で重み付け）"),
            conserved(f"{tag}conserve", res),
        ]

    res = heat_pump_network(t_run=t_end).run(t_end)
    checks += run_checks("", res, st)
    limit = carnot_cop(p["t_out"], p["t_room"])
    checks += [
        check("carnot", "📐", st["cop"] < limit, value=st["cop"], expected=limit,
              note="Carnot の上限 Th/(Th − Tc)（外気 7 °C → 室内 20 °C）より低い"),
    ]
    for item in random_scales(seed, 4, n_min=-2, n_max=-1):  # 🎲 温度を振る: Tc 10〜999 °C、Th = Tc + 1〜99.9 K
        tc = item["value"]
        th = tc + item["A"] * 10 ** (-item["N"] - 1)
        cop, eta = carnot_cop(tc, th), carnot(tc, th)
        checks.append(close(f"reverse:{item['label']}", "🔁", cop * eta, 1.0, 1e-12,
                            note="Carnot のヒートポンプは Carnot の熱機関を逆に回したもの: COP × η = 1。冷房の COP は 1 小さい"))
    for item in random_scales(seed + 1, 2, n_min=1, n_max=1):  # 🎲 どの τ も c 倍（0.1〜0.999）
        c = item["value"]
        r = heat_pump_network(slow=c, t_run=t_end).run(t_end)
        checks += run_checks(f"slow={item['label']}:", r, heat_pump_steady(slow=c, dt=p["dt"]))
    # 🔥 → ❄️ 発電所の電気で回す: 燃料（HHV）1 から部屋に届く熱 = 送電端の効率（HHV）× COP。直接燃やすと 1 を超えない
    chain = {}
    for kind, plant in THERMAL["plants"].items():
        ts = thermal_steady(kind)
        eta = ts["grid"] / ts["fuel_hhv"]
        chain[kind] = {"eta_hhv": eta, "heat_per_fuel": eta * st["cop"], "break_even_cop": 1.0 / eta}
        checks.append(check(f"chain:{kind}", plant["icon"] + "❄️", eta * st["cop"] > 1.0, value=eta * st["cop"],
                            expected=1.0, note=f"燃料 1 → 電気 {eta:.3f} → 部屋に {eta * st['cop']:.2f}。COP が {1.0 / eta:.2f} を超えれば"
                                               "、燃料を部屋で直接燃やす（どれだけ上手に燃やしても 1 まで）より暖かい"))
    # 🫧 2 段目: 冷媒（R32）の輪。冷媒は回り、熱は冷媒に乗る。外気が下がると蒸発温度が下がって COP も下がる
    q = p["r32"]
    loop = steam_loop_plan("r32")
    t_loop = 64.0 * max(w["tau"] for w in loop["water"].values())
    net = steam_loop_network("r32", t_fuel=t_loop)
    checks += _loop_run_checks("r32:loop:", loop, net, net.run(t_loop))
    cat, h3 = p["heat"] / p["power"], q["h"]["condenser"]

    def cop_ratio(eta: float | None = None) -> tuple[float, float]:  # 外気 2 °C の COP と、7 °C に対する比（除霜込み）
        low = r32_h(low=True, eta=eta)
        c = (low["compressor"] - h3) / (low["compressor"] - low["evaporator"])
        return c, c / _refrigerant_plan(eta)["cop"] * q["defrost"]

    cop_low, ratio = cop_ratio()
    (q7, w7), (q2, w2) = q["jsrae"]
    want = (q2 / w2) / (q7 / w7)
    lo_e, use_e, hi_e = q["eta"]
    swept = {e: cop_ratio(e)[1] for e in sweep_points(seed + 4, q["eta"])}
    m = loop["water"]["condenser"]["flow"]
    checks += [
        check("r32:loop:cop", "📐", abs(loop["cop"] / cat - 1.0) <= 0.01, value=loop["cop"], expected=cat, tol=0.01,
              note=f"凝縮 {q['t_cnd']:g} °C・蒸発 {p['t_out'] - (0.100 * q['t_cnd'] + 2.95):.2f} °C で回す冷媒の COP {loop['cop']:.3f}"
                   f" = カタログ 3.6 kW / 705 W（{cat:.3f}）。凝縮温度はこれに合うように選んだ（送風機の分は入れない）"),
        check("r32:loop:outlet", "🌡️", q["t_cnd"] > q["outlet"], value=q["t_cnd"], expected=q["outlet"],
              note=f"冷媒（凝縮 {q['t_cnd']:g} °C）は部屋へ吹き出す空気（測って {q['outlet']:g} °C）より熱い。熱は熱い側からしか流れない"),
        check("r32:cop:outdoor", "📚", abs(ratio / want - 1.0) <= 0.03, value=ratio, expected=want, tol=0.03,
              note=f"外気 7 → 2 °C で冷媒の COP は {loop['cop']:.2f} → {cop_low:.2f}（{cop_low / loop['cop']:.3f} 倍、蒸発温度が下がる）。"
                   f"除霜の補正 {q['defrost']:g} を掛けて {ratio:.3f} 倍。2.8 kW 級（M社）のカタログ 3.6 kW / 910 W → "
                   f"3.5 kW / 1,320 W の {want:.3f} 倍と 3 % 以内"),
        check("r32:cop:outdoor:eta", "🎚️", all(abs(r / want - 1.0) <= 0.03 for r in swept.values()),
              value=max(abs(r / want - 1.0) for r in swept.values()), expected=0.0, tol=0.03,
              note=f"圧縮機の効率を {lo_e:g}〜{hi_e:g} で振っても、7 → 2 °C の COP の比は "
                   + "・".join(f"{e:.2f} で {r:.3f}" for e, r in swept.items())
                   + f"（カタログ {want:.3f}）でどれも 3 % 以内。効率は COP そのものには効くが、比ではほぼ打ち消す"
                     f"（凝縮温度は {use_e:g} で合わせた {q['t_cnd']:g} °C のまま）"),
    ]
    # 🏠 3 段目: 外気で落ちる最大出力 × 部屋から逃げる熱 → バランス点。逃げた熱は外気へ戻る（熱も外気から来て外気へ帰る）
    mx, h = p["max"], p["house"]
    for t, cap in mx["catalog"]:
        got = heat_pump_max(t)
        checks.append(check(f"max:{t:+g}c", "📚", abs(got / cap - 1.0) <= 0.10, value=got, expected=cap, tol=0.10,
                            note=f"外気 {t:g} °C の最大暖房出力: 建築研究所の式 × 除霜 {p['r32']['defrost']:g} = {got:.2f} kW、"
                                 f"カタログ {cap:g} kW（{got / cap - 1:+.1%}）。多くの機種から作った式なので 1 割以内を目安に"
                                 + ("。−15 °C は着霜の条件が書かれていないので除霜ありとして" if t < -10 else "")))
    t_b = house_balance()
    (t1, c1), (t2, c2) = mx["catalog"][1:]
    ua = house_load(h["t100"]) / (h["t0"] - h["t100"])
    t_cat = bisection(lambda t: c1 + (c2 - c1) * (t - t1) / (t2 - t1) < house_load(t), -40.0, t1)  # カタログの 2 点を直線で
    t_iso = house_balance(iso=True)
    checks += [
        check("house:balance", "🏠", abs(t_b - t_cat) <= 1.0, value=t_b, expected=t_cat, tol=1.0,
              note=f"JIS の家（0 °C で {house_load(0.0):.2f} kW、UA {ua * 1e3:.0f} W/K）を 1 台で 20 °C に保てるのは外気 "
                   f"{t_b:.1f} °C まで（建築研究所の式）。カタログの −7・−15 °C を直線で結ぶと {t_cat:.1f} °C で 1 K 以内"),
        check("house:balance:iso", "📚", abs(t_iso - t_b) <= 1.0, value=t_iso, expected=t_b, tol=1.0,
              note=f"ISO 16358-2 のまま（0 °C で 0.82 × 定格暖房 = {house_load(0.0, iso=True):.2f} kW）でもバランス点は {t_iso:.1f} °C。"
                   "JIS の 冷房 × 1.25 との違いは 1 K 以内で、どちらの読みでも結論は同じ"),
    ]
    t_room = 64.0 * h["tau"]
    net = house_network(t_end=t_room)
    res = net.run(t_room)
    for t in h["outs"]:
        q_in = min(house_load(t), heat_pump_max(t)) + ua * (p["t_room"] - h["t0"])
        got = t + from_q(res["volumes_q"][f"room{t:+g}"]) / (ua * h["tau"])
        want = t + q_in / ua  # 一定の入り + 量に比例する出口は、エンジンでも連続の式どおり
        short = t < t_b
        checks.append(close(f"house:room{t:+g}", "🥶" if short else "🏠", got, want, 1e-9,
                            note=f"外気 {t:g} °C の部屋は {t + q_in / ua:.2f} °C で落ち着く（入る熱 {q_in:.2f} kW = 逃げる熱 UA·(室温 − 外気)）"
                                 + ("。バランス点より冷たいので最大出力でも 20 °C に届かない" if short else "。要る分だけ出して 20 °C")))
    checks.append(conserved("house:conserve", res))
    for item in random_scales(seed + 3, 2, n_min=0, n_max=1):  # 🟢 部屋の τ を c 倍（0.1〜9.99）: 落ち着く温度は同じ
        c = item["value"]
        r = house_network(slow=c, t_end=c * t_room, dt=60.0).run(c * t_room)
        worst = max(abs(t + from_q(r["volumes_q"][f"room{t:+g}"]) / (ua * c * h["tau"])
                        - t - (min(house_load(t), heat_pump_max(t)) + ua * (p["t_room"] - h["t0"])) / ua) for t in h["outs"])
        checks.append(check(f"house:slow={item['label']}", "🟢", worst <= 1e-9, value=worst, tol=1e-9,
                            note=f"部屋の τ（仮置きの {h['tau']:g} s）を {c:.3g} 倍にしても、落ち着く室温は同じ（1e-9 K 以内）。"
                                 "τ は温まる速さだけを決め、入る熱 = 逃げる熱 の温度は UA だけで決まる"))
    summary_house = {"ua_kw_k": ua, "balance_c": t_b, "balance_catalog_c": t_cat, "balance_iso_c": t_iso,
                     "room_c": {f"{t:g}": t + (min(house_load(t), heat_pump_max(t)) + ua * (p["t_room"] - h["t0"])) / ua
                                for t in h["outs"]},
                     "max_kw": {f"{t:g}": heat_pump_max(t) for t, _c in mx["catalog"]}}
    summary = {"cop": st["cop"], "carnot_cop": limit, "share_of_carnot": st["cop"] / limit, "lag_s": st["lag"],
               "outdoor_share": 1.0 - 1.0 / st["cop"], "chain": chain, "house": summary_house,
               "r32": {"flow_g_s": m, "loop_s": q["charge"] * 1e3 / m, "cop": loop["cop"], "cop_2c": cop_low,
                       "cop_2c_defrost": cop_low * q["defrost"], "outdoor_kw": loop["fuel"]["outdoor"][1],
                       "compressor_kw": loop["work"]["compressor"], "taus": {k: w["tau"] for k, w in loop["water"].items()}}}
    return {"params": {"heat_pump": p, "seed": seed}, "summary": summary, "checks": checks,
            "runs": {"heat-pump": heat_pump_network().run(p["scene"][0], sample_every=p["scene"][1])["series"],
                     "balance": balance_rows()}}


# --- E33 脳脊髄液（🧠 血液 → 脳室 → くも膜下腔 → 血液） ---------------------------------------

# 時間は分、量は mL。脈絡叢が血液から 0.35 mL/min（約 500 mL/日）作り、脳室 25 mL → くも膜下腔 125 mL を通って
# くも膜顆粒から静脈（血液）へ戻る。作る速さは圧によらない（ポンプ）、吸収は圧による（量に比例する弁）。
# scene は draw の (分, 何ステップごとに描くか): 吸収がふつうの頭 ｜ 止まった頭 を並べる。
CSF = {"production": 0.35, "volumes": {"ventricles": 25.0, "subarachnoid": 125.0}, "blood": 5000.0, "dt": 1.0,
       "scene": (1440.0, 36)}


def csf_taus() -> dict:
    """滞在時間 = 量 / 作る速さ（ふつうの状態で決め、吸収の通りやすさを変えても脳室は同じ）。"""
    return {k: v / CSF["production"] for k, v in CSF["volumes"].items()}


def csf_steady(*, opening: float = 1.0, dt: float | None = None) -> dict:
    """行き渡った後: 脳室 = 作る速さ × τ、くも膜下腔 = 作る速さ × τ / 吸収の開き具合（steady_tank）。"""
    f, tau = CSF["production"], csf_taus()
    vent, base = steady_tank(1.0 / tau["ventricles"], f, f, dt)
    sub, _base = steady_tank(opening / tau["subarachnoid"], f, base, dt)
    return {"ventricles": vent, "subarachnoid": sub, "turnover_per_day": 1440.0 * f / (vent + sub)}


def csf_network(openings=(1.0,), *, steady: bool = False, dt: float | None = None) -> Network:
    """🩸 血液 → 🧠 脳室 → くも膜下腔 → 🩸 血液（頭ごとに -1, -2 …）。opening は吸収（くも膜顆粒）の開き具合、0 なら吸収しない。
    steady=True なら、どの頭もその開き具合の離散の定常量から始める（最初から回っている状態。0 の頭は止まる直前のふつうの量）。"""
    p, tau = CSF, csf_taus()
    net = Network(dt=p["dt"] if dt is None else dt)
    empty = {"ventricles": 0.0, "subarachnoid": 0.0}
    starts = [csf_steady(opening=o if o > 0 else 1.0, dt=net.dt) if steady else empty for o in openings]
    net.tank("blood", 1.0, volume=p["blood"] - sum(st[k] for st in starts for k in p["volumes"]))
    for i, (o, start) in enumerate(zip(openings, starts), 1):
        vent, sub = f"ventricles-{i}", f"subarachnoid-{i}"
        net.tank(vent, 1.0, volume=start["ventricles"]).tank(sub, 1.0, volume=start["subarachnoid"])
        net.pump("blood", vent, rate=p["production"])
        net.valve(vent, sub, tau=tau["ventricles"])
        if o > 0:
            net.valve(sub, "blood", tau=tau["subarachnoid"], opening=o)
    return net


def e33_csf(seed: int = 33) -> dict:
    """🧠 脳脊髄液は血液から作られ、脳室とくも膜下腔を通って血液へ戻る（脳の中の小さな循環）。頭の中では 入る = 出る（±0）、
    量 = 作る速さ × 滞在時間、1 日に 3〜4 回入れ替わる。吸収が止まると、作る速さのまま溜まる。"""
    p, tau, checks = CSF, csf_taus(), []
    dt = p["dt"]
    total = sum(p["volumes"].values())
    t_end = 40.0 * tau["subarachnoid"]
    turn = total / p["production"]
    net = csf_network()
    net.tracer("dye", "ventricles-1", release=t_end - 4.0 * turn, stop="blood")  # 🧪 脳室に入れた色素（先入れ先出し）
    res = net.run(t_end)
    st = csf_steady(dt=dt)
    got = {k: from_q(res["volumes_q"][f"{k}-1"]) for k in p["volumes"]}
    worst = max(rel_err(got[k], st[k]) for k in got)
    rates = node_rates(res)
    made, absorbed = rates["ventricles-1"]["in"], rates["subarachnoid-1"]["out"]
    dye = res["tracers"][0]
    lap = None if dye["arrival"] is None else dye["arrival"] - dye["path"][0]["t"]
    checks += [
        check("steady", "🎯", worst <= 1e-9, value=worst, tol=1e-9,
              note="脳室 = 作る速さ × τ、くも膜下腔も同じ（離散の定常: エンジンの 1 ステップをそのまま解いた式）"),
        close("steady:volume", "📐", sum(got.values()), total, dt / tau["ventricles"],
              note="合計 150 mL（脳室 25 + くも膜下腔 125）。連続の式との差は刻みの分まで"),
        close("skull:±0", "🧠", absorbed, made, 1e-9, note="頭の中は 作る速さ = 吸収する速さ（入る = 出る）"),
        close("turnover", "🔄", 1440.0 * absorbed / sum(got.values()), st["turnover_per_day"], 1e-9,
              note="1 日に入れ替わる回数 = 1 日に作る量 / 量（504 / 150 = 3.36 回）"),
        check("dye", "🧪", lap is not None and abs(lap - turn) <= 2 * dt, value=lap, expected=turn, tol=2 * dt,
              note="脳室に入れた色素が血液に出るまで = 量 / 作る速さ（7.1 時間）"),
        conserved("conserve", res),
    ]
    # 🚫 吸収が止まると、作った分がそのまま溜まる（1 日で 504 mL、整数で 0 ulp）
    day = 1440.0
    blocked = csf_network((0.0,), steady=True).run(day)
    gained = sum(blocked["volumes_q"][f"{k}-1"] for k in p["volumes"]) - sum(to_q(st[k]) for k in p["volumes"])
    pump = next(ln for ln in blocked["links"] if ln["kind"] == "pump")
    checks += [
        check("blocked:gain", "🚫", gained == pump["volume_moved_q"], value=from_q(gained),
              expected=p["production"] * day, tol=0,
              note="吸収が止まると 1 日に作る 504 mL がそのまま溜まる（作るのは圧によらない。Q 整数で 0 ulp）"),
        conserved("blocked:conserve", blocked),
    ]
    for item in random_scales(seed, 3, n_min=1, n_max=1):  # 🎲 吸収の通りやすさ o 倍（0.1〜0.999）
        o = item["value"]
        want = csf_steady(opening=o, dt=dt)
        r = csf_network((o,), steady=True).run(day)  # 定常の点に置いて、1 日たっても動かないか
        sub = from_q(r["volumes_q"]["subarachnoid-1"])
        checks += [
            close(f"open={item['label']}:subarachnoid", "🎲🌊", sub, want["subarachnoid"], 1e-9,
                  note="吸収の通りやすさが o 倍になると、くも膜下腔の量は 1/o 倍で釣り合う（そこから動かない）"),
            close(f"open={item['label']}:closed-form", "📐", want["subarachnoid"], p["volumes"]["subarachnoid"] / o,
                  dt / tau["ventricles"], note="連続の式 125 mL / o との差は刻みの分まで"),
            conserved(f"open={item['label']}:conserve", r),
        ]
    summary = {"steady_ml": got, "turnover_per_day": st["turnover_per_day"], "residence_h": turn / 60.0,
               "dye_lap_min": lap, "blocked_gain_ml_per_day": from_q(gained)}
    return {"params": {"csf": p, "seed": seed}, "summary": summary, "checks": checks,
            "runs": {"csf": res["series"][-50:]}}


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
    "E17": ("viscosity", "🍯", "粘度のつまみ（直接・倍率、仮置きの流体、粘い小孔）", e17_viscosity),
    "E18": ("rain-runoff", "🌧️", "斜面に降る雨（運動波の閉じた式との精度比べ・区間数での収束）", e18_rain_runoff),
    "E19": ("shishi-odoshi", "🎍", "鹿威し（溜まる → 傾く → 一気に落ちる、周期 V_tip/q + t_dump）", e19_shishi_odoshi),
    "E20": ("pythagorean-cup", "🍷", "ピタゴラスの杯（頂まで溜まるとサイフォンが吸い切る）", e20_pythagorean_cup),
    "E21": ("mariotte", "🫗", "マリオットびん（空気管の下端まで流量が一定）", e21_mariotte),
    "E22": ("water-clock", "🏺", "水時計（断面が高さで変わる器: 4 次の放物面なら水位が一定の速さで下がる）", e22_water_clock),
    "E23": ("inertia-limits", "🚧", "慣性が要る題材（U 字管の揺れ・ダム決壊・水撃）は再現しない — 外れ方を数で記す",
            e23_inertia_limits),
    "E24": ("distillation", "⚗️", "ワインの蒸留（Rayleigh 式・共沸点・コニャック式の 2 回蒸留）", e24_distillation),
    "E25": ("angels-share", "👼", "天使の取り分（湿度で度数の上下が入れ替わる・40 % を切る年数）", e25_angels_share),
    "E26": ("money-flow", "💴", "お金の流れ（中央銀行 → メガ銀行 → 証券 → 市中: 量 = 出す速さ × 滞在時間・T+2・決済の手元）",
            e26_money_flow),
    "E27": ("gold-flow", "🪙", "金に置き換える（総量が決まっている: 取り分は滞在時間の比、速く回しても増えない）", e27_gold_flow),
    "E28": ("rare-metal", "⛏️", "レアメタル（産出国 → 船 → 加工国 → 市場: Bateman の鎖、在庫の山の時刻と高さ）", e28_rare_metal),
    "E29": ("outlet-count", "🚿", "出口の数（別々に取ると n 本で刻みの 1 次ずれる: (1 − 1/n) に比例。槽ごとに厳密に解く既定ではずれない）",
            e29_outlet_count),
    "E30": ("circulation", "🫀", "血液の循環（心臓 → 器官 → 心臓: 器官は入る = 出る、V = F·T、一周 60 s、心臓は戻ってきた分しか送れない）",
            e30_circulation),
    "E31": ("thermal-power", "🔥", "火力発電（燃料 → 熱が溜まる段 → 電気・煙突・海: 効率は割合の積、遅れは τ の和、Carnot・Curzon–Ahlborn と比べる。💧 水・蒸気の循環、📉 部分負荷、🌅 起動）",
            e31_thermal_power),
    "E32": ("heat-pump", "❄️", "ヒートポンプ（🔥 の裏返し: 電気 + 外気の熱 → 部屋、COP = 1/η_Carnot が上限、発電所の電気で回すと燃料 1 から部屋に η·COP。🫧 冷媒の輪、外気で変わる COP、🏠 バランス点）",
            e32_heat_pump),
    "E33": ("csf", "🧠", "脳脊髄液（血液 → 脳室 → くも膜下腔 → 血液: 頭の中は入る = 出る、1 日 3〜4 回入れ替わる、吸収が止まると作る速さで溜まる）",
            e33_csf),
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
               "🧱": "noodle-clog", "🍽": "chopsticks", "🌧": "rain-runoff", "🎍": "shishi-odoshi",
               "🍷": "pythagorean-cup", "🫗": "mariotte", "🏺": "water-clock", "🚧": "inertia-limits",
               "⚗": "distillation", "👼": "angels-share", "💴": "money-flow", "🪙": "gold-flow", "⛏": "rare-metal",
               "🚿": "outlet-count", "🫀": "circulation",
               "🪨": "coal-power", "🌀": "combined-cycle", "❄": "heat-pump", "🧠": "csf", "💧": "steam-loop", "♨": "hrsg-loop", "📉": "part-load", "🌅": "startup", "🌡": "metal", "🫧": "r32-loop", "🏠": "house", "🧬": "rna-turnover", "🧫": "rna-competition"}


def fluid_arg(text: str) -> str:
    """--fluid: キーでも 💧 でも 🟦 でもよい。"""
    return fluid_aliases().get(_plain(text), text)


def scene_arg(text: str) -> str:
    """draw の場面: 名前でも 🏛️ 🥢 🎋 🫙 🌊 🕳️ 🧊 🪵 🍜 🧱 🍽️ 🌧️ 🎍 🍷 🫗 🏺 🚧 ⚗️ 👼 💴 🪙 ⛏️ 🚿 🫀 🪨 🌀 ❄️ 🧠 💧 ♨️ 📉 🌅 🌡️ 🫧 🏠 でもよい。"""
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
            f"  ν={f['nu'] * 1e6:.4g} mm²/s" + (f"  🏷️ 仮置き（粘度は水の {PLACEHOLDER_MU_SCALE[f['key']]:g} 倍）"
                                                if f["placeholder"] else "")
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


# --- 📒 記録（静止の図の元になる行）----------------------------------------------------------------
# FFSS が描くのは時間で動く GIF（断面図・🕰️ 時空図・ドット絵）まで。サンキー図・T–Q 線図・トルネード図のような静止の図は
# DRAgoWing が描く: ここでは実験の runs に 1 行 = 1 つの値として入れ、`verify --jsonl` の記録として渡すだけ。

def energy_links(kind: str = "coal") -> tuple[tuple, tuple]:
    """🔀 燃料 100 の行き先（1 段目の割合のまま。エンジンは thermal_network でこれを ±0 で運ぶ）。
    (ノード (名前, 組), リンク (元, 先, 量))。組は heat・water・power・loss。"""
    p = thermal_plan(kind)
    st = {name: outs for name, _t, outs in p["stages"]}
    gross = 100.0 * p["efficiency"]
    aux = gross * p["aux"]
    if kind == "coal":
        steam = 100.0 * st["boiler"]["steam"]
        nodes = (("燃料", "heat"), ("ボイラ", "heat"), ("煙突", "loss"), ("蒸気", "water"), ("発電端", "power"),
                 ("海（復水器）", "loss"), ("送電端", "power"), ("所内", "loss"))
        return nodes, ((0, 1, 100.0), (1, 2, 100.0 - steam), (1, 3, steam), (3, 4, gross), (3, 5, steam - gross),
                       (4, 6, gross - aux), (4, 7, aux))
    exhaust = 100.0 * st["gt"]["hrsg"]
    steam = exhaust * st["hrsg"]["steam"]
    st_power = steam * st["steam"]["power"]
    nodes = (("燃料", "heat"), ("ガスタービン", "heat"), ("排熱回収ボイラ", "heat"), ("煙突", "loss"), ("蒸気", "water"),
             ("海（復水器）", "loss"), ("発電端", "power"), ("送電端", "power"), ("所内", "loss"))
    return nodes, ((0, 1, 100.0), (1, 6, 100.0 - exhaust), (1, 2, exhaust), (2, 3, exhaust - steam), (2, 4, steam),
                   (4, 6, st_power), (4, 5, steam - st_power), (6, 7, gross - aux), (6, 8, aux))


def blood_links() -> tuple[tuple, tuple]:
    """🔀 🫀 心拍出量 [L/min] の行き先（教科書の器官の割合）→ 静脈 → 肺。形は energy_links と同じ。"""
    p = CIRCULATION
    names = {"brain": "脳", "coronary": "冠動脈", "kidney": "腎臓", "gut-liver": "消化管・肝臓", "muscle": "筋肉",
             "skin": "皮膚", "other": "その他"}
    out, organs = p["output"] * 60.0, list(p["flows"])
    nodes = (("心臓（左）", "heart"), *((names[o], "organ") for o in organs), ("静脈", "vein"), ("肺", "lung"))
    vein = len(organs) + 1
    return nodes, (*((0, i + 1, out * p["flows"][o]) for i, o in enumerate(organs)),
                   *((i + 1, vein, out * p["flows"][o]) for i, o in enumerate(organs)), (vein, vein + 1, out))


def flow_rows(nodes, links, unit: str) -> list[dict]:
    """リンク 1 本 = 1 行（元・先の名前と組、量、単位）。"""
    return [{"source": nodes[a][0], "target": nodes[b][0], "source_group": nodes[a][1], "target_group": nodes[b][1],
             "value": v, "unit": unit} for a, b, v in links]


def tq_lines(pinch: float | None = None) -> dict:
    """🌡️ 🌀 排熱回収ボイラの T–Q（煙突の側から渡した熱 [MW] と温度 [°C]）。ガスは一直線（比熱一定）、水・蒸気は 3 圧の
    合成: 蒸発は平らで量は h_g − h_f から正確、節炭・過熱はつなぐ直線（比熱一定の近似）。復水の温度は h_f ÷ 4.18。"""
    q = STEAM_LOOP["lng"]
    pl = steam_loop_plan("lng", pinch=pinch)
    w = pl["water"]
    m_hp, m_ip, m_lp = (w[n]["flow"] for n in ("hp-hrsg", "ip-hrsg", "lp-hrsg"))
    (t_h, f_h, g_h), (t_i, f_i, g_i), (t_l, f_l, g_l) = q["drums"]
    h_c, stack, q_all = q["saturation"][0], q["plant"]["stack"], pl["q_in"]
    c_g = q_all / (pl["t_exh"] - stack)
    q_lp = (m_hp + m_ip + m_lp) * (f_l - h_c) / 1e3
    q_ip, q_hp = c_g * (t_i + pl["pinch"] - stack), c_g * (t_h + pl["pinch"] - stack)
    water = [(0.0, h_c / 4.18), (q_lp, t_l), (q_lp + m_lp * (g_l - f_l) / 1e3, t_l), (q_ip, t_i),
             (q_ip + m_ip * (g_i - f_i) / 1e3, t_i), (q_hp, t_h), (q_hp + m_hp * (g_h - f_h) / 1e3, t_h),
             (q_all, q["t_steam"])]
    pinches = [(x, t, stack + x / c_g) for x, t in ((q_lp, t_l), (q_ip, t_i), (q_hp, t_h))]
    return {"gas": [(0.0, stack), (q_all, pl["t_exh"])], "water": water, "pinches": pinches,
            "lp_join": q_lp + m_lp * (g_l - f_l) / 1e3 + (m_hp + m_ip) * (f_i - f_l) / 1e3}


def tq_rows(pinch: float | None = None) -> list[dict]:
    """T–Q の点 1 つ = 1 行（stream は gas・water・pinch。pinch の行はガスと水の温度と差 ΔT、どの圧か）。"""
    tq = tq_lines(pinch)
    rows = [{"stream": k, "q_mw": x, "t_c": t} for k in ("gas", "water") for x, t in tq[k]]
    return rows + [{"stream": "pinch", "drum": d, "q_mw": x, "t_c": t_w, "t_gas_c": t_g, "dt_k": t_g - t_w}
                   for (x, t_w, t_g), d in zip(tq["pinches"], ("lp", "ip", "hp"))]


def startup_rows(at: float = 0.5) -> list[dict]:
    """📅 🌅 起動の手順: 金属によらない手順 → 🧱 暖機 → 最低負荷から全負荷へ [分]（hot = 8 h・cold = 48 h 止まった後）。
    agora は出典の起動時間（並列まで）の範囲。"""
    rows = []
    for kind, p in THERMAL["plants"].items():
        for start in ("hot", "cold"):
            plan = thermal_startup(kind, start, at=at)
            warm = metal_warmup(kind, THERMAL["metal"]["off"][start]) / 60.0
            sync, full = plan[1][0] / 60.0, plan[-1][0] / 60.0
            lo, hi = p["startup"][start]
            for phase, a, b in (("steps", 0.0, sync - warm), ("warmup", sync - warm, sync), ("ramp", sync, full)):
                rows.append({"plant": kind, "start": start, "phase": phase, "from_min": a, "to_min": b,
                             "agora_lo_min": lo, "agora_hi_min": hi, "at": at})
    return rows


def sensitivity_rows(samples: int = 33) -> list[dict]:
    """🎚️ 振った仮定ごとに、🔮 のずれ ÷ 許容（1 以下なら 🔮 のまま）を幅の中の samples 点で。use は使った値の行。"""
    q, hp, m = STEAM_LOOP["lng"], HEAT_PUMP, THERMAL["metal"]
    nn = [f / sum(q["plant"]["steam"]) for f in q["plant"]["steam"]]

    def pinch(x: float) -> float:
        fl = [steam_loop_plan("lng", pinch=x)["water"][n]["flow"] for n in ("hp-hrsg", "ip-hrsg", "lp-hrsg")]
        return max(abs(f / sum(fl) - b) for f, b in zip(fl, nn)) / 0.03

    def metal(tau: float, rate: float) -> float:
        out = 0.0
        for kind, p in THERMAL["plants"].items():
            lo, hi = p["startup"]["cold"]
            for at in (0.0, 1.0):
                c = thermal_startup(kind, "cold", at=at, tau=tau, rate=rate)[1][0] / 60.0
                out = max(out, abs(c - (lo + hi) / 2) / ((hi - lo) / 2))
        return out

    (q7, w7), (q2, w2) = hp["r32"]["jsrae"]
    want, h3 = (q2 / w2) / (q7 / w7), hp["r32"]["h"]["condenser"]

    def eta(e: float) -> float:
        low = r32_h(low=True, eta=e)
        c = (low["compressor"] - h3) / (low["compressor"] - low["evaporator"])
        return abs(c / _refrigerant_plan(e)["cop"] * hp["r32"]["defrost"] / want - 1.0) / 0.03

    share = STEAM_LOOP["power"] / steam_loop_plan("coal")["q_in"]
    eff = thermal_plan("coal")["efficiency"]
    spans = (("pinch", "K", pinch, q["pinch"]), ("metal-rate", "K/min", lambda r: metal(m["tau"][1], r), m["rate"]),
             ("metal-tau", "h", lambda t: metal(t, m["rate"][1]), m["tau"]), ("compressor-eta", "", eta, hp["r32"]["eta"]),
             ("boiler", "", lambda b: abs(share - eff / b) / 0.03, (0.90, 0.94, 0.96)))
    rows = []
    for name, unit, f, (lo, use, hi) in spans:  # ボイラの幅だけは出典がなく仮（source_range = False）
        for x, used in [(lo + (hi - lo) * k / (samples - 1), False) for k in range(samples)] + [(use, True)]:
            rows.append({"assumption": name, "x": x, "unit": unit, "metric": f(x), "use": used, "lo": lo, "hi": hi,
                         "source_range": name != "boiler"})
    return rows


def rate_window_rows() -> list[dict]:
    """🧱 τ × 温め直しの速さ: cold の起動時間の Agora の範囲からのずれ ÷ 半幅（1 以下が中）。"""
    lo_hi = {k: p["startup"]["cold"] for k, p in THERMAL["plants"].items()}

    def worst(tau: float, rate: float) -> float:
        return max(abs(thermal_startup(k, "cold", at=at, tau=tau, rate=rate)[1][0] / 60.0 - (lo + hi) / 2) / ((hi - lo) / 2)
                   for k, (lo, hi) in lo_hi.items() for at in (0.0, 1.0))

    return [{"tau_h": 60.0 + 2.0 * i, "rate_k_min": 1.0 + 0.1 * j, "metric": worst(60.0 + 2.0 * i, 1.0 + 0.1 * j)}
            for i in range(11) for j in range(47)]


def balance_rows() -> list[dict]:
    """✖️ 🏠 外気ごとの最大暖房出力（建築研究所の式 × 除霜）と家の暖房負荷（JIS・ISO）[kW]。catalog はカタログの低温の点。"""
    rows = [{"kind": "curve", "t_out": t, "max_kw": heat_pump_max(t), "load_jis_kw": house_load(t),
             "load_iso_kw": house_load(t, iso=True)} for t in (-20.0 + 0.5 * k for k in range(75))]
    return rows + [{"kind": "catalog", "t_out": t, "max_kw": q} for t, q in HEAT_PUMP["max"]["catalog"]] + \
        [{"kind": "balance", "t_out": house_balance(), "iso": False}, {"kind": "balance", "t_out": house_balance(iso=True),
                                                                      "iso": True}]


def rna_turnover_scene():
    """同じ量の RNA から、合成速度と半減期を個別に変える。単位は任意量と時間。"""
    net = Network(dt=1 / 32)
    initial = 10.0
    base_rate = initial * math.log(2) / 2
    # 有限の材料槽から作り、分解した量は回収槽へ送る。
    net.tank("supply", 1.0, volume=100).tank("degraded", 1.0)
    columns = []
    for name, synthesis, half_life in (("baseline", base_rate, 2),
                                       ("double-synthesis", 2 * base_rate, 2),
                                       ("double-half-life", base_rate, 4)):
        columns.append(name)
        net.tank(name, 1.0, volume=initial)
        net.pump("supply", name, rate=synthesis, until=4)
        net.valve(name, "degraded", tau=half_life / math.log(2))
    result = net.run(8, sample_every=4)
    result["series"].insert(0, {"t": 0.0, **{f"h:{name}": initial for name in columns}})
    result["caps"] = dict.fromkeys(columns, 20.0)
    return result, columns


@register_law("mirna-turnover")
def _law_mirna_turnover(net, link, step):
    """活性のある同一miRNA経路を共有する標的に、有限の分解流量を配分する。"""
    targets = sum(net.nodes[name]["v"] for name in link["targets"])
    share = qdiv(net.nodes[link["src"]]["v"], link["half_saturation_q"] + targets)
    return qmul(link["capacity_per_step_q"], share)


def rna_competition_scene(added_percent: float = 10.0):
    """血漿中のAとBをpMで追跡する。投入直後の B/(A+B) を条件とする。"""
    if not 0 <= added_percent < 100:
        raise ValueError("added_percent must be in [0, 100)")
    net = Network(dt=1 / 64, inflows="end")
    initial = 10.0
    columns = []
    for condition in ("unshared", "shared"):
        original, added = f"{condition}-A", f"{condition}-B"
        columns.extend((original, added))
        targets = (original, added) if condition == "shared" else (original,)
        for name in (original, added):
            net.tank(name, 1.0, volume=initial if name == original else 0.0)
            net.tank(f"{name}-basal", 1.0).tank(f"{name}-mirna", 1.0)
            net.valve(name, f"{name}-basal", tau=4 / math.log(2))
            if name in targets:
                net.link("mirna-turnover", name, f"{name}-mirna", targets=targets,
                         half_saturation_q=to_q(5), capacity_per_step_q=to_q(2 * net.dt))
    # 投入前の軌跡を同じ法則で求め、その時点のAを基準に投入量を逆算する。
    before = deepcopy(net)
    before.run(2)
    pulse = before.head("unshared-A") * added_percent / (100 - added_percent)
    for condition in ("unshared", "shared"):
        supply = f"{condition}-supply"
        net.tank(supply, 1.0, volume=pulse)
        net.pump(supply, f"{condition}-B", rate=pulse / net.dt, start=2 - net.dt, until=2)
    start = {"t": 0.0, **{f"h:{name}": net.head(name) for name in net.nodes}}
    result = net.run(8, sample_every=8)
    result["series"].insert(0, start)
    result["pulse_pm"] = pulse
    # GIFの満水位はAが初期濃度、Bが投入濃度。4列とも各RNAの残存率を示す。
    result["caps"] = {name: initial if name.endswith("-A") else (pulse or initial) for name in columns}
    return result, columns


def plasma_rna_report(percentages, blood_l: float, plasma_percent: float):
    """1 pM = 1 fmol/mL。仮定した血漿中濃度から血漿全体の量(pmol)を換算する。"""
    plasma_l = blood_l * plasma_percent / 100
    rows = []
    for percent in percentages:
        result, _ = rna_competition_scene(percent)
        for row in result["series"]:
            for condition in ("unshared", "shared"):
                a, b = (row[f"h:{condition}-{rna}"] for rna in ("A", "B"))
                for rna, concentration in (("A", a), ("B", b)):
                    name = f"{condition}-{rna}"
                    supplied = result["caps"][name] if rna == "A" else (result["pulse_pm"] if row["t"] >= 2 else 0)
                    rows.append({"added_percent": percent, "hour": row["t"], "condition": condition, "rna": rna,
                                 "fmol_per_ml": concentration, "total_pmol": concentration * plasma_l,
                                 "percent_of_pair": 100 * concentration / (a + b) if a + b else None,
                                 "percent_remaining": 100 * concentration / supplied if supplied else None,
                                 "basal_loss_fmol_per_ml": row[f"h:{name}-basal"],
                                 "mirna_loss_fmol_per_ml": row[f"h:{name}-mirna"]})
    return {"blood_l": blood_l, "plasma_percent": plasma_percent, "plasma_l": plasma_l, "rows": rows}


def rig_section(scene: str, fluid: str = "water", temp_c: float | None = None, mu: float | None = None,
                mu_scale: float | None = None) -> tuple[dict, list[str]]:
    """断面図で見る試験台（run の結果, 列）。moniac-rig は頭・貯める・戻りの管・横の槽（下の槽は水源なので省く）、
    somen-flume は竹 3 本と受け皿・床、cold-water は 10 °C の細い流れ、bamboo-wear は節と割れ目のある竹。
    temp_c は水温（somen-flume・cold-water・bamboo-wear。cold-water は省くと 10 °C）、mu・mu_scale は粘度（🍯）。
    moniac-rig の弁は時定数で決めているので、粘度を変えても流れは変わらない。"""
    visc = {"mu": mu, "mu_scale": mu_scale}
    if scene == "moniac-rig":
        return moniac_rig(fluid=fluid).run(1200.0), ["header", "store", "store>sump", "side"]
    if scene == "somen-flume":
        net = somen_layout(5e-5, fluid=fluid, temp_c=temp_c, **visc)
        for i in range(3):
            net.tracer(f"somen-{i + 1}", "take/0:0", release=20.0 + 2 * i, stop="basin", size=NOODLES["somen"])
        return net.run(40.0, sample_every=20), ["take/0", "take/1", "take/2", "basin", "floor"]
    if scene == "noodle-length":
        net = somen_layout(5e-5, fluid=fluid, temp_c=temp_c, **visc)
        net.tracer("long", "take/0:0", release=20.0, stop="basin", length=0.3, **STRAND)
        return _from(net.run(32.0, sample_every=10), 20.0), ["take/0", "take/1", "take/2"]
    if scene == "noodle-clog":
        net = somen_layout(3e-6, law="laminar", beta=3.0, dt=0.05, fluid=fluid, temp_c=temp_c, **visc)
        for i in range(16):
            net.tracer(f"bunch-{i + 1:02d}", "take/0:0", release=80.0 + 5.0 * i, stop="basin", length=0.2,
                       strands=10, **STRAND)
        return _from(net.run(200.0, sample_every=20), 100.0), ["take/1", "take/2", "basin"]
    if scene == "chopsticks":
        return _from(catch_run(5e-5, noodles=10), 24.0), ["take/1", "take/2", "basin"]
    if scene == "cold-water":
        net = somen_layout(3e-6, law="laminar", beta=3.0, dt=0.05, fluid=fluid,
                           temp_c=10.0 if temp_c is None else temp_c, **visc)
        net.tracer("hiyamugi", "take/0:0", release=40.0, stop="basin", size=NOODLES["hiyamugi"])
        net.tracer("somen", "take/0:0", release=45.0, stop="basin", size=NOODLES["somen"])
        return net.run(110.0, sample_every=10), ["take/0", "take/1", "take/2", "basin"]
    if scene == "bamboo-wear":
        net = bamboo_layout(3e-6, law="laminar", beta=3.0, dt=0.05, fluid=fluid, temp_c=temp_c, **visc)
        for i, name in enumerate(("udon", "somen")):
            net.tracer(name, "take/0:0", release=40.0 + i, stop="basin", size=NOODLES[name])
        return net.run(100.0, sample_every=10), ["take/0", "take/0b", "basin", "floor"]
    if scene == "rain-runoff":  # 🌧️ 雨が板に行き渡って（上流ほど浅い）、やむと下流から引いていく
        _net, res, _eq = rain_plane(RAIN["manning"], 20, fluid=fluid, sample_every=10, **visc)
        return res, ["plane", "basin"]
    if scene == "shishi-odoshi":  # 🎍 筒がのこぎり歯、受け皿が階段
        net = shishi_layout(2e-6, dt=0.05, t_dump=1.5, dump_rate=1e-4, fluid=fluid)
        return net.run(3.2 * tip_period(TUBE["v_tip"], 2e-6, 1.5), sample_every=4), ["tube", "basin"]
    if scene == "pythagorean-cup":  # 🍷 杯が頂まで溜まって吸い切られ、鉢が階段
        net = cup_layout(5e-6, dt=0.02, fluid=fluid)
        return net.run(3.0 * 71.0, sample_every=10), ["cup", "bowl"]
    if scene == "mariotte":  # 🫗 びんは直線で下がり、空気管の下端から先は Torricelli で鈍る
        net = bottle_layout(BOTTLE["h_air"], dt=0.02, fluid=fluid)
        return net.run(250.0, sample_every=25), ["bottle", "basin"]
    if scene == "water-clock":
        h0, r_top, a = CLOCK["h0"], CLOCK["r_top"], math.pi * (CLOCK["orifice_d"] / 2) ** 2
        k = clock_rate(h0, r_top, CLOCK["orifice_d"], CLOCK["cd"])
        net = Network((h0 / k) / 400.0, fluid, **visc).tank("sea", 1.0)
        net.tank("clepsydra", h0, head=h0, profile=clepsydra_area(h0, r_top))
        net.tank("cylinder", math.pi * r_top * r_top, head=h0)
        net.tank("cone", h0, head=h0, profile=cone_area(h0, r_top, 0.5 * r_top))
        for name in ("clepsydra", "cylinder", "cone"):
            net.orifice(name, "sea", area=a, cd=CLOCK["cd"])
        return net.run(h0 / k, sample_every=1), ["clepsydra", "cylinder", "cone"]
    if scene == "distillation":  # ⚗️ 釜が減り、前留・ハート・セコンド・後留の受け器が順に溜まる
        ne, nw = moles_for_abv(30.0, 0.010)
        net = still_network(ne, nw, dt=2.0, receivers=("heads", "heart", "secondes", "tails"),
                            cuts=[{"to": "heads", "volume": 1e-4}, {"to": "heart", "abv_below": STILL["heart_cut"]},
                                  {"to": "secondes", "abv_below": STILL["secondes_cut"]}, {"to": "tails"}])
        return net.run(3000.0, sample_every=15), ["pot", "heads", "heart", "secondes", "tails"]
    if scene == "inertia-limits":
        # ダム決壊（勾配 1 %）。水位は奥の壁から下がり始め、先端は尖らずにだらだら伸びる（Ritter とは別物）
        net, _info = dam_break(0.1, slope=1e-2, dt=0.02, fluid=fluid)
        return net.run(3.0, sample_every=3), ["bed"]
    if scene == "angels-share":  # 👼 4 つの蔵の樽が 30 年かけて減っていく（台湾は速く、スコットランドはゆっくり）
        rates = angel_rates()
        climates = {k: {**c, "f": (1.0 if k == "scotland" else fit_temperature(rates, c["rh"], c["loss"]))}
                    for k, c in CLIMATES.items()}
        res = angel_cask(rates, climates, steps_per_year=12).run(30 * YEAR, sample_every=12)
        for row in res["series"]:
            row["t"] = round(row["t"] / YEAR, 6)  # ⏱️ は年で出す
        return res, list(CLIMATES)
    if scene in HUB_SCENES:  # 💴 🪙 ⛏️ 同じ配管に流す物を置き換えた 3 つ（⏱️ は日）
        kind, days, every = HUB_SCENES[scene]
        net = hub_network(kind)
        res = net.run(days, sample_every=every)
        res["caps"] = hub_caps(kind)  # 縦いっぱいを槽どうしで揃える（どこにどれだけあるかを比べられる）
        return res, hub_columns(kind, net)
    if scene == "outlet-count":  # 🚿 別々に取る 1・2・4・8 本 ｜ 厳密に解く 8 本（刻みを粗くしてずれを見せる）
        o = OUTLET_SCENE
        net = outlet_network([[1.0 / n] * n for n in o["counts"]], tau=o["tau"], dt=o["dt"], modes=o["modes"])
        res = net.run(o["t_end"], sample_every=1)
        cols = [f"pool-{i}" for i in range(1, len(o["counts"]) + 1)]
        res["caps"] = dict.fromkeys(cols, OUTLETS["volume"])
        return res, cols
    if scene == "circulation":  # 🫀 最初は全部心臓に。送り出して器官・静脈・肺を通って戻る（⏱️ は秒、縦いっぱいは槽ごと）
        t_end, every = CIRCULATION["scene"]  # 🛑 stop で拍動を止め、全部心臓に戻るまで（往路 ＋ 復路で一巡）
        res = circulation_network(stop=CIRCULATION["stop"]).run(t_end, sample_every=every)
        res["gif_frames"] = round(pw.SECTION_GIF_FRAMES * t_end / CIRCULATION["stop"])  # 1 コマの秒は往路だけのときと同じ
        return res, circulation_columns()
    if scene in THERMAL_SCENES:  # 🪨 🌀 燃料を燃やし切るまで。出口（送電・所内・煙突・海）は縦いっぱいを揃える（⏱️ は秒）
        kind = THERMAL_SCENES[scene]
        t_end, every = THERMAL["scene"]
        res = thermal_network(kind).run(t_end, sample_every=every)
        sinks = ("grid", "aux", "stack", "sea")
        res["caps"] = dict.fromkeys(sinks, max(from_q(res["volumes_q"][n]) for n in sinks))
        return res, thermal_columns(kind)
    if scene == "heat-pump":  # ❄️ 電気と外気の熱を使い切るまで。🔌 🌬️ 🏠 は縦いっぱいを揃える（部屋 = 電気 + 外気）
        t_end, every = HEAT_PUMP["scene"]
        res = heat_pump_network().run(t_end, sample_every=every)
        res["caps"] = dict.fromkeys(("power", "outdoor", "room"), from_q(res["volumes_q"]["room"]))
        return res, heat_pump_columns()
    if scene == "rna-turnover":
        return rna_turnover_scene()
    if scene == "rna-competition":
        return rna_competition_scene()
    if scene == "csf":  # 🧠 吸収がふつうの頭 ｜ 止まった頭（どちらも回っている状態から 1 日。⏱️ は分）。くも膜下腔は縦いっぱいを揃える
        t_end, every = CSF["scene"]
        res = csf_network((1.0, 0.0), steady=True).run(t_end, sample_every=every)
        cols = ["ventricles-1", "subarachnoid-1", "ventricles-2", "subarachnoid-2"]
        res["caps"] = dict.fromkeys(cols[1::2], max(from_q(res["volumes_q"][c]) for c in cols[1::2]))
        return res, cols
    if scene == "house":  # 🏠 外気 2 °C ｜ −15 °C の部屋（高さ = 室温 − 外気、縦いっぱい = 20 °C、⏱️ は h）
        t_end, every = HEAT_PUMP["house"]["scene"]
        res = house_network().run(t_end, sample_every=every)
        for row in res["series"]:
            row["t"] = round(row["t"] / 3600.0, 6)
        outs = HEAT_PUMP["house"]["outs"]
        res["caps"] = {f"room{t:+g}": HEAT_PUMP["t_room"] - t for t in outs}
        return res, [f"room{t:+g}" for t in outs]
    if scene == "metal":  # 🌡️ 🪨 止まって 8 h ｜ 48 h で温め直す厚い金属（高さ = 外気より何 K 熱いか、⏱️ は h）
        t_end, every = THERMAL["metal"]["scene"]
        offs = tuple(THERMAL["metal"]["off"].values())
        res = metal_network("coal", offs, t_end=t_end).run(t_end, sample_every=every)
        for row in res["series"]:
            row["t"] = round(row["t"] / 3600.0, 6)
        cols = [f"off-{off:g}h" for off in offs]
        res["caps"] = dict.fromkeys(cols, THERMAL["plants"]["coal"]["startup"]["metal"] - THERMAL["t_cold"])
        return res, cols
    if scene in ("part-load", "startup"):  # 📉 🪨 全負荷で燃やし、途中で 50 % へ（無負荷の燃料・所内電力の固定の分は減らない）
        t_end, every = THERMAL["scene"]  # 🌅 止まった状態から hot 起動: 暖機（無負荷の燃料は海へ）→ 25 % で並列 → 全負荷へ
        load = ((0.0, 1.0), (t_end / 2, 0.5))
        if scene == "startup":
            load = thermal_startup("coal")
            t_end, every = load[-1][0] + 4.0 * max(tau for _n, tau, _o in THERMAL["plants"]["coal"]["stages"]), 2 * every
        res = thermal_network("coal", load=load, willans=True, t_fuel=t_end).run(t_end, sample_every=every)
        sinks = ("grid", "aux", "stack", "sea")
        res["caps"] = dict.fromkeys(sinks, max(from_q(res["volumes_q"][n]) for n in sinks))
        return res, thermal_columns("coal")
    if scene in ("steam-loop", "hrsg-loop", "r32-loop"):  # 💧 ♨️ 冷えた水から火を入れる（🫧 冷媒は凝縮した液から）。
        kind = {"steam-loop": "coal", "hrsg-loop": "lng", "r32-loop": "r32"}[scene]  # 高さは比エンタルピー、縦いっぱいは揃える
        t_end, every = STEAM_LOOP["scene"] if kind != "r32" else HEAT_PUMP["r32"]["scene"]
        res = steam_loop_network(kind).run(t_end, sample_every=every)
        cols = steam_loop_columns(kind)
        res["caps"] = dict.fromkeys(cols, max(w["h"] for w in steam_loop_plan(kind)["water"].values()) / 1e3)
        return res, cols
    raise KeyError(f"unknown rig scene {scene!r}")


def _from(res: dict, t0: float) -> dict:
    """断面図に出す時刻を t0 から後に絞る（水が行き渡って麺を流し始めてから）。"""
    return {**res, "series": [row for row in res["series"] if row["t"] >= t0]}


HUB_SCENES = {"money-flow": ("money", 64.0, 4), "gold-flow": ("gold", 64.0, 4), "rare-metal": ("rare-metal", 320.0, 8)}
THERMAL_SCENES = {"coal-power": "coal", "combined-cycle": "lng"}
RIG_SCENES = ("moniac-rig", "somen-flume", "noodle-length", "noodle-clog", "chopsticks", "cold-water",
              "bamboo-wear", "rain-runoff", "shishi-odoshi", "pythagorean-cup", "mariotte", "water-clock",
              "inertia-limits", "distillation", "angels-share", *HUB_SCENES, "outlet-count", "circulation", *THERMAL_SCENES, "heat-pump", "csf", "steam-loop", "hrsg-loop", "part-load", "startup", "metal", "r32-loop", "house", "rna-turnover", "rna-competition")


def _visc(args) -> dict:
    """--mu（mPa·s）と --viscosity-scale を Network の mu（Pa·s）・mu_scale へ。"""
    return {"mu": None if args.mu is None else args.mu * 1e-3, "mu_scale": args.viscosity_scale}


FILE_ICONS = {"gif": "🎞️", "png": "🖼️", "json": "🧾"}


def section_run(scene: str, fluid: str = "water", temp_c: float | None = None, mu: float | None = None,
                mu_scale: float | None = None) -> tuple[dict, list[str]]:
    """断面図で見る場面（試験台 RIG_SCENES と pixelwater の SECTION_CASES）を回す: (run の結果, 列)。"""
    if scene in RIG_SCENES:
        return rig_section(scene, fluid, temp_c, mu, mu_scale)
    return pw.section_case(scene, fluid, mu=mu, mu_scale=mu_scale)


def _draw(args) -> int:
    files = args.view in ("files", "kymo")
    if args.scene in RIG_SCENES or args.scene in pw.SECTION_CASES:
        res, cols = section_run(args.scene, args.fluid, args.temp_c, **_visc(args))
        if files:
            rows = args.rows or pw.SECTION_GIF_ROWS
            frames = args.frames or res.get("gif_frames", pw.SECTION_GIF_FRAMES)
            replay = {k: v for k, v in {"temp_c": args.temp_c, **_visc(args)}.items() if v is not None}
            out = pw.draw_section(res, cols, args.scene, args.out_dir, frames=frames, rows=rows,
                                  cap_mode=args.cap_mode, fill=args.fill, scale=args.scale,
                                  formats=_formats(args), file_prefix=args.file_prefix,
                                  file_suffix=args.file_suffix, plan=replay, kymo=args.view == "kymo")
            for kind, path in out["paths"].items():
                print(f"{FILE_ICONS[kind]} {path}")
            print(f"🎬 {out['frames']}  📐 {out['canvas']}  📏 ±{out['half_block']:.3g} m  "
                  f"⚖️ {out['conservation_error_q']}")
            return 0 if out["conservation_error_q"] == 0 else 1
        frames = 1 if args.view == "frame" else args.frames or 6
        print(pw.render_frames(pw.network_frames(res, cols, frames=frames, rows=args.rows or pw.SECTION_ROWS,
                                                 cap_mode=args.cap_mode, fill=args.fill)))
        return 0 if res["conservation_error_q"] == 0 else 1
    if args.view == "kymo":
        print("🕰️ 時空図は断面図の場面だけ（ドット絵の場面は frames・frame・files）", file=sys.stderr)
        return 2
    if files:
        out = pw.draw(args.scene, args.fluid, args.out_dir, ticks=args.ticks, seed=args.seed, scale=args.scale,
                      frame_every=args.frame_every, formats=_formats(args), file_prefix=args.file_prefix,
                      file_suffix=args.file_suffix, **_visc(args))
        for kind, path in out["paths"].items():
            print(f"{FILE_ICONS[kind]} {path}")
        print(f"🔢 {out['pixels']} px  🎬 {out['frames']}  📐 {out['canvas']}  ⏱️ {out['settle_tick']}  ⚖️ {out['count_error']}")
        return 0 if out["count_error"] == 0 else 1
    res = pw.run_pixels(args.scene, args.fluid, ticks=args.ticks, seed=args.seed, frame_every=args.frame_every,
                        **_visc(args))
    frames = 1 if args.view == "frame" else args.frames or 6
    print(pw.render_frames(pw.pixel_frames(res, frames=frames, fill=args.fill)))
    return 0 if res["count_error"] == 0 else 1


def _formats(args) -> tuple[str, ...]:
    return ("gif", "png") if args.format == "both" else (args.format,)


# --- 🎞️ README のデモ（場面ごとに GIF 1 枚） ----------------------------------------

# (場面, 流体, 塗り[, 見せ方])。塗りは色だけ（物理は流体で決まる）。⚗️ はエタノールの 🟨、👼 は琥珀の 🟧、🫀 は血の 🟥、
# 🪨 🌀 は熱の 🟧 で塗る。見せ方 "kymo" は 🕰️ 時空図（同じ場面を、横 = 列・縦 = 時刻で。既定は断面図）
DEMOS = (("moniac-rig", "water", None), ("moniac-drain", "water", None), ("flume", "water", None),
         ("flume-hopper", "water", None), ("somen-flume", "water", None), ("noodle-length", "water", None),
         ("noodle-clog", "water", None), ("chopsticks", "water", None), ("cold-water", "water", None),
         ("bamboo-wear", "water", None), ("rain-runoff", "water", None), ("shishi-odoshi", "water", None),
         ("pythagorean-cup", "poison-water", None), ("mariotte", "water", None), ("water-clock", "water", None),
         ("dam-break", "namekian-water", None), ("inertia-limits", "water", None), ("distillation", "water", "🟨"),
         ("angels-share", "water", "🟧"), ("money-flow", "water", "🟩"), ("gold-flow", "water", "🟨"),
         ("rare-metal", "water", "🟫"), ("outlet-count", "water", None), ("circulation", "water", "🟥"),
         ("coal-power", "water", "🟧"), ("combined-cycle", "water", "🟧"),
         ("heat-pump", "water", "🟧"), ("csf", "water", None),
         ("steam-loop", "water", "🟧"), ("hrsg-loop", "water", "🟧"),
         ("part-load", "water", "🟧"), ("startup", "water", "🟧"),
         ("metal", "water", "🟥"), ("r32-loop", "water", "🟧"), ("house", "water", "🟥"),
         ("flume", "water", None, "kymo"), ("rain-runoff", "water", None, "kymo"),
         ("rna-turnover", "water", None), ("rna-competition", "water", None))
DEMOS_SCHEMA = "fluidfullspeedsimulation-demos/1"
_NEEDS_VS16 = set("🏛🕳🍽🌧⚗⛏❄♨🌡")  # 絵文字で出すには U+FE0F が要るもの


def scene_icon(scene: str) -> str:
    icon = next((k for k, v in SCENE_ICONS.items() if v == scene), "🎞")
    return icon + "️" if icon in _NEEDS_VS16 else icon


def demo_gif(scene: str, fluid: str = "water", fill: str | None = None, *, scale: int = 8,
             kymo: bool = False) -> tuple[bytes, dict]:
    """デモの GIF 1 枚（既定の引数だけ）。断面図は 16 段・40 時刻（場面が gif_frames を持てばその数）、
    ドット絵は draw --view files と同じ。kymo=True なら同じ段・時刻の 🕰️ 時空図。"""
    if scene in pw.SCENES:
        res = pw.run_pixels(scene, fluid)
        data = pw.gif_bytes(res, scale=scale)
        return data, {"frames": len(res["frames"]), "canvas": f"{res['width'] * scale}x{res['height'] * scale}",
                      "kind": "pixel", "conserved": res["count_error"] == 0}
    res, cols = section_run(scene, fluid)
    view = pw.network_frames(res, cols, frames=res.get("gif_frames", pw.SECTION_GIF_FRAMES), rows=pw.SECTION_GIF_ROWS, fill=fill)
    lay = (pw.kymo_layout if kymo else pw.section_layout)(view, scale)
    return (pw.kymo_gif_bytes if kymo else pw.section_gif_bytes)(view, scale=scale), {
        "frames": len(view["frames"]), "canvas": f"{lay['width']}x{lay['height']}", "kind": "kymo" if kymo else "section",
        "columns": view["columns"], "t": [view["times"][0], view["times"][-1]],
        "conserved": res["conservation_error_q"] == 0}


def demo_ids() -> list[str]:
    return [f"D{n:02d}" for n in range(1, len(DEMOS) + 1)]


def demo_arg(text: str) -> str:
    """--only: 'D01'・'d1'・場面の名前・場面の絵文字（🏛️ など）のどれでも D<NN> に直す。"""
    k = _plain(text)
    scene = SCENE_ICONS.get(k, k)
    for did, (name, *_rest) in zip(demo_ids(), DEMOS):
        if scene == name or k.upper() == did or (k[:1] in "dD" and k[1:].isdigit() and int(k[1:]) == int(did[1:])):
            return did
    return text


def write_demos(out_dir: str = "docs/demos", *, check: bool = False, scale: int = 8, file_prefix: str = "",
                file_suffix: str = "", only=None) -> dict:
    """DEMOS を D01〜 の GIF に書く（名前は <prefix>D<NN>_<場面>_<流体><suffix>.gif、README から貼る）。

    隣の demos.json に何をどう描いたかと bytes・sha256（only で一部だけ書いたときは、その分だけ差し替える）。
    check=True なら書かずに、今の GIF とバイト単位で同じかだけ確かめる（README のデモがコードとずれていないか）。
    only は D<NN> の列（既定は全部）。Node から ``--call write_demos`` で呼べる。
    """
    out = Path(out_dir)
    wanted = None if only is None else {demo_arg(k) for k in only}
    items = []
    for did, (scene, fluid, fill, *how) in zip(demo_ids(), DEMOS):
        if wanted is not None and did not in wanted:
            continue
        n, kymo = int(did[1:]), how == ["kymo"]
        data, info = demo_gif(scene, fluid, fill, scale=scale, kymo=kymo)
        name = f"{file_prefix}D{n:02d}_{scene}{'-kymo' if kymo else ''}_{fluid}{file_suffix}.gif"
        path = out / name
        status = "written"
        if check:
            status = "missing" if not path.exists() else "same" if path.read_bytes() == data else "differs"
        else:
            out.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
        items.append({"id": f"D{n:02d}", "scene": scene, "icon": ("🕰️" if kymo else "") + scene_icon(scene),
                      "fluid": fluid, "fill": fill,
                      "file": name, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest(),
                      "status": status, **info})
    if not check:
        mpath = out / f"{file_prefix}demos{file_suffix}.json"
        kept = {}
        if mpath.exists() and wanted is not None:
            kept = {d["id"]: d for d in json.loads(mpath.read_text(encoding="utf-8")).get("demos", [])}
        kept.update({it["id"]: {k: v for k, v in it.items() if k != "status"} for it in items})
        manifest = {"schema": DEMOS_SCHEMA, "scale": scale, "demos": [kept[k] for k in sorted(kept)]}
        mpath.write_text(json.dumps(jsonable(manifest), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    ok = all(it["conserved"] for it in items) and (not check or all(it["status"] == "same" for it in items))
    return jsonable({"ok": ok, "out_dir": str(out), "demos": items})


def _demos(args) -> int:
    out = write_demos(args.out_dir, check=args.check, scale=args.scale, file_prefix=args.file_prefix,
                      file_suffix=args.file_suffix, only=args.only)
    if args.json:
        print(json.dumps(out, ensure_ascii=False, indent=2))
        return 0 if out["ok"] else 1
    marks = {"written": "🎞️", "same": "✅", "differs": "❌", "missing": "❓"}
    for it in out["demos"]:
        print(f"{marks[it['status']]} {it['id']} {it['icon']} {it['scene']:<16} {pw.EMOJI.get(it['fluid'], '🟦')} "
              f"{it['file']}  📐 {it['canvas']}  🎬 {it['frames']}  📦 {it['bytes'] / 1024:.1f} KB")
    return 0 if out["ok"] else 1


def _positive(text: str) -> float:
    try:
        x = float(text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("正の数") from exc
    if not math.isfinite(x) or x <= 0:
        raise argparse.ArgumentTypeError("正の数")
    return x


def build_parser() -> argparse.ArgumentParser:
    """引数は argparse の choices / type で宣言する（vendor/cli_args.py の make_parser。省略形は受け付けない）。
    実験・場面・流体は名前のほか絵文字をそのまま書ける: ``run 🥢``、``draw 🏛️ 🟩``、``draw 🥢 🟥``。"""
    cli = vendored("cli_args")
    p = cli.make_parser(suggest_unknown=True, prog="fluidfullspeedsimulation", description="🌊 流れ・遅れ・検証")
    sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("list", help="🗂️ 実験の一覧", allow_abbrev=False)
    plasma = sub.add_parser("rna-plasma", help="血漿中のRNA濃度・構成比・残存率の時系列をJSONで出力", allow_abbrev=False)
    plasma.add_argument("--added-percent", type=float, nargs="+", default=[1, 10, 50], help="投入直後のB/(A+B) [%%]")
    plasma.add_argument("--blood-litres", type=_positive, default=5.0, help="成人全血量 [L]")
    plasma.add_argument("--plasma-percent", type=_positive, default=55.0, help="全血に占める血漿の割合 [%%]")
    sub.add_parser("fluids", help="🧪 流体の物性", allow_abbrev=False)
    run = sub.add_parser("run", help="▶️ 実験を回す", allow_abbrev=False)
    run.add_argument("experiment", type=experiment_arg, choices=sorted(EXPERIMENTS),
                     help="ID・識別子・アイコン（例 E10 / somen-flume / 🥢）")
    ver = sub.add_parser("verify", help="✅ 全実験を回して検証", allow_abbrev=False)
    drw = sub.add_parser("draw", help="🎨 🟦 を積む断面図で溜まる・流れる・遅れるを見る（ドット絵・GIF・PNG も）",
                         allow_abbrev=False)
    scenes = RIG_SCENES + pw.SECTION_CASES + pw.SCENES
    drw.add_argument("scene", type=scene_arg, choices=scenes,
                     help="場面（名前か 🏛️ 🥢 🍜 🧱 🍽️ 🧊 🪵 🌧️ 🎍 🍷 🫗 🏺 🚧 ⚗️ 👼 🎋 🫙 🌊 🕳️）")
    drw.add_argument("block", nargs="?", type=block_arg,
                     help="流体のブロック・アイコン（🟦 🟨 🟧 🟩 🟪 / 💧 🍶 🧴 🟢 ☠️）か、塗りたい好きな絵文字")
    drw.add_argument("--view", choices=("frames", "frame", "files", "kymo"), default="frames",
                     help="frames = 時刻ごとの断面図 / frame = 最後の 1 枚 / files = GIF・PNG・受領記録 JSON / "
                          "kymo = 🕰️ 時空図（横 = 列・縦 = 時刻）の GIF・PNG・JSON")
    drw.add_argument("--fluid", type=fluid_arg, choices=ALL_FLUID_KEYS, default=None,
                     help="作動流体（キーか絵文字。🏷️ 仮置きも選べる。既定 water）")
    drw.add_argument("--fill", help="水のブロック（色は自由。既定は流体ごと: 🟦 水 🟨 エタノール 🟧 メタノール "
                                    "🟩 Namekian Water 🟪 Poison Water）")
    drw.add_argument("--temp-c", type=float, help="🧊 水温 °C（🥢 somen-flume・🧊 cold-water・🪵 bamboo-wear。表の範囲内）")
    visc = drw.add_mutually_exclusive_group()
    visc.add_argument("--mu", type=_positive, help="🍯 粘度 mPa·s を直接（水は 1.0016、温度より後に効く）")
    visc.add_argument("--viscosity-scale", type=_positive, help="🍯 粘度を何倍にするか（例 5 でどろっと）")
    drw.add_argument("--frames", type=int, default=None,
                     help=f"並べる時刻の数（既定 6、断面図の files は {pw.SECTION_GIF_FRAMES}）")
    drw.add_argument("--rows", type=int, default=None,
                     help=f"断面図の縦の段数（既定 {pw.SECTION_ROWS}、files は {pw.SECTION_GIF_ROWS}）")
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
    dem = sub.add_parser("demos", help="🎞️ README のデモ GIF（場面ごとに 1 枚）を書く・--check で確かめる",
                         allow_abbrev=False)
    dem.add_argument("--out-dir", default="docs/demos", help="書き出し先（既定 docs/demos）")
    dem.add_argument("--check", action="store_true", help="書かずに、今の GIF とバイト単位で同じか確かめる（違えば 1）")
    dem.add_argument("--only", nargs="+", type=demo_arg, choices=demo_ids(), metavar="DEMO",
                     help="一部だけ（D01・場面の名前・🏛️ などの絵文字）")
    dem.add_argument("--scale", type=int, default=8, help="1 マスのピクセル数（既定 8）")
    dem.add_argument("--file-prefix", default="", help="ファイル名の頭に付ける文字列")
    dem.add_argument("--file-suffix", default="", help="ファイル名の末尾（拡張子の前）に付ける文字列")
    dem.add_argument("--json", action="store_true", help="結果を JSON で出す")
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
    if args.command == "rna-plasma":
        if args.plasma_percent > 100 or any(not 0 <= p < 100 for p in args.added_percent):
            parser.error("血漿割合は0より大きく100以下、投入割合は0以上100未満です")
        print(json.dumps(plasma_rna_report(args.added_percent, args.blood_litres, args.plasma_percent), ensure_ascii=False))
        return 0
    if args.command == "list":
        for e in list_experiments():
            print(f"{e['icon']} {e['id']} {e['slug']}  {e['title']}")
        return 0
    if args.command == "fluids":
        print(render_fluids())
        return 0
    if args.command == "draw":
        return _draw(args)
    if args.command == "demos":
        return _demos(args)
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
