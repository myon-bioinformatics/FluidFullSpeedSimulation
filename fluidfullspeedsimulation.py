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
from fluids import FLUIDS, RE_LAMINAR, Section, channel_reynolds, fluid_table, get_fluid, normal_depth
from fluids import poiseuille_head, poiseuille_tau, torricelli_empty_time, torricelli_head, transport_delay
from network import Network, critical_gain, inverse_table, jsonable, peak_to_peak, settling_time, window

FLUID_KEYS = tuple(FLUIDS)
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
    rows = [f"{f['icon']} {f['key']:<9} ρ={f['rho']:g} kg/m³  μ={f['mu'] * 1e3:g} mPa·s  ν={f['nu'] * 1e6:.4g} mm²/s"
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


def build_parser() -> argparse.ArgumentParser:
    choices = sorted(EXPERIMENTS) + sorted(SLUGS)
    p = argparse.ArgumentParser(prog="fluidfullspeedsimulation", description="🌊 流れ・遅れ・検証")
    sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("list", help="🗂️ 実験の一覧")
    sub.add_parser("fluids", help="🧪 流体の物性")
    run = sub.add_parser("run", help="▶️ 実験を回す")
    run.add_argument("experiment", choices=choices, help="ID か識別子")
    ver = sub.add_parser("verify", help="✅ 全実験を回して検証")
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
