"""🧪 fluids — 作動流体（水・エタノール・メタノール）の物性と、水槽・管・樋の流れの閉じた式。

ここは閉じた式（float）だけ。整数シミュレーションは network.py が受け持ち、両者のずれ（δ）をそのまま出す。

- Torricelli（小孔・非粘性）: dh/dt = -(Cd·a/A)·√(2gh)。ρ・μ が消えるので流体によらない。
- Hagen–Poiseuille（細管・層流）: Q = π r⁴ ρ g h / (8 μ L) = (π r⁴ g / (8 ν L))·h。
  h は時定数 τ = 8 ν L A / (π r⁴ g) で指数的に減る。速さは動粘度 ν = μ/ρ だけで決まる。
- 層流の前提は Re = v·D/ν < 2300 で確かめる。入口・出口損失や温度変化は入れない。
- 樋（開水路）: Manning（乱流、粘度によらない）と層流の薄い流れ（ν に効く）。断面は長方形と半円（竹）。
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass


G = 9.80665  # 標準重力加速度 [m/s²]
RE_LAMINAR = 2300.0


@dataclass(frozen=True)
class Fluid:
    key: str
    name: str
    icon: str
    rho: float  # 密度 [kg/m³]
    mu: float  # 粘度 [Pa·s]
    mu_range: tuple[float, float]  # 文献間の幅 [Pa·s]
    temp_c: float
    sources: tuple[str, ...]
    placeholder: bool = False  # 仮置き（名前と色だけ。物性はまだ決めていない）

    @property
    def nu(self) -> float:
        """動粘度 [m²/s]。"""
        return self.mu / self.rho

    def to_dict(self) -> dict:
        d = asdict(self)
        d["nu"] = self.nu
        d["mu_range"] = list(self.mu_range)
        d["sources"] = list(self.sources)
        return d


FLUIDS: dict[str, Fluid] = {
    "water": Fluid(
        "water", "Water", "💧", 998.21, 1.0016e-3, (1.0016e-3, 1.002e-3), 20.0,
        ("Wikipedia: List of viscosities (1.0016 mPa·s, 20 °C)",
         "chem-casts water table (998.2 kg/m³, 1.002 mPa·s)"),
    ),
    "ethanol": Fluid(
        "ethanol", "Ethanol", "🍶", 789.45, 1.20e-3, (1.144e-3, 1.20e-3), 20.0,
        ("Wikipedia: Ethanol data page (Lange 0.78945 g/cm³ at 20 °C)",
         "Wikipedia infobox 1.2 mPa·s at 20 °C; CRC 44th ed table 1.144 mPa·s"),
    ),
    "methanol": Fluid(
        "methanol", "Methanol", "🧴", 791.1, 0.587e-3, (0.587e-3, 0.59e-3), 20.0,
        ("Flottweg viscosity table (0.587 mPa·s, 0.742 mm²/s at 20 °C; ρ = μ/ν)",),
    ),
}


# 🏷️ 仮置きの流体: 名前と色だけ先に置く（物性は水と同じ）。実験の流体一覧（REAL_FLUIDS）には入れない
_W = FLUIDS["water"]
for _key, _name, _icon in (("namekian-water", "Namekian Water", "🟢"), ("poison-water", "Poison Water", "☠️")):
    FLUIDS[_key] = Fluid(_key, _name, _icon, _W.rho, _W.mu, _W.mu_range, _W.temp_c,
                         ("仮置き: 物性は水と同じ（名前と色だけ。中身は後で決める）",), placeholder=True)
REAL_FLUIDS = tuple(k for k, f in FLUIDS.items() if not f.placeholder)


def get_fluid(key: str) -> Fluid:
    try:
        return FLUIDS[key.lower()]
    except KeyError:
        raise KeyError(f"unknown fluid {key!r}; choose from {sorted(FLUIDS)}") from None


def fluid_table() -> list[dict]:
    return [f.to_dict() for f in FLUIDS.values()]


# --- 閉じた式 -------------------------------------------------------------

def poiseuille_tau(nu: float, tank_area: float, pipe_radius: float, pipe_length: float, g: float = G) -> float:
    """細管で排水する水槽の時定数 τ [s]。"""
    return 8.0 * nu * pipe_length * tank_area / (math.pi * pipe_radius ** 4 * g)


def poiseuille_head(h0: float, t: float, tau: float) -> float:
    return h0 * math.exp(-t / tau)


def pipe_mean_velocity(nu: float, head: float, pipe_radius: float, pipe_length: float, g: float = G) -> float:
    """細管の平均流速 v = r² g h / (8 ν L) [m/s]。"""
    return pipe_radius ** 2 * g * head / (8.0 * nu * pipe_length)


def reynolds(nu: float, velocity: float, diameter: float) -> float:
    return velocity * diameter / nu


def transport_delay(nu: float, head: float, pipe_radius: float, pipe_length: float, g: float = G) -> float:
    """水が細管を通り抜けるまでの遅れ L / v [s]（栓流近似）。"""
    return pipe_length / pipe_mean_velocity(nu, head, pipe_radius, pipe_length, g)


def torricelli_empty_time(h0: float, tank_area: float, orifice_area: float, cd: float, g: float = G) -> float:
    return tank_area / (cd * orifice_area) * math.sqrt(2.0 * h0 / g)


def torricelli_head(h0: float, t: float, tank_area: float, orifice_area: float, cd: float, g: float = G) -> float:
    s = math.sqrt(h0) - cd * orifice_area / tank_area * math.sqrt(g / 2.0) * t
    return s * s if s > 0 else 0.0


# --- 開水路（樋・竹の半円） ---------------------------------------------------

class Section:
    """樋の断面。rect（幅 size の長方形）と semicircle（半径 size の半円＝割った竹）。"""

    def __init__(self, shape: str, size: float):
        if shape not in ("rect", "semicircle"):
            raise ValueError("section shape must be 'rect' or 'semicircle'")
        if size <= 0:
            raise ValueError("section size must be positive")
        self.shape, self.size = shape, size

    def to_dict(self) -> dict:
        return {"shape": self.shape, "size": self.size}

    def area(self, h: float) -> float:
        if h <= 0:
            return 0.0
        if self.shape == "rect":
            return self.size * h
        r = self.size
        h = min(h, r)
        th = 2.0 * math.acos(1.0 - h / r)
        return r * r * (th - math.sin(th)) / 2.0

    def perimeter(self, h: float) -> float:
        if h <= 0:
            return 0.0
        if self.shape == "rect":
            return self.size + 2.0 * h
        r = self.size
        return r * 2.0 * math.acos(1.0 - min(h, r) / r)

    def top_width(self, h: float) -> float:
        """水面の幅。半円は 2√(h(2r − h))。"""
        if h <= 0:
            return 0.0
        if self.shape == "rect":
            return self.size
        r = self.size
        h = min(h, r)
        return 2.0 * math.sqrt(h * (2.0 * r - h))

    def max_area(self, max_depth: float | None = None) -> float:
        if self.shape == "semicircle":
            return math.pi * self.size ** 2 / 2.0
        if max_depth is None:
            raise ValueError("rect section needs max_depth")
        return self.size * max_depth

    def depth(self, area: float) -> float:
        """面積 → 水深。長方形は割り算、半円は二分法。"""
        if area <= 0:
            return 0.0
        if self.shape == "rect":
            return area / self.size
        lo, hi = 0.0, self.size
        if area >= self.area(hi):
            return hi
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if self.area(mid) < area else (lo, mid)
        return 0.5 * (lo + hi)


def hydraulic_radius(section: Section, h: float) -> float:
    p = section.perimeter(h)
    return section.area(h) / p if p > 0 else 0.0


def channel_flow(section: Section, h: float, slope: float, law: str, *, nu: float, n_manning: float = 0.010,
                 beta: float | None = None, g: float = G) -> float:
    """樋の流量 Q [m³/s]。

    - manning（乱流・粗さ支配）: Q = A·R^(2/3)·√S / n。粘度によらない
    - laminar（層流の薄い流れ）: v = g·S·R² / (β·ν)。β は長方形の広い流れで 3、半円（半分満ちた管）で 2
    """
    a = section.area(h)
    if a <= 0:
        return 0.0
    r = hydraulic_radius(section, h)
    if law == "manning":
        return a * r ** (2.0 / 3.0) * math.sqrt(slope) / n_manning
    if law == "laminar":
        if beta is None:
            beta = 3.0 if section.shape == "rect" else 2.0
        return a * g * slope * r * r / (beta * nu)
    raise ValueError("channel law must be 'manning' or 'laminar'")


def normal_depth(section: Section, q: float, slope: float, law: str, *, nu: float, h_max: float,
                 n_manning: float = 0.010, beta: float | None = None, g: float = G) -> float:
    """一様な流れの水深（等流水深）を二分法で求める。"""
    lo, hi = 0.0, h_max
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if channel_flow(section, mid, slope, law, nu=nu, n_manning=n_manning, beta=beta, g=g) < q:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def channel_reynolds(section: Section, h: float, q: float, nu: float) -> float:
    """開水路のレイノルズ数 Re = v·R/ν（目安: < 500 層流、> 2000 乱流）。"""
    a = section.area(h)
    return 0.0 if a <= 0 else (q / a) * hydraulic_radius(section, h) / nu


# --- 樋の継ぎ目・曲がり（流しそうめん） ------------------------------------------

def fall_time(drop: float, g: float = G) -> float:
    """段差 drop [m] を自由落下する時間 √(2·drop/g)。"""
    return math.sqrt(2.0 * drop / g) if drop > 0 else 0.0


def lip_flow(section: Section, h: float, lip: float, cd: float = 0.62, g: float = G,
             width: float | None = None) -> float:
    """継ぎ目の段（次の竹の縁が lip だけ高い）を越える流量。幅 = 縁の高さでの水面幅の長方形堰。
    width を渡すとその幅で計算する（縁に麺が積もって堰の頂が上がっても、幅は元の縁のまま、とみなすとき）。"""
    over = h - lip
    if over <= 0:
        return 0.0
    w = section.top_width(lip) if width is None else width
    return cd * (2.0 / 3.0) * w * math.sqrt(2.0 * g) * over ** 1.5


def lip_head(section: Section, q: float, lip: float, cd: float = 0.62, g: float = G,
             width: float | None = None) -> float:
    """lip_flow の逆: 流量 q が段を越えるときの水深（縁 + 越流水深）。"""
    w = section.top_width(lip) if width is None else width
    k = cd * (2.0 / 3.0) * w * math.sqrt(2.0 * g)
    return lip + (q / k) ** (2.0 / 3.0)


def superelevation(velocity: float, top_width: float, radius: float, g: float = G) -> float:
    """曲がりで外側の水面が内側より高くなる量 Δh = v²·B / (g·R)。"""
    return velocity * velocity * top_width / (g * radius)


def side_spill(cell_length: float, h_outer: float, rim: float, cd: float = 0.62, g: float = G) -> float:
    """外側の水面が縁 rim を越えた分が、セルの長さ全体から横へあふれる流量（横越流堰）。"""
    over = h_outer - rim
    if over <= 0:
        return 0.0
    return cd * (2.0 / 3.0) * cell_length * math.sqrt(2.0 * g) * over ** 1.5


__all__ = [
    "FLUIDS", "Fluid", "G", "REAL_FLUIDS", "RE_LAMINAR", "fluid_table", "get_fluid",
    "pipe_mean_velocity", "poiseuille_head", "poiseuille_tau", "reynolds",
    "torricelli_empty_time", "torricelli_head", "transport_delay",
    "Section", "channel_flow", "channel_reynolds", "hydraulic_radius", "normal_depth",
    "fall_time", "lip_flow", "lip_head", "side_spill", "superelevation",
]
