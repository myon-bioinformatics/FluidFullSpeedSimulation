"""💧 fixedq — 整数の固定小数点（Q 形式）とビット演算で、水を 1 滴も失わずに流す。

- 量はすべて ``int``。``1.0 == ONE == 1 << FRAC_BITS``。
- 掛け算は ``(a * b) >> FRAC_BITS``（床）。誤差は 1 ulp 未満で、常に「少なめ」側に出る。
- ``take`` は取り分を切り出し、端数は残りへ回す。``part + rest == amount`` が常に成り立つので、
  水槽の総量（保存則）がビット単位で一致する。
- ``lag_step`` は 1 次遅れ（水槽 1 段）。``k = 2**-s`` なら ``(x - y) >> s`` のシフトだけになる。
- ``qsqrt`` は ``math.isqrt`` による整数平方根（Torricelli の √h 用）。

標準ライブラリのみ。Python 3.11+。
"""

from __future__ import annotations

import bisect
import math
import random
from collections import deque
from fractions import Fraction

FRAC_BITS = 48  # 1 ulp = 2^-48 ≈ 3.6e-15（m³ なら 3.6e-6 mm³）。床の偏りは n ステップで高々 n ulp
ONE = 1 << FRAC_BITS


def to_q(x) -> int:
    """数値（int・float・str・Fraction）を Q 整数へ。str は 10 進のまま厳密に読む。"""
    if isinstance(x, bool):
        raise TypeError("bool is not a quantity")
    f = Fraction(x) if not isinstance(x, str) else Fraction(x.strip())
    return round(f * ONE)


def from_q(q: int) -> float:
    return q / ONE


def qmul(a: int, b: int) -> int:
    """Q × Q → Q（床）。"""
    return (a * b) >> FRAC_BITS


def qdiv(a: int, b: int) -> int:
    """Q ÷ Q → Q（床）。"""
    if b == 0:
        raise ZeroDivisionError("qdiv by zero")
    return (a << FRAC_BITS) // b


def qsqrt(a: int) -> int:
    """√a（Q）。負は 0 とみなさず ValueError。"""
    if a < 0:
        raise ValueError("qsqrt of negative")
    return math.isqrt(a << FRAC_BITS)


def take(amount: int, share: int) -> tuple[int, int]:
    """amount から share（Q、0..ONE）の割合を取り出す。戻り値 (part, rest) は合計が amount と完全一致。"""
    if not 0 <= share <= ONE:
        raise ValueError("share must be within 0..ONE")
    part = qmul(amount, share)
    return part, amount - part


def clamp(x: int, lo: int, hi: int) -> int:
    return lo if x < lo else hi if x > hi else x


def lag_gain(dt: float, tau: float) -> int:
    """1 次遅れ dy/dt = (x - y)/tau の厳密離散化ゲイン k = 1 - exp(-dt/tau)（Q）。"""
    if dt <= 0 or tau <= 0:
        raise ValueError("dt and tau must be positive")
    return to_q(-math.expm1(-dt / tau))


def lag_step(y: int, x: int, k: int) -> int:
    """水槽 1 段の遅れ: y ← y + k(x - y)。"""
    return y + qmul(x - y, k)


def shift_lag_step(y: int, x: int, s: int) -> int:
    """k = 2**-s のときの lag_step。掛け算なしのビットシフトだけ。"""
    return y + ((x - y) >> s)


class DelayLine:
    """輸送の遅れ（パイプ内の水）。steps 回の push の後に、入れた量がそのまま出てくる。"""

    def __init__(self, steps: int, fill: int = 0):
        if steps < 0:
            raise ValueError("steps must be >= 0")
        self.steps = steps
        self._q = deque([fill] * steps)

    def push(self, x: int) -> int:
        if self.steps == 0:
            return x
        self._q.append(x)
        return self._q.popleft()

    @property
    def in_transit(self) -> int:
        return sum(self._q)


class Cam:
    """交換可能な関数カム（区分線形の表）。x の外側は端の値で止まる。整数だけで補間する（二分探索）。"""

    def __init__(self, points):
        pts = sorted((to_q(x), to_q(y)) for x, y in points)
        if len(pts) < 1:
            raise ValueError("cam needs at least one point")
        if any(pts[i][0] == pts[i + 1][0] for i in range(len(pts) - 1)):
            raise ValueError("cam x values must be distinct")
        self.points = pts
        self._xs = [x for x, _ in pts]

    def __call__(self, x: int) -> int:
        pts = self.points
        if x <= pts[0][0]:
            return pts[0][1]
        if x >= pts[-1][0]:
            return pts[-1][1]
        i = bisect.bisect_left(self._xs, x)
        (x0, y0), (x1, y1) = pts[i - 1], pts[i]
        return y0 + (y1 - y0) * (x - x0) // (x1 - x0)

    @property
    def x_max(self) -> int:
        return self.points[-1][0]

    def to_list(self) -> list[list[float]]:
        return [[from_q(x), from_q(y)] for x, y in self.points]


def random_scales(seed: int, count: int, n_min: int = 1, n_max: int = 6, digits: int = 3) -> list[dict]:
    """試す係数を A×10^-N で返す。A（1 ≤ A < 10、有効 digits 桁）と N（n_min..n_max）を乱数で選ぶ。

    10^-3 刻みのような人間都合の決め打ちを避けるための道具。seed で再現できる。
    """
    if count < 0 or n_min > n_max or digits < 1:
        raise ValueError("bad random_scales arguments")
    rng = random.Random(seed)
    out = []
    for _ in range(count):
        a = round(rng.uniform(1.0, 10.0), digits - 1)
        if a >= 10.0:
            a = 9.99
        n = rng.randint(n_min, n_max)
        out.append({"A": a, "N": n, "value": a * 10.0 ** -n, "label": f"{a}e-{n}"})
    return out
