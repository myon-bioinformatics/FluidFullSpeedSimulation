import random
import unittest
from fractions import Fraction

from fixedq import (
    FRAC_BITS, ONE, Cam, DelayLine, from_q, lag_gain, lag_step, qdiv, qmul, qsqrt, random_scales, shift_lag_step,
    take, to_q,
)


class QArithmetic(unittest.TestCase):
    def test_round_trip_and_exact_decimal(self):
        self.assertEqual(to_q(1), ONE)
        self.assertEqual(to_q("0.5"), ONE // 2)
        self.assertEqual(to_q(Fraction(1, 3)), round(Fraction(ONE, 3)))
        self.assertAlmostEqual(from_q(to_q(0.123456)), 0.123456, places=12)
        with self.assertRaises(TypeError):
            to_q(True)

    def test_qmul_floors_and_qdiv_inverts(self):
        a, b = to_q(3.25), to_q(-1.5)
        self.assertEqual(qmul(a, b), to_q(-4.875))
        self.assertEqual(qmul(-1, 1), -1)  # 床: −2^-48 は −1 ulp に落ちる
        self.assertEqual(qdiv(to_q(1), to_q(4)), to_q(0.25))
        with self.assertRaises(ZeroDivisionError):
            qdiv(1, 0)

    def test_qsqrt(self):
        self.assertEqual(qsqrt(to_q(9)), to_q(3))
        self.assertLessEqual(abs(from_q(qsqrt(to_q(2))) - 2 ** 0.5), 2.0 ** -FRAC_BITS * 2)
        with self.assertRaises(ValueError):
            qsqrt(-1)

    def test_take_is_exact_for_random_amounts(self):
        rng = random.Random(7)
        for _ in range(2000):
            amount = rng.randrange(0, ONE * 10 ** 6)
            share = rng.randrange(0, ONE + 1)
            part, rest = take(amount, share)
            self.assertEqual(part + rest, amount)
            self.assertGreaterEqual(part, 0)
            self.assertGreaterEqual(rest, 0)
        with self.assertRaises(ValueError):
            take(10, ONE + 1)


class Lags(unittest.TestCase):
    def test_shift_lag_matches_multiplying_lag(self):
        rng = random.Random(3)
        for _ in range(3000):
            y, x, s = rng.randrange(-ONE * 50, ONE * 50), rng.randrange(-ONE * 50, ONE * 50), rng.randint(0, 40)
            self.assertEqual(lag_step(y, x, ONE >> s), shift_lag_step(y, x, s))

    def test_lag_gain_is_exact_discretisation(self):
        for item in random_scales(11, 6, n_min=1, n_max=3):
            dt = item["value"]
            k = from_q(lag_gain(dt, 1.0))
            self.assertAlmostEqual(k, 1 - 2.718281828459045 ** -dt, places=12)
        with self.assertRaises(ValueError):
            lag_gain(0, 1)

    def test_delay_line_returns_input_after_steps(self):
        for steps in (0, 1, 5):
            line = DelayLine(steps, fill=0)
            out = [line.push(v) for v in range(1, 12)]
            self.assertEqual(out[steps:], list(range(1, 12 - steps)))
            self.assertEqual(line.in_transit, sum(range(12 - steps, 12)))
        with self.assertRaises(ValueError):
            DelayLine(-1)


class CamTable(unittest.TestCase):
    def test_interpolates_and_clamps(self):
        cam = Cam([[0, 0], [1, 10], [3, 30]])
        got = [from_q(cam(to_q(x))) for x in (-1, 0, 0.5, 1, 2, 3, 4)]
        self.assertEqual(got, [0, 0, 5, 10, 20, 30, 30])
        self.assertEqual(cam.to_list()[-1], [3.0, 30.0])

    def test_binary_search_matches_linear_scan(self):
        rng = random.Random(5)
        pts = sorted({(round(rng.uniform(0, 10), 6), round(rng.uniform(-5, 5), 6)) for _ in range(60)})
        xs_seen = set()
        pts = [p for p in pts if not (p[0] in xs_seen or xs_seen.add(p[0]))]
        cam = Cam(pts)
        qpts = cam.points
        for _ in range(500):
            x = to_q(rng.uniform(-1, 11))
            if x <= qpts[0][0]:
                ref = qpts[0][1]
            elif x >= qpts[-1][0]:
                ref = qpts[-1][1]
            else:
                ref = next(y0 + (y1 - y0) * (x - x0) // (x1 - x0)
                           for (x0, y0), (x1, y1) in zip(qpts, qpts[1:]) if x <= x1)
            self.assertEqual(cam(x), ref)

    def test_rejects_duplicate_x(self):
        with self.assertRaises(ValueError):
            Cam([[1, 1], [1, 2]])


class RandomScales(unittest.TestCase):
    def test_a_and_n_are_random_and_reproducible(self):
        a = random_scales(42, 50, n_min=2, n_max=7)
        self.assertEqual(a, random_scales(42, 50, n_min=2, n_max=7))
        self.assertNotEqual(a, random_scales(43, 50, n_min=2, n_max=7))
        for item in a:
            self.assertTrue(1.0 <= item["A"] < 10.0)
            self.assertTrue(2 <= item["N"] <= 7)
            self.assertAlmostEqual(item["value"], item["A"] * 10.0 ** -item["N"])
        self.assertGreater(len({i["N"] for i in a}), 1)
        self.assertGreater(len({i["A"] for i in a}), 10)


if __name__ == "__main__":
    unittest.main()
