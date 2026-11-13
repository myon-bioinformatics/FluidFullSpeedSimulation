"""Independent analytic oracles; keep the default rounding policy."""
import math
import unittest
from fixedq import from_q
from network import Network

class SingleValveReview(unittest.TestCase):
    def circuit(self, dt, a, count=1, incoming=0, outgoing=0, mode='exact', inflows='uniform', rounding='q48'):
        net = Network(dt, outlets=mode, inflows=inflows, rounding=rounding)
        net.tank('a', 1, volume=8).tank('sink', 1).tank('feed', 1, volume=128)
        for _ in range(count):
            net.valve('a', 'sink', tau=1, opening=a/count)
        if incoming:
            net.pump('feed', 'a', rate=incoming)
        if outgoing:
            net.pump('a', 'sink', rate=outgoing)
        return net

    def test_analytic_matrix(self):
        for a in (0, .25, .5, .75, 1):
            for dt in (1/8, 1/64):
                for count in (1, 2):
                    for incoming, outgoing in ((0, 0), (.25, 0), (0, .125), (.25, .125)):
                        with self.subTest(a=a, dt=dt, count=count, incoming=incoming, outgoing=outgoing):
                            net = self.circuit(dt, a, count, incoming, outgoing)
                            self.assertEqual(net.rounding, 'q48')
                            r = net.run(2)
                            want = 8+(incoming-outgoing)*2 if not a else 8*math.exp(-2*a)+(incoming-outgoing)*(-math.expm1(-2*a))/a
                            self.assertAlmostEqual(from_q(r['volumes_q']['a']), want, delta=2e-10)
                            self.assertEqual(r['conservation_error_q'], 0)
                            self.assertTrue(all(v >= 0 for v in r['volumes_q'].values()))

    def test_legacy_is_not_euler(self):
        dt, a, n = 1/8, .25, 16
        r = self.circuit(dt, a, mode='independent').run(n*dt)
        want = 8*(1-a*(-math.expm1(-dt)))**n
        self.assertAlmostEqual(from_q(r['volumes_q']['a']), want, delta=2e-10)
        self.assertGreater(want, 8*math.exp(-a*n*dt))
        self.assertLess(8*(1-a*dt)**n, 8*math.exp(-a*n*dt))

    def test_end_inflow(self):
        dt, a, incoming, n = 1/8, .25, .25, 16
        r = self.circuit(dt, a, incoming=incoming, inflows='end').run(n*dt)
        decay = math.exp(-a*dt)
        want = 8*decay**n+incoming*dt*(1-decay**n)/(1-decay)
        self.assertAlmostEqual(from_q(r['volumes_q']['a']), want, delta=2e-10)
        self.assertEqual(r['conservation_error_q'], 0)

    def test_rounding_modes(self):
        for rounding in ('q48', 'carry64'):
            r = self.circuit(1/8, .25, rounding=rounding).run(2)
            self.assertAlmostEqual(from_q(r['volumes_q']['a']), 8*math.exp(-.5), delta=2e-10)
            self.assertEqual(r['conservation_error_q'], 0)

    def test_starvation(self):
        r = self.circuit(1/8, .5, outgoing=100).run(2)
        self.assertEqual(r['conservation_error_q'], 0)
        self.assertTrue(all(v >= 0 for v in r['volumes_q'].values()))

if __name__ == '__main__':
    unittest.main()
