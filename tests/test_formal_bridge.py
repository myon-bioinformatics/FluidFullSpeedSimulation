"""Concrete checks connecting the Lean accounting model to the Python engine.

These bounded checks are regression/counterexample tests, not formal proofs of
Python. Fractions create exact Q48 inputs, including values beyond the exact
integer precision of floats. No Lean executable or third-party package is needed.
"""

import itertools
import random
import unittest
from fractions import Fraction

from fixedq import ONE, DelayLine
from network import Network


def quantity(integer):
    """An exact physical quantity whose Q48 representation is integer."""
    return Fraction(integer, ONE)


class FormalModelBridge(unittest.TestCase):
    def test_bounded_single_transfers_match_the_model(self):
        # Exhaustive small states find empty-source, full-drain and cap errors;
        # large cases detect accidental conversion through a floating number.
        cases = itertools.chain(
            itertools.product(range(7), range(4), range(9)),
            ((ONE - 1, ONE + 1, ONE),
             (2 ** 80 + 17, 2 ** 70 + 3, 2 ** 79 + 1),
             (2 ** 80 + 17, 0, 2 ** 81)),
        )
        for source, destination, requested in cases:
            with self.subTest(source=source, destination=destination, requested=requested):
                net = Network(1.0, inflows="end")
                net.tank("source", 1.0, volume=quantity(source))
                net.tank("destination", 1.0, volume=quantity(destination))
                net.pump("source", "destination", rate=quantity(requested))
                moved, = net.step(0)
                self.assertEqual(moved, min(source, requested))
                self.assertEqual(net.nodes["source"]["v"], source - moved)
                self.assertEqual(net.nodes["destination"]["v"], destination + moved)
                self.assertEqual(net.total_q(), source + destination)
                self.assertGreaterEqual(net.nodes["source"]["v"], 0)
                self.assertGreaterEqual(net.nodes["destination"]["v"], 0)

    def test_simultaneous_outlets_preserve_integer_remainders(self):
        # In particular, 2 quanta shared among 3 equal requests must not vanish.
        for available in (0, 1, 2, 5, 13, ONE + 1, 2 ** 70 + 3):
            for requests in ((1, 1, 1), (0, 3, 7), (3, 5, 7),
                             (2 ** 70, 2 ** 69, 2 ** 68)):
                with self.subTest(available=available, requests=requests):
                    net = Network(1.0, inflows="end")
                    net.tank("source", 1.0, volume=quantity(available))
                    for index, requested in enumerate(requests):
                        net.tank(str(index), 1.0)
                        net.pump("source", str(index), rate=quantity(requested))
                    moved = net.step(0)
                    expected_total = min(available, sum(requests))
                    self.assertEqual(sum(moved), expected_total)
                    self.assertEqual(net.nodes["source"]["v"], available - expected_total)
                    self.assertEqual(sum(net.nodes[str(i)]["v"] for i in range(3)), expected_total)
                    self.assertEqual(net.total_q(), available)
                    self.assertTrue(all(type(value) is int and value >= 0 for value in moved))

    def test_queue_push_conservation_and_order(self):
        inputs = [0, 1, ONE - 1, 2 ** 80 + 17, 2, 0, 5]
        for steps, fill in itertools.product((0, 1, 2, 5), (0, 3, 2 ** 70 + 1)):
            with self.subTest(steps=steps, fill=fill):
                line = DelayLine(steps, fill=fill)
                expected = ([fill] * steps + inputs)[:len(inputs)]
                outputs = []
                for sent in inputs:
                    before = line.in_transit
                    released = line.push(sent)
                    outputs.append(released)
                    self.assertEqual(line.in_transit + released, before + sent)
                    self.assertGreaterEqual(line.in_transit, 0)
                    self.assertGreaterEqual(released, 0)
                self.assertEqual(outputs, expected)

    def test_delayed_link_counts_queue_and_pending_separately(self):
        for delay in (1, 2, 5):
            with self.subTest(delay=delay):
                net = Network(1.0, inflows="end")
                net.tank("source", 1.0, volume=quantity(7))
                net.tank("destination", 1.0)
                link = net.pump("source", "destination", rate=quantity(7), until=1.0, delay=delay)
                saw_pending = False
                saw_queued = False
                for step in range(delay + 2):
                    source = net.nodes["source"]["v"]
                    destination = net.nodes["destination"]["v"]
                    queued = link["line"].in_transit
                    pending = link["pending"]
                    sent, = net.step(step)
                    released = link["pending"]
                    self.assertEqual(net.nodes["source"]["v"], source - sent)
                    self.assertEqual(net.nodes["destination"]["v"], destination + pending)
                    self.assertEqual(link["line"].in_transit, queued + sent - released)
                    self.assertEqual(net.nodes["destination"]["v"], 7 if step >= delay else 0)
                    self.assertEqual(net.total_q(), 7)
                    self.assertEqual(
                        sum(node["v"] for node in net.nodes.values()) +
                        link["line"].in_transit + link["pending"],
                        7,
                    )
                    saw_pending |= link["pending"] > 0
                    saw_queued |= link["line"].in_transit > 0
                self.assertTrue(saw_pending)
                self.assertEqual(saw_queued, delay > 1)

    def test_seeded_networks_compose_transfers_without_losing_stock(self):
        rng = random.Random(2718)
        for case in range(20):
            with self.subTest(case=case):
                net = Network(1.0, inflows="uniform" if case % 2 else "end",
                              rounding="carry64" if case % 3 == 0 else "q48")
                initial = [rng.randrange(51) for _ in range(5)]
                if case == 0:
                    initial[0] = 2 ** 80 + 17
                for index, value in enumerate(initial):
                    net.tank(str(index), 1.0, volume=quantity(value))
                # Always include a cycle; extra edges add competing outlets.
                pairs = [(i, (i + 1) % 5) for i in range(5)]
                pairs += [tuple(rng.sample(range(5), 2)) for _ in range(5)]
                for source, destination in pairs:
                    net.pump(str(source), str(destination), rate=quantity(rng.randrange(101)),
                             delay=rng.randrange(4))
                for step in range(24):
                    net.step(step)
                    # Sum actual storage directly, independently of total_q().
                    stocks = [node["v"] for node in net.nodes.values()]
                    stocks += [link["pending"] for link in net.links]
                    stocks += [link["line"].in_transit for link in net.links if link["line"] is not None]
                    self.assertTrue(all(type(value) is int and value >= 0 for value in stocks))
                    self.assertEqual(sum(stocks), sum(initial))
                    self.assertEqual(net.total_q(), sum(initial))


if __name__ == "__main__":
    unittest.main()
