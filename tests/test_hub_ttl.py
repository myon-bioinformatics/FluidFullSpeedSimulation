import unittest
from fluidfullspeedsimulation import HUB_FLOWS, hub_network
from hub_ttl import to_ttl, from_ttl

class HubTTLTest(unittest.TestCase):
    def test_roundtrip(self):
        ttl = to_ttl(HUB_FLOWS)
        restored = from_ttl(ttl)
        self.assertEqual(restored, HUB_FLOWS)
        self.assertEqual(to_ttl(restored), ttl)
        for name, fields in HUB_FLOWS.items():
            for key, value in fields.items():
                self.assertIs(type(restored[name][key]), type(value))

    def test_simulation(self):
        restored = from_ttl(to_ttl(HUB_FLOWS))
        for name, fields in restored.items():
            with self.subTest(name=name):
                original = hub_network(name).run(64)
                rebuilt = hub_network(name, **fields).run(64)
                self.assertEqual(original["volumes_q"], rebuilt["volumes_q"])
                self.assertEqual(original["links"], rebuilt["links"])
                self.assertEqual(rebuilt["conservation_error_q"], 0)

    def test_reject(self):
        for text in ("", "invalid", to_ttl(HUB_FLOWS)*2):
            with self.assertRaises(ValueError):
                from_ttl(text)
