"""Restricted RDF profile for FFSS hub parameters; no runtime RDF dependency."""
import json
import math
import re

PREFIX = "urn:ffss:hub:"
FIELD = "urn:ffss:field:"
# JSON lexical form preserves ordered tuples (the only sequence type in HUB_FLOWS).
DATATYPE = "urn:ffss:parameter-json"
LINE = re.compile(r'<urn:ffss:hub:([a-z-]+)> <urn:ffss:field:([a-z_]+)> ("(?:[^"\\\\]|\\\\.)*")\^\^<urn:ffss:parameter-json> \\.')

def _value(value):
    if value is None or type(value) in (str, int):
        return value
    if type(value) is float and math.isfinite(value):
        return value
    if type(value) is tuple:
        return [_value(v) for v in value]
    raise ValueError("unsupported parameter type")

def _restore(value):
    if isinstance(value, list):
        return tuple(_restore(v) for v in value)
    return _value(value)

def to_ttl(models):
    lines = []
    for name, fields in sorted(models.items()):
        if not re.fullmatch(r"[a-z-]+", name) or not fields:
            raise ValueError("invalid model")
        for key, value in sorted(fields.items()):
            if not re.fullmatch(r"[a-z_]+", key):
                raise ValueError("invalid field")
            lexical = json.dumps(_value(value), ensure_ascii=True, allow_nan=False)
            # ASCII JSON escapes are also valid Turtle string escapes.
            literal = json.dumps(lexical, ensure_ascii=True)
            lines.append(f"<{PREFIX}{name}> <{FIELD}{key}> {literal}^^<{DATATYPE}> .")
    return "\n".join(lines) + "\n"

def from_ttl(text):
    models = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        match = LINE.fullmatch(line)
        if not match:
            raise ValueError("unsupported FFSS TTL profile")
        name, key, literal = match.groups()
        fields = models.setdefault(name, {})
        if key in fields:
            raise ValueError("duplicate parameter")
        fields[key] = _restore(json.loads(json.loads(literal)))
    if not models:
        raise ValueError("empty model")
    return models

if __name__ == "__main__":
    import argparse
    import sys
    from fluidfullspeedsimulation import HUB_FLOWS
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("export", "verify"))
    args = parser.parse_args()
    if args.mode == "export":
        sys.stdout.write(to_ttl(HUB_FLOWS))
    else:
        if from_ttl(sys.stdin.read()) != HUB_FLOWS:
            parser.exit(1, "model mismatch\n")
        print("OK")
