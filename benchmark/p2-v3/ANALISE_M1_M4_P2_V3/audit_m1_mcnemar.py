#!/usr/bin/env python3
"""Rebuild the M1 paired block table and exact two-sided McNemar test.

Uses only the published post-unblinding pair outcomes in
RESULTADOS_M1_M4_P2_V3.json. No third-party packages are required.
"""
import argparse
import csv
import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
DEFAULT_INPUT = HERE / "RESULTADOS_M1_M4_P2_V3.json"
DEFAULT_OUTPUT = HERE / "M1_BLOCK_AUDIT.csv"


def exact_mcnemar(a_pass, b_pass):
    # b: A pass/B fail; c: A fail/B pass.
    b = sum(a and not b for a, b in zip(a_pass, b_pass))
    c = sum(not a and b for a, b in zip(a_pass, b_pass))
    n = b + c
    if n == 0:
        return b, c, 1.0
    tail = sum(math.comb(n, k) for k in range(min(b, c) + 1)) / (2 ** n)
    return b, c, min(1.0, 2.0 * tail)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    data = json.loads(args.input.read_text(encoding="utf-8"))
    pairs = data["pairs"]
    if len(pairs) != 8 or len({p["pair_id"] for p in pairs}) != 8:
        raise SystemExit("Expected exactly eight unique paired blocks")

    rows = []
    a_pass, b_pass = [], []
    for p in pairs:
        a, b = p["A_M1"], p["B_M1"]
        if not isinstance(a, bool) or not isinstance(b, bool):
            raise SystemExit(f"M1 outcomes must be booleans: {p['pair_id']}")
        a_pass.append(a)
        b_pass.append(b)
        rows.append({
            "pair_id": p["pair_id"],
            "A_candidate": p["A_candidate"],
            "A_M1": "PASS" if a else "FAIL",
            "B_candidate": p["B_candidate"],
            "B_M1": "PASS" if b else "FAIL",
            "discordant": "yes" if a != b else "no",
            "direction": "B" if not a and b else "A" if a and not b else "tie",
        })

    b, c, p_value = exact_mcnemar(a_pass, b_pass)
    if sum(a_pass) != 7 or sum(b_pass) != 8 or (b, c, p_value) != (0, 1, 1.0):
        raise SystemExit(f"Unexpected M1 result: A={sum(a_pass)}, B={sum(b_pass)}, b={b}, c={c}, p={p_value}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    print(f"M1 contingency (rows A, columns B): [[{sum(a and b for a,b in zip(a_pass,b_pass))}, {b}], [{c}, {sum(not a and not b for a,b in zip(a_pass,b_pass))}]]")
    print(f"Exact two-sided McNemar: b={b}, c={c}, discordant={b+c}, p={p_value:.1f}")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
