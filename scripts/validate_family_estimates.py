"""Family-aware estimation vs. treating relatives as unrelated, on simulated families.

For each replicate, simulates families with known founder haplotypes and
compares three estimates with the true haplotype frequencies:

* founders only, using the phase the children reveal (estimate_families)
* founders only, unphased (estimate_families with use_family_phase=False)
* every family member treated as an unrelated person (estimate)

Error is the total variation distance, 0.5 * sum |estimate - truth|.

    uv run python scripts/validate_family_estimates.py --families 1000 --replicates 3
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import hlahaplotools as h

ROOT = Path(__file__).resolve().parent.parent


def tvd(est: h.Estimate, truth: dict[str, float]) -> float:
    e = dict(est.frequencies(cut=0.0, renormalize=False))
    return 0.5 * sum(abs(e.get(k, 0.0) - truth.get(k, 0.0)) for k in set(e) | set(truth))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--families", type=int, default=1000)
    ap.add_argument("--children", default="1-4")
    ap.add_argument("--missing-parent", type=float, default=0.1)
    ap.add_argument("--missing", default="DQB1:0.2", help="per-locus untyped rate")
    ap.add_argument("--replicates", type=int, default=3)
    args = ap.parse_args()

    release = json.loads((ROOT / "data" / "manifest.json").read_text())["imgthla"]["release"]
    nom = h.Nomenclature(str(ROOT / "data" / "derived" / release))
    print("| replicate | founders, family phase | founders, unphased | everyone as unrelated |")
    print("|---|---|---|---|")
    for rep in range(args.replicates):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            subprocess.run(
                [
                    sys.executable,
                    str(ROOT / "scripts" / "simulate.py"),
                    str(tmp),
                    "--format",
                    "families",
                    "--samples",
                    str(args.families),
                    "--children",
                    args.children,
                    "--missing-parent",
                    str(args.missing_parent),
                    "--missing",
                    args.missing,
                    "--seed",
                    str(11 + rep),
                ],
                check=True,
                capture_output=True,
            )
            truth = {}
            for line in (tmp / "truth.tsv").read_text().splitlines():
                hap, f = line.split("\t")
                truth[hap] = float(f)
            fam = h.read_family_table(tmp / "families.csv")
            phased, _ = h.estimate_families(fam, nom, "asItIs", epsilon=1e-7)
            unphased, _ = h.estimate_families(fam, nom, "asItIs", use_family_phase=False, epsilon=1e-7)
            everyone = h.estimate(fam.dataset, nom, "asItIs", epsilon=1e-7)
        print(
            f"| {rep + 1} | {tvd(phased, truth):.4f} ({phased.n_samples} people) "
            f"| {tvd(unphased, truth):.4f} ({unphased.n_samples}) "
            f"| {tvd(everyone, truth):.4f} ({everyone.n_samples}) |",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
