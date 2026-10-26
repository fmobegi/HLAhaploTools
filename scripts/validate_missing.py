"""Check the missing-data handling against simulated truth.

For each untyped rate, simulates a population, blanks loci at random, and
compares three estimates with the true haplotype frequencies:

* all samples, untyped loci handled by hlahaplotools  ("with missing")
* only fully typed samples, as Hapl-o-Mat would use    ("drop incomplete")
* the same samples before blanking                     ("no missing", best case)

Error is the total variation distance, 0.5 * sum |estimate - truth|, over
all haplotypes (0 = perfect, 1 = disjoint).

    uv run python scripts/validate_missing.py --samples 5000 --rates 0.05,0.2,0.4 --replicates 3
"""

from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

import hlahaplotools as h

ROOT = Path(__file__).resolve().parent.parent


def read_truth(path: Path) -> dict[str, float]:
    return {hap: float(f) for hap, f in (line.split("\t") for line in path.read_text().splitlines())}


def tvd(est: list[tuple[str, float]], truth: dict[str, float]) -> float:
    e = dict(est)
    return 0.5 * sum(abs(e.get(k, 0.0) - truth.get(k, 0.0)) for k in set(e) | set(truth))


def run(path: Path, nom: h.Nomenclature) -> tuple[h.Estimate, int]:
    """Estimate at the typed (2-field) resolution; return the estimate and sample count."""
    ds = h.read_mac(str(path), half_typed="unknown")  # --half drops alleles at random
    return h.estimate(ds, nom, "asItIs", epsilon=1e-7), len(ds)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", type=int, default=5000)
    ap.add_argument("--loci", default="A,B,C,DRB1,DQB1")
    ap.add_argument("--blank", default="C,DQB1", help="loci to blank at random (default C,DQB1)")
    ap.add_argument("--rates", default="0.05,0.2,0.4")
    ap.add_argument("--half-rate", type=float, default=0.0, help="also make one allele untyped at this rate")
    ap.add_argument("--replicates", type=int, default=3)
    ap.add_argument("--out", type=Path, help="write a markdown table here")
    args = ap.parse_args()

    manifest_release = json.loads((ROOT / "data" / "manifest.json").read_text())["imgthla"]["release"]
    nom = h.Nomenclature(str(ROOT / "data" / "derived" / manifest_release))

    rows = []
    header = (
        f"| untyped rate ({args.blank}) | samples kept if dropped | error: with missing | "
        f"error: drop incomplete | error: no missing |"
    )
    print(header)
    print("|---|---|---|---|---|")
    for rate in [float(r) for r in args.rates.split(",")]:
        errs = {"missing": [], "drop": [], "full": []}
        kept = []
        for rep in range(args.replicates):
            with tempfile.TemporaryDirectory() as tmp:
                tmp = Path(tmp)
                spec = ",".join(f"{l}:{rate}" for l in args.blank.split(","))
                cmd = [
                    sys.executable,
                    str(ROOT / "scripts" / "simulate.py"),
                    str(tmp),
                    "--samples",
                    str(args.samples),
                    "--loci",
                    args.loci,
                    "--missing",
                    spec,
                    "--seed",
                    str(100 + rep),
                ]
                if args.half_rate:
                    cmd += ["--half", ",".join(f"{l}:{args.half_rate}" for l in args.blank.split(","))]
                subprocess.run(cmd, check=True, capture_output=True)
                truth = read_truth(tmp / "truth.tsv")
                est_m, _ = run(tmp / "population.dat", nom)
                est_d, n_kept = run(tmp / "population_complete_only.dat", nom)
                est_f, _ = run(tmp / "population_full.dat", nom)
                errs["missing"].append(tvd(est_m.frequencies(cut=0.0, renormalize=False), truth))
                errs["drop"].append(tvd(est_d.frequencies(cut=0.0, renormalize=False), truth))
                errs["full"].append(tvd(est_f.frequencies(cut=0.0, renormalize=False), truth))
                kept.append(n_kept)
        fmt = lambda v: f"{statistics.mean(v):.4f} ± {statistics.stdev(v):.4f}" if len(v) > 1 else f"{v[0]:.4f}"  # noqa: E731
        row = (
            f"| {rate:.0%} | {statistics.mean(kept):.0f} / {args.samples} | {fmt(errs['missing'])} | "
            f"{fmt(errs['drop'])} | {fmt(errs['full'])} |"
        )
        print(row, flush=True)
        rows.append(row)
    if args.out:
        args.out.write_text("\n".join([header, "|---|---|---|---|---|", *rows]) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
