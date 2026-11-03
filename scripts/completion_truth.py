"""Linked versus exhaustive completion, judged against the true frequencies.

Simulated populations at the 8 classical loci with DQA1 and DPA1 untyped in some people.
Compared: exhaustive completion, linked completion, fully typed people only, and the
complete data before blanking (the best achievable). Several replicates.

    uv run python scripts/completion_truth.py [--n 3000] [--rate 0.3] [--reps 3] [--out bench/completion_truth.json]

Usually run through `benchmark.py completion`.
"""

import argparse
import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import hlahaplotools as h

ROOT = Path(__file__).resolve().parent.parent
LOCI = ["A", "B", "C", "DRB1", "DQB1", "DPB1", "DQA1", "DPA1"]


def tvd(est: dict, truth: dict) -> float:
    keys = set(est) | set(truth)
    return 0.5 * sum(abs(est.get(k, 0.0) - truth.get(k, 0.0)) for k in keys)


def canon(hap: str) -> str:
    return "~".join(sorted(hap.split("~"), key=lambda a: a.split("*")[0]))


def run(nom, path, completion):
    ds = h.read_mac(str(path), half_typed="unknown", drb345=False)
    t = time.perf_counter()
    e = h.estimate(ds, nom, "asItIs", starts=10, completion=completion, diplotypes=False, epsilon=1e-6)
    return {canon(k): v for k, v in e.frequencies(cut=0.0)}, time.perf_counter() - t, e


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=3000)
    ap.add_argument("--rate", type=float, default=0.3)
    ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--out", type=Path, default=ROOT / "bench" / "completion_truth.json")
    args = ap.parse_args()
    release = json.loads((ROOT / "data" / "manifest.json").read_text())["imgthla"]["release"]
    nom = h.Nomenclature(str(ROOT / "data" / "derived" / release))
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        for rep in range(1, args.reps + 1):
            d = Path(tmp) / f"r{rep}"
            subprocess.run(
                [sys.executable, str(ROOT / "scripts/simulate.py"), str(d), "--samples", str(args.n),
                 "--loci", ",".join(LOCI), "--haplotypes", "300", "--alleles", "20", "--seed", str(rep),
                 "--missing", f"DQA1:{args.rate},DPA1:{args.rate}"],
                check=True, capture_output=True,
            )  # fmt: skip
            truth = {}
            for line in (d / "truth.tsv").read_text().splitlines():
                k, f = line.split("\t")
                truth[canon(k)] = truth.get(canon(k), 0.0) + float(f)
            out = {"rep": rep}
            for name, path, comp in [
                ("complete data (before blanking)", d / "population_full.dat", "linked"),
                ("fully typed only", d / "population_complete_only.dat", "linked"),
                ("everyone, linked", d / "population.dat", "linked"),
                ("everyone, all", d / "population.dat", "all"),
            ]:
                f, secs, e = run(nom, path, comp)
                out[name] = {"tvd": tvd(f, truth), "seconds": secs, "people": e.n_samples,
                             "candidates": e.n_candidates, "logl": e.log_likelihood}  # fmt: skip
                print(f"rep {rep}  {name:32} TVD {out[name]['tvd']:.4f}  people {e.n_samples:5}  "
                      f"candidates {e.n_candidates:>9,}  {secs:6.1f}s", flush=True)  # fmt: skip
            rows.append(out)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(rows, indent=1) + "\n")
    print("mean TVD:")
    for name in [k for k in rows[0] if k != "rep"]:
        v = [r[name]["tvd"] for r in rows]
        print(f"  {name:32} {sum(v) / len(v):.4f}  (range {min(v):.4f}–{max(v):.4f})")


if __name__ == "__main__":
    main()
