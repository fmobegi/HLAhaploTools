"""haplo.stats (haplo.em) on cases already run by scripts/benchmark.py.

Each case folder from `benchmark.py run` or `extended` holds the input
(input.dat), the reference estimate (freqs_reference.tsv) and the other tools'
results (case.json). This runs haplo.em on the same input, scores it with the
benchmark's own log-likelihood and adds it to a table beside the stored tools.
A missing-typing case (--missing) is scored against the true frequencies only.

Usually run through `benchmark.py haplostats` (all of it, see run_all); by hand:

    python3 scripts/bench_haplostats.py --out bench/haplostats \\
        --cases bench/paper/cases/sim_n1000_2loci,bench/paper/cases/real_2f_B_C_n5000 \\
        [--missing bench/paper/sim_ext/n2000_l10_mic_missing] [--n-try 10] \\
        [--timeout 1800] [--mem-gb 6] [--rscript Rscript]

Needs R with haplo.stats (R_LIBS_USER may point at a user library) and Linux
(peak memory comes from getrusage; the memory cap uses prlimit).
"""

from __future__ import annotations

import argparse
import json
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark import (  # noqa: E402
    ROOT,
    canon,
    compare,
    log_likelihood,
    prepare_real,
    read_hfs,
    read_mac,
    run_ours,
    save_freqs,
    timed,
    write_mac,
)

R_SCRIPT = Path(__file__).resolve().parent / "haplostats_em.R"
RSCRIPT = ["Rscript"]  # command that runs R; set from --rscript (e.g. "module load R && Rscript" is not a command)


def run_haplostats(work: Path, data: Path, n_try: int, seed: int, timeout: float, mem_gb: float | None) -> dict:
    out = work / "haplostats"
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    cmd = [*RSCRIPT, str(R_SCRIPT), str(data), str(out), str(n_try), str(seed)]
    if mem_gb and shutil.which("prlimit"):
        cmd = ["prlimit", f"--as={int(mem_gb * 2**30)}", *cmd]
    r = timed(cmd, out, timeout)
    if r["status"] == "ok" and (out / "hfs.tsv").exists():
        r["freqs"] = read_hfs(out / "hfs.tsv")
        meta = dict(l.split("\t") for l in (out / "run.tsv").read_text().splitlines() if l)
        r["reported_logl"] = float(meta["lnlike"])
        r["converged"] = meta["converge"] == "1"
        r["subjects"] = int(meta["subjects"])
    elif r["status"] == "ok":
        r["status"] = "no output"
    if r["status"] != "ok":
        # Last R error line, e.g. memory or max.haps.limit.
        err = [l for l in r.get("stderr", "").splitlines() if l.strip() and not l.startswith("BENCH")]
        r["error"] = " | ".join(err[-3:])[-500:]
    r.pop("stdout", None)
    if r["status"] == "ok":
        r.pop("stderr", None)
    return r


def read_truth(path: Path) -> dict[str, float]:
    out = {}
    for line in path.read_text().splitlines():
        if line.strip():
            h, f = line.split("\t")[:2]
            out[canon(h)] = out.get(canon(h), 0.0) + float(f)
    return out


def complete_case(case_dir: Path, args) -> dict:
    loci, rows = read_mac(case_dir / "input.dat")
    stored = case_dir / "case.json"
    if stored.exists():
        case = json.loads(stored.read_text())
    else:  # an interrupted benchmark run: reference only, no other tools
        case = {"tools": {}, "reference": {"logl": log_likelihood(rows, read_hfs(case_dir / "freqs_reference.tsv"))}}
    work = Path(args.out).resolve() / case_dir.name
    work.mkdir(parents=True, exist_ok=True)
    print(f"\n== {case_dir.name}: {len(rows)} individuals, {'~'.join(loci)}", flush=True)
    r = run_haplostats(work, (case_dir / "input.dat").resolve(), args.n_try, args.seed, args.timeout, args.mem_gb)
    freqs = r.pop("freqs", None)
    ref = read_hfs(case_dir / "freqs_reference.tsv")
    if freqs is not None:
        save_freqs(work / "freqs_haplo_stats.tsv", freqs)
        total = sum(freqs.values())
        r["haplotypes"] = len(freqs)
        r["vs_reference"] = compare(freqs, ref)
        r["logl"] = log_likelihood(rows, {h: f / total for h, f in freqs.items()})
        r["delta_logl"] = r["logl"] - case["reference"]["logl"]
        truth = case_dir / "freqs_truth.tsv"
        if truth.exists():
            r["vs_truth"] = compare(freqs, read_truth(truth))
    out = {"case": case_dir.name, "loci": loci, "n": len(rows), "haplo.stats": r,
           "stored": {k: {f: v.get(f) for f in ("status", "delta_logl", "seconds", "peak_mb")}
                      for k, v in case["tools"].items()}}  # fmt: skip
    (work / "haplostats.json").write_text(json.dumps(out, indent=2) + "\n")
    show(r)
    return out


def missing_case(sim: Path, args) -> dict:
    """haplo.em on everyone, untyped alleles left missing; scored against the truth."""
    loci, rows = read_mac(sim / "population.dat")
    truth = read_truth(sim / "truth.tsv") if (sim / "truth.tsv").exists() else None  # none for real data
    name = f"missing_{sim.name}"
    work = Path(args.out).resolve() / name
    work.mkdir(parents=True, exist_ok=True)
    untyped = sum(1 for _, per in rows if not all(a and b for a, b in per))
    print(f"\n== {name}: {len(rows)} individuals ({untyped} incompletely typed), {'~'.join(loci)}", flush=True)
    r = run_haplostats(work, (sim / "population.dat").resolve(), args.n_try, args.seed, args.timeout, args.mem_gb)
    freqs = r.pop("freqs", None)
    # hlahaplotools on the same individuals: --missing-ours, or ours.json + freqs_hlahaplotools.tsv in the folder.
    ours = Path(args.missing_ours) if args.missing_ours else sim / "freqs_hlahaplotools.tsv"
    ours_freqs = read_hfs(ours) if ours.exists() else None
    if freqs is not None:
        save_freqs(work / "freqs_haplo_stats.tsv", freqs)
        r["haplotypes"] = len(freqs)
        if truth is not None:
            r["vs_truth"] = compare(freqs, truth)
        if ours_freqs is not None:
            r["vs_hlahaplotools"] = compare(freqs, ours_freqs)
    out = {"case": name, "loci": loci, "n": len(rows), "incomplete": untyped, "haplo.stats": r, "stored": {}}
    if ours_freqs is not None:
        s = json.loads((sim / "ours.json").read_text()) if (sim / "ours.json").exists() else {}
        if truth is not None:
            s["tvd_truth"] = compare(ours_freqs, truth)["tvd"]
        out["stored"]["hlahaplotools, everyone"] = s
    (work / "haplostats.json").write_text(json.dumps(out, indent=2) + "\n")
    show(r)
    return out


def show(r: dict) -> None:
    line = f"  haplo.stats  {r['status']}"
    if "seconds" in r:
        line += f"  {r['seconds']:.1f}s  {r.get('peak_mb', 0):.0f} MB"
    if r.get("delta_logl") is not None:
        line += f"  dlogL {r['delta_logl']:.3g}  TVD {r['vs_reference']['tvd']:.2e}"
    if "vs_truth" in r:
        line += f"  TVD to truth {r['vs_truth']['tvd']:.4f}"
    if "vs_hlahaplotools" in r:
        line += f"  TVD to hlahaplotools {r['vs_hlahaplotools']['tvd']:.4f}"
    if r.get("error"):
        line += f"  ({r['error']})"
    print(line, flush=True)


def report(results: list[dict], path: Path) -> None:
    lines = ["case\tn\ttool\tstatus\tdelta_logl\tseconds\tpeak_mb\ttvd_truth"]
    for c in results:
        r = c["haplo.stats"]
        lines.append("\t".join(map(str, [c["case"], c["n"], "haplo.stats", r["status"], r.get("delta_logl", ""),
                                         r.get("seconds", ""), r.get("peak_mb", ""),
                                         r.get("vs_truth", {}).get("tvd", "")])))  # fmt: skip
        for tool, s in c["stored"].items():
            cells = [c["case"], c["n"], tool, s.get("status", "ok"), s.get("delta_logl", ""), s.get("seconds", ""),
                     s.get("peak_mb", ""), s.get("tvd_truth", "")]  # fmt: skip
            lines.append("\t".join(map(str, cells)))
    path.write_text("\n".join(lines) + "\n")
    print(f"\nWrote {path}")


# ------------------------------------------------------------------ sets with untyped loci

C8 = ["A", "B", "C", "DRB1", "DQB1", "DPB1", "DQA1", "DPA1"]
SIM_SETS = {  # name: (loci, the two loci left untyped, share of individuals untyped at each)
    "l4_r30": (C8[:4], C8[2:4], 0.3),
    "l6_r30": (C8[:6], C8[4:6], 0.3),
    "l8_r30": (C8, C8[6:8], 0.3),
    "l10_r30": (C8 + ["MICA", "MICB"], ["MICA", "MICB"], 0.3),
    "l6_r10": (C8[:6], C8[4:6], 0.1),
    "l8_r10": (C8, C8[6:8], 0.1),
    "l10_r10": (C8 + ["MICA", "MICB"], ["MICA", "MICB"], 0.1),
}
REAL_SETS = {  # name: (loci, loci that may be untyped, as in the cohort)
    "DQA1_DQB1": (["DQA1", "DQB1"], {"DQA1"}),
    "DRB1_DQA1_DQB1": (["DRB1", "DQA1", "DQB1"], {"DQA1"}),
    "class2": (["DRB1", "DQA1", "DQB1", "DPA1", "DPB1"], {"DQA1", "DPA1"}),
    "classical8": (C8, {"DQA1", "DPA1"}),
}


def _ours(sim: Path, data: Path, loci, threads: int, timeout: float) -> None:
    """hlahaplotools (10 starts) on a set, next to it: ours.json and freqs_hlahaplotools.tsv."""
    r = run_ours(sim, data, loci, 1e-6, threads, timeout, starts=10)
    freqs = r.pop("freqs", None)
    r.pop("stdout", None)
    if freqs is not None:
        save_freqs(sim / "freqs_hlahaplotools.tsv", freqs)
        if (sim / "truth.tsv").exists():
            r["tvd_truth"] = compare(freqs, read_hfs(sim / "truth.tsv"))["tvd"]
        r.pop("stderr", None)
    (sim / "ours.json").write_text(json.dumps(r, indent=2) + "\n")
    print(
        f"  hlahaplotools {sim.name}: {r.get('status')} {r.get('seconds', 0):.1f}s {r.get('peak_mb', 0):.0f} MB",
        flush=True,
    )


def build_missing_sim(out: Path, n: int = 2000, seed: int = 1, threads: int = 8, timeout: float = 3600) -> list[Path]:
    """Simulated populations with two loci untyped in some individuals, and hlahaplotools on each."""
    dirs = []
    for name, (loci, miss, rate) in SIM_SETS.items():
        sim = out / name
        sim.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "simulate.py"), str(sim), "--samples", str(n),
             "--loci", ",".join(loci), "--haplotypes", "300", "--alleles", "25", "--seed", str(seed),
             "--missing", ",".join(f"{l}:{rate}" for l in miss)],
            check=True, capture_output=True,
        )  # fmt: skip
        _ours(sim, sim / "population.dat", loci, threads, timeout)
        dirs.append(sim)
    return dirs


def build_missing_real(real: Path, out: Path, n: int = 5000, seed: int = 1, threads: int = 8,
                       timeout: float = 3600) -> list[Path]:  # fmt: skip
    """Samples of the real cohort with its own untyped DQA1/DPA1, and hlahaplotools on each."""
    import hlahaplotools as h

    manifest = json.loads((ROOT / "data" / "manifest.json").read_text())
    nom = h.Nomenclature(str(ROOT / "data" / "derived" / manifest["imgthla"]["release"]))
    dirs = []
    for name, (loci, untyped) in REAL_SETS.items():
        rows, info = prepare_real(real, loci, "2f", nom, n, seed, may_be_untyped=untyped)
        sim = out / name
        sim.mkdir(parents=True, exist_ok=True)
        write_mac(sim / "population.dat", loci, rows)
        _ours(sim, sim / "population.dat", loci, threads, timeout)
        r = json.loads((sim / "ours.json").read_text())
        r["info"] = {**info, "incompletely_typed": sum(1 for _, per in rows if not all(a for a, _ in per))}
        (sim / "ours.json").write_text(json.dumps(r, indent=2) + "\n")
        dirs.append(sim)
    return dirs


def run_all(root: Path, real: Path | None, case_dirs: list[Path], *, n_try=10, seed=1, timeout=1800.0, mem_gb=6.0,
            retry_gb=12.0, threads=8) -> int:  # fmt: skip
    """Every haplo.stats run of the paper, in the folder layout tables_and_figures.py reads (under root):
    haplostats/ (fully typed cases), missing_typing/ and real_missing/ (sets with untyped loci, with
    hlahaplotools), haplostats_missing/ and haplostats_real_missing/ (haplo.stats at mem_gb), and
    haplostats_12gb/ (the runs that failed, again at retry_gb)."""
    args = argparse.Namespace(out=str(root / "haplostats"), n_try=n_try, seed=seed, timeout=timeout, mem_gb=mem_gb,
                              missing_ours="")  # fmt: skip
    cases = sorted(d for base in case_dirs if base.exists() for d in base.iterdir()
                   if (d / "input.dat").exists() and (d / "freqs_reference.tsv").exists())  # fmt: skip
    Path(args.out).mkdir(parents=True, exist_ok=True)
    done = [complete_case(c, args) for c in cases]
    (Path(args.out) / "haplostats.json").write_text(json.dumps(done, indent=2) + "\n")
    report(done, Path(args.out) / "haplostats.tsv")

    groups = [("haplostats_missing", build_missing_sim(root / "missing_typing", seed=seed, threads=threads))]
    if real:
        groups.append(("haplostats_real_missing", build_missing_real(real, root / "real_missing", seed=seed,
                                                                     threads=threads)))  # fmt: skip
    failed = []
    for name, dirs in groups:
        args.out = str(root / name)
        res = [missing_case(d, args) for d in dirs]
        (root / name / "haplostats.json").write_text(json.dumps(res, indent=2) + "\n")
        report(res, root / name / "haplostats.tsv")
        failed += [d for d, r in zip(dirs, res, strict=True) if r["haplo.stats"]["status"] != "ok"]
    if failed and retry_gb and retry_gb > mem_gb:
        args.out, args.mem_gb = str(root / "haplostats_12gb"), retry_gb
        res = [missing_case(d, args) for d in failed]
        (root / "haplostats_12gb" / "haplostats.json").write_text(json.dumps(res, indent=2) + "\n")
        report(res, root / "haplostats_12gb" / "haplostats.tsv")
    return 0


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--rscript", default="Rscript", help="command that runs R scripts (default Rscript)")
    p.add_argument("--out", required=True)
    p.add_argument("--cases", default="", help="comma-separated case folders from benchmark.py")
    p.add_argument("--missing", default="", help="comma-separated simulation folders with population.dat and truth.tsv")
    p.add_argument("--missing-ours", default="", help="hlahaplotools frequencies on everyone (freqs_everyone.tsv)")
    p.add_argument("--n-try", type=int, default=10, help="haplo.em random starts (its default is 10)")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--timeout", type=float, default=1800)
    p.add_argument("--mem-gb", type=float, default=6)
    args = p.parse_args()
    RSCRIPT[:] = shlex.split(args.rscript)
    Path(args.out).mkdir(parents=True, exist_ok=True)
    results = [complete_case(Path(c), args) for c in args.cases.split(",") if c]
    results += [missing_case(Path(m), args) for m in args.missing.split(",") if m]
    (Path(args.out) / "haplostats.json").write_text(json.dumps(results, indent=2) + "\n")
    report(results, Path(args.out) / "haplostats.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
