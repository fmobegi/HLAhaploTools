"""Benchmark hlahaplotools against Hapl-o-Mat and PyPop.

All three tools get the *same* input: people fully typed at the chosen loci,
each typing a single allele name (no ambiguity, no missing data), used as
typed. Differences therefore come from the EM, not from how the tools read
nomenclature. Each case is scored against a reference: hlahaplotools run to a
tight tolerance (1e-10), i.e. the maximum-likelihood estimate.

Per tool and case: largest frequency difference and total variation distance
(TVD) to the reference, log-likelihood of the tool's frequencies (computed here,
independently of every tool, by summing over all phases of each person), the
tool's own reported log-likelihood, run time and peak memory. Simulated cases
are also scored against the true frequencies.

    uv run python scripts/benchmark.py run --out bench/results \\
        [--real typing.dat] [--sizes 1000,5000] [--threads 1,8]
    uv run python scripts/benchmark.py multistart --out bench/results --real typing.dat
    uv run python scripts/benchmark.py bootstrap --out bench/results --real typing.dat
    uv run python scripts/benchmark.py report --out bench/results

Needs: reference/Hapl-o-Mat/haplomat (built with `make`) and PyPop in its own
environment (`uv venv bench/pypop-env && uv pip install pypop-genomics`), or
pass --haplomat / --pypop. A tool that is missing, fails or exceeds its limits
is recorded as such; the others still run.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import platform
import random
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
HAPLOMAT = ROOT / "reference" / "Hapl-o-Mat" / "haplomat"
PYPOP = ROOT / "bench" / "pypop-env" / "bin" / "pypop"

# Locus sets for the real-data cases (all within PyPop's limits).
REAL_SETS = [
    ["B", "C"],
    ["A", "B", "C"],
    ["DRB1", "DQB1"],
    ["DQA1", "DQB1"],
    ["DPA1", "DPB1"],
    ["A", "B", "C", "DRB1"],
    ["A", "B", "C", "DRB1", "DQB1"],
    ["A", "B", "C", "DRB1", "DQB1", "DPB1"],
]
SIM_LOCI = ["A", "C", "B", "DRB1", "DQB1", "DPB1"]

# ------------------------------------------------------------------ data


def read_mac(path: Path) -> tuple[list[str], list[tuple[str, list[tuple[str, str]]]]]:
    """A MAC table as (loci, [(id, [(a, b) per locus])]); names get their locus prefix."""
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    header = lines[0].split("\t")
    loci = header[1::2]
    rows = []
    for line in lines[1:]:
        cells = line.split("\t")
        cells += [""] * (len(header) - len(cells))
        per = []
        for k, locus in enumerate(loci):
            a, b = cells[1 + 2 * k].strip(), cells[2 + 2 * k].strip()
            per.append((f"{locus}*{a}" if a else "", f"{locus}*{b}" if b else ""))
        rows.append((cells[0], per))
    return loci, rows


def write_mac(path: Path, loci: list[str], rows) -> None:
    with path.open("w") as f:
        f.write("id\t" + "\t".join(f"{l}\t{l}" for l in loci) + "\n")
        for sid, per in rows:
            f.write(sid + "\t" + "\t".join(f"{a.split('*')[-1]}\t{b.split('*')[-1]}" for a, b in per) + "\n")


def write_pop(path: Path, loci: list[str], rows) -> None:
    """PyPop genotype file: population header block, then id + LOCUS_1/LOCUS_2 columns."""
    with path.open("w") as f:
        f.write("populationname\nbench\n")
        f.write("id\t" + "\t".join(f"{l.lower()}_1\t{l.lower()}_2" for l in loci) + "\n")
        for sid, per in rows:
            f.write(sid + "\t" + "\t".join(f"{a.split('*')[1]}\t{b.split('*')[1]}" for a, b in per) + "\n")


def prepare_real(
    src: Path, loci: list[str], resolution: str, nom, n: int | None, seed: int, may_be_untyped=()
):
    """People fully typed at `loci`, each typing reduced to one name at `resolution`.

    A typing that maps to more than one name (a G group, a MAC code) leaves the
    person out, so every tool sees the same unambiguous data. Loci in
    `may_be_untyped` may instead be left untyped (both alleles blank).
    """
    all_loci, rows = read_mac(src)
    idx = [all_loci.index(l) for l in loci]
    cache: dict[str, str | None] = {}

    def one(name: str) -> str | None:
        if name not in cache:
            try:
                codes = {c for c, _ in nom.translate(name, resolution)}
            except ValueError:
                codes = set()
            cache[name] = codes.pop() if len(codes) == 1 else None
        return cache[name]

    kept, dropped = [], Counter()
    for sid, per in rows:
        out = []
        for k in idx:
            a, b = per[k]
            if not a and not b:
                if all_loci[k] in may_be_untyped:
                    out.append(("", ""))
                    continue
                dropped["untyped"] += 1
                break
            b = b or a  # blank second allele: homozygous (all tools get the copy)
            ra, rb = one(a), one(b)
            if ra is None or rb is None:
                dropped["ambiguous or unknown at " + resolution] += 1
                break
            out.append((ra, rb))
        else:
            kept.append((sid, out))
    total_kept = len(kept)
    if n is not None and len(kept) > n:
        kept = random.Random(seed).sample(kept, n)
    return kept, {"fully_typed_unambiguous": total_kept, "used": len(kept), "dropped": dict(dropped)}


def simulate(out: Path, loci: list[str], n: int, seed: int) -> tuple[list, dict[str, float]]:
    """Simulated population (fully typed) and its true haplotype frequencies."""
    out.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "simulate.py"), str(out), "--samples", str(n),
         "--loci", ",".join(loci), "--haplotypes", "300", "--alleles", "25", "--seed", str(seed)],
        check=True, capture_output=True,
    )  # fmt: skip
    _, rows = read_mac(out / "population.dat")
    truth = {}
    for line in (out / "truth.tsv").read_text().splitlines():
        h, f = line.split("\t")
        truth[canon(h)] = truth.get(canon(h), 0.0) + float(f)
    return rows, truth


# ------------------------------------------------------------------ scoring


def canon(hap: str) -> str:
    """Haplotype key independent of locus order: alleles sorted by locus name."""
    return "~".join(sorted(hap.split("~"), key=lambda a: a.split("*")[0]))


def log_likelihood(rows, freqs: dict[str, float]) -> float:
    """log L = sum over people of log sum over phases of f(h) f(h') (x2 if h != h').

    Computed here from the data and a frequency table, so all tools are scored
    by the same function. Returns -inf if some person has probability zero.
    """
    total = 0.0
    for _, per in rows:
        het = [k for k, (a, b) in enumerate(per) if a != b]
        p = 0.0
        if not het:
            h = canon("~".join(a for a, _ in per))
            p = freqs.get(h, 0.0) ** 2
        else:
            # Fix the first heterozygous locus to remove mirror duplicates.
            for flips in itertools.product((False, True), repeat=len(het) - 1):
                flip = dict(zip(het[1:], flips, strict=True))
                h1 = canon("~".join(b if flip.get(k) else a for k, (a, b) in enumerate(per)))
                h2 = canon("~".join(a if flip.get(k) else b for k, (a, b) in enumerate(per)))
                p += 2 * freqs.get(h1, 0.0) * freqs.get(h2, 0.0)
        if p <= 0:
            return float("-inf")
        total += math.log(p)
    return total


def slug(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def save_freqs(path: Path, freqs: dict[str, float]) -> None:
    """Haplotype<TAB>frequency, largest first (for plotting)."""
    with path.open("w") as f:
        for h, v in sorted(freqs.items(), key=lambda x: -x[1]):
            f.write(f"{h}\t{v:.12g}\n")


def compare(freqs: dict[str, float], ref: dict[str, float]) -> dict:
    keys = set(freqs) | set(ref)
    diffs = [abs(freqs.get(h, 0.0) - ref.get(h, 0.0)) for h in keys]
    return {"max_abs_diff": max(diffs, default=0.0), "tvd": 0.5 * sum(diffs)}


# ------------------------------------------------------------------ runners


def timed(cmd: list[str], cwd: Path, timeout: float) -> dict:
    """Run `cmd`, returning wall time, peak memory (MB) and exit status."""
    probe = (
        "import resource, subprocess, sys, time\n"
        "t = time.perf_counter()\n"
        "r = subprocess.run(sys.argv[1:])\n"
        "print('BENCH', r.returncode, time.perf_counter() - t,"
        " resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024, file=sys.stderr)\n"
    )
    try:
        run = subprocess.run(
            [sys.executable, "-c", probe, *cmd], cwd=cwd, capture_output=True, text=True, timeout=timeout
        )
    except subprocess.TimeoutExpired:
        return {"status": f"timeout ({timeout:.0f}s)"}
    line = next((l for l in run.stderr.splitlines() if l.startswith("BENCH ")), None)
    if line is None:
        return {"status": "failed", "stderr": run.stderr[-2000:]}
    _, code, seconds, mb = line.split()
    out = {"status": "ok" if code == "0" else f"exit {code}", "seconds": float(seconds), "peak_mb": float(mb)}
    out["stdout"], out["stderr"] = run.stdout[-4000:], run.stderr[-4000:]
    return out


def read_hfs(path: Path) -> dict[str, float]:
    out = {}
    for line in path.read_text().splitlines():
        if line.strip():
            h, f = line.split("\t")[:2]
            out[canon(h)] = out.get(canon(h), 0.0) + float(f)
    return out


def run_ours(
    work: Path, data: Path, loci, epsilon: float, threads: int | None, timeout: float, starts: int = 1
) -> dict:
    out = work / f"ours_t{threads or 'all'}_e{epsilon:g}_s{starts}"
    cmd = [
        sys.executable, "-m", "hlahaplotools", "estimate", str(data),
        "--loci", ",".join(f"{l}:asItIs" for l in loci), "-o", str(out),
        "--epsilon", str(epsilon), "--cut", "0", "--no-drb345", "-q",
        "--data-dir", str(ROOT / "data"), "--starts", str(starts),
    ]  # fmt: skip
    if threads:
        cmd += ["--threads", str(threads)]
    r = timed(cmd, work, timeout)
    if r["status"] == "ok":
        r["freqs"] = read_hfs(out / "hfs.dat")
        info = json.loads((out / "run.json").read_text())
        em = info.get("em", info)
        r["iterations"] = em.get("iterations")
        r["reported_logl"] = em.get("log_likelihood")
    return r


def run_haplomat(work: Path, data: Path, loci, epsilon: float, binary: Path, timeout: float) -> dict:
    if not binary.exists():
        return {"status": "not installed"}
    cpp = work / "haplomat"
    (cpp / "results").mkdir(parents=True, exist_ok=True)
    data_dir = binary.parent / "data"
    if not (cpp / "data").exists():
        os.symlink(data_dir, cpp / "data")
    (cpp / "parametersMAC").write_text(
        "\n".join(
            [
                f"FILENAME_INPUT={data}",
                "FILENAME_HAPLOTYPES=results/haplotypes.dat",
                "FILENAME_GENOTYPES=results/genotypes.dat",
                "FILENAME_HAPLOTYPEFREQUENCIES=results/hfs.dat",
                "FILENAME_EPSILON_LOGL=results/epsilon.dat",
                "FILENAME_ANALYTICS=results/analytics",
                "LOCI_AND_RESOLUTIONS=" + ",".join(f"{l}:asItIs" for l in loci),
                "MINIMAL_FREQUENCY_GENOTYPES=1e-5",
                "DO_AMBIGUITYFILTER=false",
                "EXPAND_LINES_AMBIGUITYFILTER=false",
                "WRITE_GENOTYPES=false",
                "DO_ANALYTICS=false",
                "INITIALIZATION_HAPLOTYPEFREQUENCIES=numberOccurrence",
                f"EPSILON={epsilon}",
                "CUT_HAPLOTYPEFREQUENCIES=0",
                "RENORMALIZE_HAPLOTYPEFREQUENCIES=false",
                "SEED=0",
            ]
        )
        + "\n"
    )
    r = timed([str(binary), "MAC"], cpp, timeout)
    if r["status"] == "ok":
        r["freqs"] = read_hfs(cpp / "results" / "hfs.dat")
        m = re.search(r"Used\s+(\d+)\s+steps", r["stdout"])
        r["iterations"] = int(m.group(1)) if m else None
        # Last line of the epsilon/logL trace: "epsilon<TAB>logL".
        trace = (cpp / "results" / "epsilon.dat").read_text().split()
        r["reported_logl"] = float(trace[-1]) if trace else None
    return r


def run_pypop(work: Path, pop: Path, loci, binary: Path, timeout: float, init_cond: int = 50) -> dict:
    if not binary.exists():
        return {"status": "not installed"}
    py = work / "pypop"
    py.mkdir(exist_ok=True)
    fields = "\n".join([" id"] + [f" *{l.lower()}_{i}" for l in loci for i in (1, 2)])
    (py / "config.ini").write_text(
        "[General]\ndebug=0\n\n[ParseGenotypeFile]\nuntypedAllele=****\nalleleDesignator=*\n"
        f"validPopFields=populationname\nvalidSampleFields={fields.strip()}\n\n"
        f"[Emhaplofreq]\nlociToEstHaplo={':'.join(l.upper() for l in loci)}\nnumInitCond={init_cond}\n"
    )
    r = timed([str(binary), "-c", "config.ini", str(pop), "-o", str(py)], py, timeout)
    return parse_pypop(py, loci, r)


def parse_pypop(py: Path, loci, r: dict) -> dict:
    """Read PyPop's XML output. PyPop can finish the EM, write the XML and then fail
    while formatting its text report; the estimate is used and the failure noted."""
    xml = next(py.glob("*-out.xml"), None)
    if xml is None or "<haplotypefreq>" not in xml.read_text():
        if r.get("status") == "ok":
            r["status"] = "no output"
        return r
    if r.get("status") != "ok":
        r["note_exit"] = f"{r.get('status')} after writing its estimate (text report failed)"
        r["status"] = "ok"
    text = xml.read_text()
    block = re.search(r"<haplotypefreq>(.*?)</haplotypefreq>", text, re.S)
    if block is None:
        role = re.search(r'<group [^>]*role="([^"]+)"', text)
        r["status"] = f"no estimate ({role.group(1) if role else 'unknown reason'})"
        return r
    body = block.group(1)
    freqs = {}
    for name, f in re.findall(r'<haplotype name="([^"]+)"><frequency>([\d.eE+-]+)', body):
        hap = "~".join(f"{l}*{a}" for l, a in zip(loci, name.split("~"), strict=True))
        freqs[canon(hap)] = float(f)
    r["freqs"] = freqs
    m = re.search(r"<loglikelihood>([-\d.eE+]+)</loglikelihood>", body)
    r["reported_logl"] = float(m.group(1)) if m else None
    m = re.search(r"<iterConverged>(\d+)</iterConverged>", body)
    r["iterations"] = int(m.group(1)) if m else None
    r["converged"] = '<condition role="converged"/>' in body
    starts = [float(x) for x in re.findall(r"Log likelihood after\s+\d+ iterations:\s*([-\d.]+)", text)]
    if starts:
        r["start_logl"] = {"best": max(starts), "worst": min(starts), "n": len(starts)}
    r["note"] = "frequencies printed to 5 decimals by PyPop"
    return r


# ------------------------------------------------------------------ cases


def run_case(name: str, rows, loci, args, truth: dict[str, float] | None = None, info: dict | None = None) -> dict:
    work = Path(args.out).resolve() / "cases" / name
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    data, pop = work / "input.dat", work / "input.pop"
    write_mac(data, loci, rows)
    write_pop(pop, loci, rows)
    case = {"case": name, "loci": loci, "n": len(rows), "info": info or {}, "tools": {}}
    print(f"\n== {name}: {len(rows)} people, loci {'~'.join(loci)}", flush=True)

    # Reference: the best of many starts, run to a tight tolerance (the maximum-likelihood estimate).
    ref = run_ours(work, data, loci, 1e-10, None, args.timeout, starts=args.ref_starts)
    if ref["status"] != "ok":
        case["error"] = f"reference run failed: {ref['status']}: {ref.get('stderr', '')[-500:]}"
        print("  " + case["error"])
        return case
    ref_freqs = ref.pop("freqs")
    save_freqs(work / "freqs_reference.tsv", ref_freqs)
    if truth:
        save_freqs(work / "freqs_truth.tsv", truth)
    case["reference"] = {
        "logl": log_likelihood(rows, ref_freqs),
        "haplotypes": len(ref_freqs),
        "starts": args.ref_starts,
        "seconds": ref.get("seconds"),
    }

    runs = {}
    for t in args.threads:
        runs[f"hlahaplotools, 1 start, {t} thread(s)"] = run_ours(work, data, loci, args.epsilon, t, args.timeout)
    runs[f"hlahaplotools, {args.starts} starts"] = run_ours(
        work, data, loci, args.epsilon, max(args.threads), args.timeout, starts=args.starts
    )
    runs["Hapl-o-Mat"] = run_haplomat(work, data, loci, args.epsilon, Path(args.haplomat), args.timeout)
    if len(rows) <= 5000 or args.pypop_all:
        runs["PyPop"] = run_pypop(work, pop, loci, Path(args.pypop), args.timeout)
    else:
        runs["PyPop"] = {"status": "skipped (more than 5000 people, PyPop's limit)"}

    for tool, r in runs.items():
        freqs = r.pop("freqs", None)
        if freqs is not None:
            save_freqs(work / f"freqs_{slug(tool)}.tsv", freqs)
        r.pop("stdout", None)
        if r.get("status") == "ok":
            r.pop("stderr", None)
        if freqs is not None:
            total = sum(freqs.values())
            r["haplotypes"] = len(freqs)
            r["vs_reference"] = compare(freqs, ref_freqs)
            r["logl"] = log_likelihood(rows, {h: f / total for h, f in freqs.items()}) if total > 0 else None
            r["delta_logl"] = r["logl"] - case["reference"]["logl"] if r["logl"] is not None else None
            if truth:
                r["vs_truth"] = compare(freqs, truth)
        case["tools"][tool] = r
        line = f"  {tool:28} {r.get('status')}"
        if "vs_reference" in r:
            line += (
                f"  max|diff| {r['vs_reference']['max_abs_diff']:.2e}  TVD {r['vs_reference']['tvd']:.2e}"
                f"  dlogL {r['delta_logl']:.2e}  {r['seconds']:.2f}s  {r['peak_mb']:.0f} MB"
            )
        print(line, flush=True)
    if truth:
        case["reference"]["vs_truth"] = compare(ref_freqs, truth)
    (work / "case.json").write_text(json.dumps(case, indent=2) + "\n")
    return case


def environment() -> dict:
    def ver(cmd):
        try:
            return subprocess.run(cmd, capture_output=True, text=True, timeout=60).stdout.strip().splitlines()[0]
        except Exception:
            return None

    return {
        "platform": platform.platform(),
        "cpu_count": os.cpu_count(),
        "python": sys.version.split()[0],
        "hlahaplotools": ver([sys.executable, "-m", "hlahaplotools", "--version"]),
        "pypop": ver([str(PYPOP), "--version"]),
        "haplomat": str(HAPLOMAT) if HAPLOMAT.exists() else None,
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
    }


def cmd_run(args) -> int:
    import hlahaplotools as h

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cases = []
    sizes = [int(x) for x in args.sizes.split(",") if x]  # --sizes "" skips the simulated cases

    # Simulated populations: known truth, increasing size and locus count.
    for n in sizes:
        for k in (2, 4, 6):
            loci = SIM_LOCI[:k]
            rows, truth = simulate(out / "sim" / f"n{n}_l{k}", loci, n, seed=args.seed)
            cases.append(run_case(f"sim_n{n}_{k}loci", rows, loci, args, truth=truth))

    # Real data, resolved to one name per typing.
    if args.real:
        manifest = json.loads((ROOT / "data" / "manifest.json").read_text())
        nom = h.Nomenclature(str(ROOT / "data" / "derived" / manifest["imgthla"]["release"]))
        sets = [s.split("~") for s in args.real_sets.split(",")] if args.real_sets else REAL_SETS
        for loci in sets:
            rows, info = prepare_real(Path(args.real), loci, args.resolution, nom, args.real_n, args.seed)
            name = f"real_{args.resolution}_{'_'.join(loci)}_n{len(rows)}"
            cases.append(run_case(name, rows, loci, args, info=info))

    result = {"environment": environment(), "settings": vars(args), "cases": cases}
    (out / "benchmark.json").write_text(json.dumps(result, indent=2, default=str) + "\n")
    print(f"\nWrote {out / 'benchmark.json'}; run `report` for the tables.")
    return 0


def cmd_rescore_pypop(args) -> int:
    """Re-read existing PyPop outputs (no rerun) and update the scores in benchmark.json."""
    out = Path(args.out).resolve()
    path = out / "benchmark.json"
    bench = json.loads(path.read_text())
    for case in bench["cases"]:
        work = out / "cases" / case["case"]
        r = case.get("tools", {}).get("PyPop")
        if r is None or not (work / "pypop").exists() or "reference" not in case:
            continue
        loci = case["loci"]
        _, rows = read_mac(work / "input.dat")
        ref = {h: float(f) for h, f in (l.split("\t") for l in (work / "freqs_reference.tsv").read_text().splitlines())}
        truth_file = work / "freqs_truth.tsv"
        truth = (
            {h: float(f) for h, f in (l.split("\t") for l in truth_file.read_text().splitlines())}
            if truth_file.exists()
            else None
        )
        r = parse_pypop(work / "pypop", loci, dict(r))
        freqs = r.pop("freqs", None)
        if freqs is not None:
            save_freqs(work / "freqs_pypop.tsv", freqs)
            total = sum(freqs.values())
            r["haplotypes"] = len(freqs)
            r["vs_reference"] = compare(freqs, ref)
            r["logl"] = log_likelihood(rows, {h: f / total for h, f in freqs.items()})
            r["delta_logl"] = r["logl"] - case["reference"]["logl"]
            if truth:
                r["vs_truth"] = compare(freqs, truth)
        case["tools"]["PyPop"] = r
        print(f"{case['case']}: PyPop {r.get('status')} {r.get('note_exit', '')}")
    path.write_text(json.dumps(bench, indent=2, default=str) + "\n")
    return 0


# ------------------------------------------------------------------ extended loci

# The eight classical loci routinely typed (with DRB3/4/5, the 11 routine loci), then MIC and the
# non-classical loci. Cumulative sets build up in this order; class I and class II are also tested alone.
CLASSICAL8 = ["A", "B", "C", "DRB1", "DQB1", "DPB1", "DQA1", "DPA1"]
EXT_SETS = {
    "classical6": CLASSICAL8[:6],
    "classical8": CLASSICAL8,
    "classical8_mic": CLASSICAL8 + ["MICA", "MICB"],
    "extended15": CLASSICAL8 + ["MICA", "MICB", "E", "F", "G", "H", "J"],
    "class1": ["A", "B", "C"],
    "class1_extended": ["A", "B", "C", "MICA", "MICB", "E", "F", "G", "H", "J"],
    "class2": ["DRB1", "DQA1", "DQB1", "DPA1", "DPB1"],
}


def cmd_extended(args) -> int:
    """Haplotypes over more loci (MICA, MICB, E, F, G, H, J): accuracy, fit, cost.

    1. Fully typed populations for each named locus set (cumulative: 6 and 8
       classical loci, + MICA/MICB, + E/F/G/H/J; and class I and class II alone);
       every tool.
    2. MICA and MICB untyped in some people (as when MIC typing is added to a
       cohort later), on the 8 classical loci + MIC: hlahaplotools on everyone
       against every tool on the fully typed people only, scored against the truth.
    """
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cases = []
    for name in args.ext_sets.split(","):
        loci = EXT_SETS[name]
        rows, truth = simulate(out / "sim_ext" / f"n{args.n}_{name}", loci, args.n, seed=args.seed)
        case = run_case(f"ext_n{args.n}_{name}_{len(loci)}loci", rows, loci, args, truth=truth)
        case["set"] = name
        cases.append(case)

    # Incomplete MIC typing, on the 8 classical loci + MICA and MICB.
    loci = EXT_SETS["classical8_mic"]
    k = len(loci)
    sim = out / "sim_ext" / f"n{args.n}_l{k}_mic_missing"
    sim.mkdir(parents=True, exist_ok=True)
    rate = args.mic_missing
    subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "simulate.py"), str(sim), "--samples", str(args.n),
         "--loci", ",".join(loci), "--haplotypes", "300", "--alleles", "25", "--seed", str(args.seed),
         "--missing", f"MICA:{rate},MICB:{rate}"],
        check=True, capture_output=True,
    )  # fmt: skip
    truth = {}
    for line in (sim / "truth.tsv").read_text().splitlines():
        h, f = line.split("\t")
        truth[canon(h)] = truth.get(canon(h), 0.0) + float(f)
    _, all_rows = read_mac(sim / "population.dat")
    complete = [(sid, per) for sid, per in all_rows if all(a and b for a, b in per)]
    print(
        f"\n== MICA/MICB untyped in {rate:.0%}: {len(all_rows)} people, {len(complete)} fully typed, "
        f"loci {'~'.join(loci)}"
    )
    # Every tool on the fully typed people (what Hapl-o-Mat and PyPop can use).
    case = run_case(f"ext_mic_missing_complete_only_n{len(complete)}_{k}loci", complete, loci, args, truth=truth)
    # hlahaplotools on everyone, MIC left untyped where it was not typed.
    work = Path(args.out).resolve() / "cases" / f"ext_mic_missing_all_n{len(all_rows)}_{k}loci"
    work.mkdir(parents=True, exist_ok=True)
    data = sim / "population.dat"
    r = run_ours(work, data.resolve(), loci, args.epsilon, max(args.threads), args.timeout, starts=args.starts)
    freqs = r.pop("freqs", None)
    r.pop("stdout", None)
    if freqs is not None:
        r["vs_truth"] = compare(freqs, truth)
        r["haplotypes"] = len(freqs)
        save_freqs(work / "freqs_everyone.tsv", freqs)
        save_freqs(work / "freqs_truth.tsv", truth)
    case["everyone"] = {"n": len(all_rows), "tool": f"hlahaplotools, {args.starts} starts, everyone", **r}
    print(
        f"  everyone ({len(all_rows)} people), hlahaplotools: {r.get('status')}  "
        f"TVD to truth {r.get('vs_truth', {}).get('tvd', float('nan')):.4f}  {r.get('seconds', 0):.1f}s  "
        f"{r.get('peak_mb', 0):.0f} MB"
    )
    cases.append(case)
    (out / "extended.json").write_text(
        json.dumps({"environment": environment(), "settings": vars(args), "cases": cases}, indent=2, default=str) + "\n"
    )
    print(f"\nWrote {out / 'extended.json'}")
    return 0


# ------------------------------------------------------------------ multi-start and bootstrap


def _estimate_rows(h, nom, rows, loci, tmp: Path, **kw):
    write_mac(tmp, loci, rows)
    ds = h.read_mac(str(tmp), half_typed="unknown", drb345=False)
    return h.estimate(ds, nom, "asItIs", **kw)


def cmd_multistart(args) -> int:
    """Different start values: spread of log-likelihoods and frequencies."""
    import hlahaplotools as h

    manifest = json.loads((ROOT / "data" / "manifest.json").read_text())
    nom = h.Nomenclature(str(ROOT / "data" / "derived" / manifest["imgthla"]["release"]))
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    results = []
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp) / "x.dat"
        for loci in args.sets:
            rows, info = prepare_real(Path(args.real), loci, args.resolution, nom, args.real_n, args.seed)
            runs = []
            for init, seed in [("numberOccurrence", 0)] + [
                (init, s) for init in ("random", "perturbation") for s in range(1, args.starts + 1)
            ]:
                est = _estimate_rows(h, nom, rows, loci, tmp, init=init, seed=seed, epsilon=1e-10)
                runs.append(
                    {
                        "init": init,
                        "seed": seed,
                        "logl": est.log_likelihood,
                        "iterations": est.iterations,
                        "freqs": {canon(k): v for k, v in est.frequencies(cut=0.0, renormalize=False)},
                    }
                )
            best = max(runs, key=lambda r: r["logl"])
            for r in runs:
                r.update(compare(r.pop("freqs"), best["freqs"]) if r is not best else {"max_abs_diff": 0.0, "tvd": 0.0})
            best.pop("freqs", None)
            entry = {
                "loci": loci,
                "n": len(rows),
                "best_logl": best["logl"],
                "logl_range": best["logl"] - min(r["logl"] for r in runs),
                "max_tvd_to_best": max(r["tvd"] for r in runs),
                "runs": runs,
            }
            results.append(entry)
            print(
                f"{'~'.join(loci)} n={len(rows)}: {len(runs)} starts, logL range {entry['logl_range']:.2e}, "
                f"largest TVD to best {entry['max_tvd_to_best']:.2e}",
                flush=True,
            )
    (out / "multistart.json").write_text(json.dumps(results, indent=2) + "\n")
    return 0


def cmd_bootstrap(args) -> int:
    """Resample people with replacement: standard errors and 95% intervals."""
    import hlahaplotools as h

    manifest = json.loads((ROOT / "data" / "manifest.json").read_text())
    nom = h.Nomenclature(str(ROOT / "data" / "derived" / manifest["imgthla"]["release"]))
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    results = []
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp) / "x.dat"
        for loci in args.sets:
            rows, info = prepare_real(Path(args.real), loci, args.resolution, nom, None, args.seed)
            full = _estimate_rows(h, nom, rows, loci, tmp, epsilon=1e-8)
            point = {canon(k): v for k, v in full.frequencies(cut=0.0, renormalize=False)}
            top = sorted(point, key=lambda k: -point[k])[: args.top]
            rng = random.Random(args.seed)
            reps = {k: [] for k in top}
            t = time.perf_counter()
            for b in range(args.replicates):
                sample = [rows[rng.randrange(len(rows))] for _ in rows]
                sample = [(f"b{i}", per) for i, (_, per) in enumerate(sample)]
                est = _estimate_rows(h, nom, sample, loci, tmp, epsilon=1e-8)
                f = {canon(k): v for k, v in est.frequencies(cut=0.0, renormalize=False)}
                for k in top:
                    reps[k].append(f.get(k, 0.0))
                if (b + 1) % 50 == 0:
                    print(f"  {'~'.join(loci)}: {b + 1}/{args.replicates} replicates", flush=True)
            table = []
            for k in top:
                v = sorted(reps[k])
                mean = sum(v) / len(v)
                se = math.sqrt(sum((x - mean) ** 2 for x in v) / (len(v) - 1))
                lo, hi = v[int(0.025 * (len(v) - 1))], v[int(math.ceil(0.975 * (len(v) - 1)))]
                table.append({"haplotype": k, "frequency": point[k], "se": se, "ci95": [lo, hi]})
            entry = {
                "loci": loci,
                "n": len(rows),
                "replicates": args.replicates,
                "seconds": time.perf_counter() - t,
                "haplotypes": table,
            }
            results.append(entry)
            print(
                f"{'~'.join(loci)} n={len(rows)}: {args.replicates} replicates in {entry['seconds']:.0f}s", flush=True
            )
    (out / "bootstrap.json").write_text(json.dumps(results, indent=2) + "\n")
    return 0


# ------------------------------------------------------------------ report


def fmt(x, spec=".2e"):
    return "-" if x is None else format(x, spec)


def cmd_report(args) -> int:
    out = Path(args.out)
    lines = ["# Benchmark: hlahaplotools, Hapl-o-Mat and PyPop", ""]
    bench = json.loads((out / "benchmark.json").read_text()) if (out / "benchmark.json").exists() else None
    if bench:
        env = bench["environment"]
        lines += [
            f"Run {env['date']} on {env['platform']}, {env['cpu_count']} CPUs; {env['hlahaplotools']}; "
            f"{env['pypop']}; Hapl-o-Mat from {env['haplomat']}.",
            "",
            "Every tool gets the same file (fully typed, one name per typing, used as typed). The reference is "
            "hlahaplotools, best of `ref_starts` starts at epsilon 1e-10. ΔlogL is the log-likelihood of the tool's frequencies minus the "  # noqa: E501
            "reference's, computed by the benchmark itself; negative means a worse fit. PyPop prints "
            "frequencies to 5 decimals, which limits its agreement to about 1e-5.",
            "",
            "| Case | People | Tool | Status | max abs diff | TVD | ΔlogL | TVD to truth | Iterations | Time (s) | Peak MB |",  # noqa: E501
            "| --- | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
        tsv = ["case\tn\ttool\tstatus\tmax_abs_diff\ttvd\tdelta_logl\ttvd_truth\titerations\tseconds\tpeak_mb"]
        for c in bench["cases"]:
            for tool, r in c.get("tools", {}).items():
                vr, vt = r.get("vs_reference", {}), r.get("vs_truth", {})
                cells = [
                    c["case"], str(c["n"]), tool, r.get("status", "?"),
                    fmt(vr.get("max_abs_diff")), fmt(vr.get("tvd")), fmt(r.get("delta_logl"), ".3g"),
                    fmt(vt.get("tvd"), ".4f"), str(r.get("iterations") or "-"),
                    fmt(r.get("seconds"), ".2f"), fmt(r.get("peak_mb"), ".0f"),
                ]  # fmt: skip
                lines.append("| " + " | ".join(cells) + " |")
                tsv.append("\t".join(cells))
            if "reference" in c and "vs_truth" in c["reference"]:
                lines.append(
                    f"| {c['case']} | {c['n']} | reference (1e-10) | ok | 0 | 0 | 0 | "
                    f"{c['reference']['vs_truth']['tvd']:.4f} | - | - | - |"
                )
        (out / "benchmark.tsv").write_text("\n".join(tsv) + "\n")
    ms = out / "multistart.json"
    if ms.exists():
        lines += ["", "## Start values (hlahaplotools)", "",
                  "| Loci | People | Starts | logL range | Largest TVD to best |", "| --- | ---: | ---: | ---: | ---: |"]  # fmt: skip  # noqa: E501
        for e in json.loads(ms.read_text()):
            lines.append(
                f"| {'~'.join(e['loci'])} | {e['n']} | {len(e['runs'])} | {e['logl_range']:.2e} | {e['max_tvd_to_best']:.2e} |"  # noqa: E501
            )
    bs = out / "bootstrap.json"
    if bs.exists():
        for e in json.loads(bs.read_text()):
            lines += ["", f"## Bootstrap: {'~'.join(e['loci'])} ({e['n']} people, {e['replicates']} replicates)", "",
                      "| Haplotype | Frequency | SE | 95% interval |", "| --- | ---: | ---: | --- |"]  # fmt: skip
            for t in e["haplotypes"]:
                lines.append(
                    f"| `{t['haplotype']}` | {t['frequency']:.4f} | {t['se']:.4f} | {t['ci95'][0]:.4f}–{t['ci95'][1]:.4f} |"  # noqa: E501
                )
    (out / "report.md").write_text("\n".join(lines) + "\n")
    print(f"Wrote {out / 'report.md'}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p):
        p.add_argument("--out", default="bench/results")
        p.add_argument("--real", help="MAC typing file for the real-data cases")
        p.add_argument("--resolution", default="2f", help="resolution the real data is reduced to (default 2f)")
        p.add_argument("--seed", type=int, default=1)

    r = sub.add_parser("run", help="three-way comparison")
    common(r)
    r.add_argument("--sizes", default="1000,5000", help="simulated population sizes")
    r.add_argument("--real-n", type=int, default=5000, help="people sampled from the real data (PyPop allows 5000)")
    r.add_argument("--real-sets", help="real-data locus sets instead of the default list, e.g. A~B~C,DRB1~DQB1")
    r.add_argument("--threads", default="1,8", help="hlahaplotools thread counts to time")
    r.add_argument("--epsilon", type=float, default=1e-6, help="stopping rule for hlahaplotools and Hapl-o-Mat")
    r.add_argument("--starts", type=int, default=10, help="starts for the multi-start hlahaplotools row")
    r.add_argument("--ref-starts", type=int, default=40, help="starts for the reference estimate")
    r.add_argument("--timeout", type=float, default=3600)
    r.add_argument("--haplomat", default=str(HAPLOMAT))
    r.add_argument("--pypop", default=str(PYPOP))
    r.add_argument("--pypop-all", action="store_true", help="also try PyPop above 5000 people")
    r.set_defaults(func=cmd_run)

    x = sub.add_parser("extended", help="more loci (MICA, MICB, E, F, G, H, J) and incomplete MIC typing")
    x.add_argument("--out", default="bench/results")
    x.add_argument("--seed", type=int, default=1)
    x.add_argument("--n", type=int, default=2000, help="simulated people")
    x.add_argument(
        "--sets",
        dest="ext_sets",
        default=",".join(EXT_SETS),
        help=f"named locus sets, comma-separated (default: all of {', '.join(EXT_SETS)})",
    )
    x.add_argument("--mic-missing", type=float, default=0.3, help="share of people with MICA and MICB untyped")
    x.add_argument("--threads", default="8")
    x.add_argument("--starts", type=int, default=10)
    x.add_argument("--ref-starts", type=int, default=20)
    x.add_argument("--epsilon", type=float, default=1e-6)
    x.add_argument("--timeout", type=float, default=3600)
    x.add_argument("--haplomat", default=str(HAPLOMAT))
    x.add_argument("--pypop", default=str(PYPOP))
    x.add_argument("--pypop-all", action="store_true")
    x.set_defaults(func=cmd_extended)

    def sets(p):
        p.add_argument("--sets", default="A,B,C;A,B,C,DRB1;A,B,C,DRB1,DQB1,DPB1", help="locus sets, ';'-separated")

    m = sub.add_parser("multistart", help="spread over start values (hlahaplotools)")
    common(m)
    sets(m)
    m.add_argument("--real-n", type=int, default=5000)
    m.add_argument("--starts", type=int, default=10, help="seeds per random start method")
    m.set_defaults(func=cmd_multistart)

    b = sub.add_parser("bootstrap", help="bootstrap standard errors (hlahaplotools)")
    common(b)
    sets(b)
    b.add_argument("--replicates", type=int, default=200)
    b.add_argument("--top", type=int, default=10)
    b.set_defaults(func=cmd_bootstrap)

    rs = sub.add_parser("rescore-pypop", help="re-read PyPop outputs and update benchmark.json")
    rs.add_argument("--out", default="bench/results")
    rs.set_defaults(func=cmd_rescore_pypop)

    rep = sub.add_parser("report", help="write report.md and benchmark.tsv from the results")
    rep.add_argument("--out", default="bench/results")
    rep.set_defaults(func=cmd_report)

    args = ap.parse_args()
    if getattr(args, "threads", None) is not None:
        args.threads = [int(t) for t in args.threads.split(",")]
    if getattr(args, "sets", None) is not None and isinstance(args.sets, str):
        args.sets = [s.split(",") for s in args.sets.split(";")]
    if args.cmd in ("multistart", "bootstrap") and not args.real:
        sys.exit("error: --real is needed for multistart and bootstrap")
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
