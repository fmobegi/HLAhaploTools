"""Run the original C++ Hapl-o-Mat and hlahaplotools on the same input and
compare the haplotype frequencies.

Usage (from the project root, after building reference/Hapl-o-Mat with
`make` and its data with prepareData, or copying data/derived/<release> into
reference/Hapl-o-Mat/data):

    uv run python scripts/compare_haplomat.py MAC reference/Hapl-o-Mat/examplePopulations/populationData_a.dat \
        --loci A:g,B:g,DRB1:1f

Both tools get the same parameter file. Hapl-o-Mat's random initialisation
depends on hash-table order, so comparisons use numberOccurrence by default.
Exit status is 1 if any frequency differs by more than --tolerance.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
HAPLOMAT = ROOT / "reference" / "Hapl-o-Mat"


def gls_loci_order(pull: Path, glid: Path) -> str:
    """Hapl-o-Mat's LOCIORDER: the locus of each pull-file column, read from the GLID file."""
    locus_of = {}
    for line in glid.read_text().splitlines():
        if ";" in line:
            n, gl = line.split(";", 1)
            locus_of[n.strip()] = gl.strip().replace("HLA-", "").split("*")[0]
    order: list[str | None] = []
    for line in pull.read_text().splitlines():
        if ";" not in line:
            continue
        ids = line.split(";", 1)[1].split(":")
        order += [None] * (len(ids) - len(order))
        for k, n in enumerate(ids):
            if order[k] is None and n.strip() not in ("", "0"):
                order[k] = locus_of.get(n.strip())
        if all(order):
            break
    return ",".join(o or "?" for o in order)


def write_parameters(
    path: Path, fmt: str, input_file: Path, loci: str, epsilon: float, init: str, glid: Path | None = None
) -> None:
    if fmt == "GLS":
        lines = [
            f"FILENAME_PULL={input_file}",
            f"FILENAME_GLID={glid}",
            "FILENAME_HAPLOTYPES=results/haplotypes.dat",
            "FILENAME_GENOTYPES=results/genotypes.dat",
            "FILENAME_HAPLOTYPEFREQUENCIES=results/hfs.dat",
            "FILENAME_EPSILON_LOGL=results/epsilon.dat",
            f"LOCIORDER={gls_loci_order(input_file, glid)}",
            f"LOCI_AND_RESOLUTIONS={loci}",
            "MINIMAL_FREQUENCY_GENOTYPES=1e-5",
            "DO_AMBIGUITYFILTER=false",
            "EXPAND_LINES_AMBIGUITYFILTER=false",
            "RESOLVE_MISSING_GENOTYPES=false",
            "WRITE_GENOTYPES=false",
            f"INITIALIZATION_HAPLOTYPEFREQUENCIES={init}",
            f"EPSILON={epsilon}",
            "CUT_HAPLOTYPEFREQUENCIES=1e-6",
            "RENORMALIZE_HAPLOTYPEFREQUENCIES=true",
            "SEED=0",
        ]
        path.write_text("\n".join(lines) + "\n")
        return
    lines = [
        f"FILENAME_INPUT={input_file}",
        "FILENAME_HAPLOTYPES=results/haplotypes.dat",
        "FILENAME_GENOTYPES=results/genotypes.dat",
        "FILENAME_HAPLOTYPEFREQUENCIES=results/hfs.dat",
        "FILENAME_EPSILON_LOGL=results/epsilon.dat",
    ]
    if fmt == "MAC":
        lines.append("FILENAME_ANALYTICS=results/analytics")
    lines += [
        f"LOCI_AND_RESOLUTIONS={loci}",
        "MINIMAL_FREQUENCY_GENOTYPES=1e-5",
        "DO_AMBIGUITYFILTER=false",
        "EXPAND_LINES_AMBIGUITYFILTER=false",
        "WRITE_GENOTYPES=false",
    ]
    if fmt == "MAC":
        lines.append("DO_ANALYTICS=false")
    lines += [
        f"INITIALIZATION_HAPLOTYPEFREQUENCIES={init}",
        f"EPSILON={epsilon}",
        "CUT_HAPLOTYPEFREQUENCIES=1e-6",
        "RENORMALIZE_HAPLOTYPEFREQUENCIES=true",
        "SEED=0",
    ]
    path.write_text("\n".join(lines) + "\n")


def read_hfs(path: Path) -> dict[str, float]:
    out = {}
    for line in path.read_text().splitlines():
        if line.strip():
            h, f = line.split("\t")[:2]
            out[h] = float(f)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("format", choices=["MAC", "GLSC", "GLS"])
    ap.add_argument("input", type=Path, help="MAC/GLSC file, or the GLS pull file")
    ap.add_argument("--glid", type=Path, help="GLS: GLID file (default: the pull file with a .glid suffix)")
    ap.add_argument("--loci", required=True, help="LOCI_AND_RESOLUTIONS, e.g. A:g,B:g,DRB1:1f")
    ap.add_argument("--epsilon", type=float, default=1e-6)
    ap.add_argument("--init", default="numberOccurrence")
    ap.add_argument("--tolerance", type=float, default=1e-6)
    ap.add_argument("--keep", action="store_true", help="keep the working folder")
    args = ap.parse_args()

    binary = HAPLOMAT / "haplomat"
    if not binary.exists():
        sys.exit(f"error: {binary} not found; run `make` in {HAPLOMAT}")
    work = Path(tempfile.mkdtemp(prefix="hh_compare_"))
    input_file = args.input.resolve()
    glid = (args.glid or args.input.with_suffix(".glid")).resolve() if args.format == "GLS" else None

    # C++: needs ./data, ./parametersMAC|GLSC and ./results in its working folder.
    cpp = work / "cpp"
    (cpp / "results").mkdir(parents=True)
    os.symlink(HAPLOMAT / "data", cpp / "data")
    write_parameters(
        cpp / f"parameters{args.format}", args.format, input_file, args.loci, args.epsilon, args.init, glid
    )
    t = time.perf_counter()
    run = subprocess.run([str(binary), args.format], cwd=cpp, capture_output=True, text=True)
    cpp_seconds = time.perf_counter() - t
    if run.returncode != 0:
        print(run.stdout, run.stderr)
        sys.exit("error: Hapl-o-Mat failed")
    steps = next((l.split()[1] for l in run.stdout.splitlines() if "Used" in l and "steps" in l), "?")

    # Ours: same parameter file, same data tables.
    ours = work / "ours"
    (ours / "results").mkdir(parents=True)
    write_parameters(ours / "params", args.format, input_file, args.loci, args.epsilon, args.init, glid)
    t = time.perf_counter()
    run2 = subprocess.run(
        [
            sys.executable,
            "-m",
            "hlahaplotools",
            "haplomat",
            args.format,
            "params",
            "--nomenclature",
            str(HAPLOMAT / "data"),
            "--mac-codes",
            str(ROOT / "data" / "mac" / "alpha.v3.txt"),
            "-q",
        ],
        cwd=ours,
        capture_output=True,
        text=True,
    )
    our_seconds = time.perf_counter() - t
    if run2.returncode not in (0, 2):
        print(run2.stdout, run2.stderr)
        sys.exit("error: hlahaplotools failed")

    a = read_hfs(cpp / "results" / "hfs.dat")
    b = read_hfs(ours / "results" / "hfs.dat")
    only_cpp = sorted(set(a) - set(b), key=lambda h: -a[h])
    only_ours = sorted(set(b) - set(a), key=lambda h: -b[h])
    diffs = sorted(((abs(a[h] - b[h]), h) for h in set(a) & set(b)), reverse=True)
    worst = max([d for d, _ in diffs] + [a[h] for h in only_cpp] + [b[h] for h in only_ours] + [0.0])

    print(f"input: {args.input}  loci: {args.loci}  init: {args.init}  epsilon: {args.epsilon}")
    print(f"Hapl-o-Mat:    {len(a):6} haplotypes, {steps} steps, {cpp_seconds:.2f}s")
    print(f"hlahaplotools: {len(b):6} haplotypes, {our_seconds:.2f}s (includes Python start-up)")
    print(f"largest difference: {worst:.3e}  (only in Hapl-o-Mat: {len(only_cpp)}, only in ours: {len(only_ours)})")
    for _, h in diffs[:3]:
        print(f"  {h}: {a[h]:.10f} vs {b[h]:.10f}")
    for h in only_cpp[:3]:
        print(f"  only Hapl-o-Mat: {h} {a[h]:.3e}")
    for h in only_ours[:3]:
        print(f"  only ours: {h} {b[h]:.3e}")
    if args.keep:
        print(f"kept {work}")
    else:
        shutil.rmtree(work)
    return 0 if worst <= args.tolerance else 1


if __name__ == "__main__":
    sys.exit(main())
