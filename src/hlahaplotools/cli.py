"""Command-line interface: `hlahaplotools {estimate,families,haplomat,summary}`."""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from hlahaplotools._core import Nomenclature, estimate, read_glsc, read_mac

INITS = ["numberOccurrence", "perturbation", "equal", "random"]


def _version() -> str:
    try:
        return version("hlahaplotools")
    except PackageNotFoundError:
        return "unknown"


# ---------------------------------------------------------------- data lookup


def find_data(data_dir: str | None) -> tuple[Path | None, dict]:
    """Locate the data folder written by scripts/build_nomenclature.py.

    Order: --data-dir, $HLAHAPLOTOOLS_DATA, ./data. Returns (folder, manifest).
    """
    for candidate in (data_dir, os.environ.get("HLAHAPLOTOOLS_DATA"), "data"):
        if candidate and (Path(candidate) / "manifest.json").exists():
            folder = Path(candidate)
            return folder, json.loads((folder / "manifest.json").read_text())
    return None, {}


def resolve_data_paths(args) -> tuple[Path, Path | None, dict]:
    """Nomenclature folder and MAC table, from explicit options or the manifest."""
    folder, manifest = find_data(args.data_dir)
    nom = Path(args.nomenclature) if args.nomenclature else None
    mac = Path(args.mac_codes) if args.mac_codes else None
    if nom is None and folder is not None and "derived" in manifest:
        nom = folder / "derived" / manifest["imgthla"]["release"]
    if mac is None and folder is not None and "mac" in manifest:
        mac = folder / "mac" / "alpha.v3.txt"
    if nom is None:
        sys.exit(
            "error: no nomenclature tables found. Run `uv run python scripts/build_nomenclature.py`, "
            "or pass --nomenclature DIR (a data/derived/<release> folder or Hapl-o-Mat's data/)."
        )
    return nom, mac, manifest


# ---------------------------------------------------------------- running


def read_input(path: str, fmt: str | None, mac: Path | None, loci: list[str] | None):
    fmt = fmt or ("glsc" if path.lower().endswith((".glc", ".glsc")) else "mac")
    reader = read_glsc if fmt == "glsc" else read_mac
    return reader(path, mac_codes=str(mac) if mac and mac.exists() else None, loci=loci), fmt


def run_estimate(
    *,
    input_path: str,
    fmt: str | None,
    loci_res: list[tuple[str, str]],
    nom_dir: Path,
    mac: Path | None,
    manifest: dict,
    out_freqs: Path,
    out_history: Path | None,
    epsilon: float,
    init: str,
    seed: int,
    cut: float,
    renormalize: bool,
    min_genotype_frequency: float,
    max_completions: int,
    threads: int | None,
    quiet: bool = False,
) -> int:
    say = (lambda *a: None) if quiet else print
    loci = [l for l, _ in loci_res]
    t0 = time.perf_counter()
    ds, fmt = read_input(input_path, fmt, mac, loci)
    nom = Nomenclature(str(nom_dir))
    t1 = time.perf_counter()
    say(f"Read {len(ds)} samples from {input_path} ({fmt.upper()}); {len(ds.issues)} lines unreadable")
    missing, half = ds.missing_counts(), ds.one_allele_counts()
    for locus, res in loci_res:
        say(f"  {locus:6} -> {res:6} untyped: {missing[locus]:6}  one allele typed: {half[locus]:6}")

    est = estimate(
        ds,
        nom,
        [r for _, r in loci_res],
        epsilon=epsilon,
        init=init,
        seed=seed,
        min_genotype_frequency=min_genotype_frequency,
        max_completions=max_completions,
        threads=threads,
    )
    t2 = time.perf_counter()
    say(
        f"EM: {est.n_samples} samples ({est.n_samples_incomplete} incomplete), {est.n_phenotypes} genotypes, "
        f"{est.n_candidates} candidate haplotypes ({est.n_completions} filled in), "
        f"{est.iterations} iterations, converged={est.converged}, logL={est.log_likelihood:.6f}"
    )
    if est.issues:
        say(f"  {len(est.issues)} samples left out of the estimate:")
        for _, sid, reason in est.issues[:10]:
            say(f"    {sid}: {reason}")
        if len(est.issues) > 10:
            say(f"    ... see {out_freqs.with_name('run.json')}")

    out_freqs.parent.mkdir(parents=True, exist_ok=True)
    est.write_frequencies(str(out_freqs), cut=cut, renormalize=renormalize)
    if out_history:
        out_history.parent.mkdir(parents=True, exist_ok=True)
        est.write_history(str(out_history))
    kept = est.frequencies(cut=cut, renormalize=renormalize)

    run_info = {
        "tool": "hlahaplotools",
        "version": _version(),
        "input": str(input_path),
        "format": fmt,
        "loci_and_resolutions": dict(loci_res),
        "imgthla_release": manifest.get("imgthla", {}).get("release") or nom.source,
        "nomenclature_dir": str(nom_dir),
        "mac_table": str(mac) if mac else None,
        "mac_last_updated": manifest.get("mac", {}).get("last_updated"),
        "parameters": {
            "epsilon": epsilon,
            "init": init,
            "seed": seed,
            "cut": cut,
            "renormalize": renormalize,
            "min_genotype_frequency": min_genotype_frequency,
            "max_completions": max_completions,
        },
        "samples_read": len(ds),
        "samples_used": est.n_samples,
        "samples_incomplete": est.n_samples_incomplete,
        "genotypes": est.n_phenotypes,
        "candidate_haplotypes": est.n_candidates,
        "filled_in_haplotypes": est.n_completions,
        "haplotypes_written": len(kept),
        "iterations": est.iterations,
        "converged": est.converged,
        "log_likelihood": est.log_likelihood,
        "seconds": {
            "read": round(t1 - t0, 3),
            "estimate": round(t2 - t1, 3),
            **{k: round(v, 3) for k, v in est.timings.items()},
        },
        "unreadable_lines": [{"line": l, "id": i, "reason": r} for l, i, r in ds.issues],
        "samples_left_out": [{"id": i, "reason": r} for _, i, r in est.issues],
    }
    run_path = out_freqs.with_name("run.json")
    run_path.write_text(json.dumps(run_info, indent=2) + "\n")
    say(f"Wrote {len(kept)} haplotypes to {out_freqs} (details in {run_path})")
    return 0 if est.converged else 2


# ---------------------------------------------------------------- subcommands


def parse_loci(spec: str) -> list[tuple[str, str]]:
    """'A:g,B:g,DRB1:1f' -> [('A','g'), ('B','g'), ('DRB1','1f')]."""
    out = []
    for item in spec.split(","):
        locus, _, res = item.strip().partition(":")
        if not locus or not res:
            sys.exit(f"error: bad locus spec {item!r}; expected LOCUS:RESOLUTION, e.g. A:g")
        out.append((locus, res))
    return out


def cmd_estimate(args) -> int:
    nom, mac, manifest = resolve_data_paths(args)
    out_dir = Path(args.output)
    return run_estimate(
        input_path=args.input,
        fmt=args.format,
        loci_res=parse_loci(args.loci),
        nom_dir=nom,
        mac=mac,
        manifest=manifest,
        out_freqs=out_dir / "hfs.dat",
        out_history=out_dir / "epsilon.dat",
        epsilon=args.epsilon,
        init=args.init,
        seed=args.seed,
        cut=args.cut,
        renormalize=not args.no_renormalize,
        min_genotype_frequency=args.min_genotype_frequency,
        max_completions=args.max_completions,
        threads=args.threads,
        quiet=args.quiet,
    )


def read_parameters(path: str) -> dict[str, str]:
    params = {}
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#") and "=" in line:
            key, _, value = line.partition("=")
            params[key.strip()] = value.strip()
    return params


def cmd_haplomat(args) -> int:
    """Run from a Hapl-o-Mat parameter file (MAC or GLSC)."""
    p = read_parameters(args.parameters)
    for key, what in [
        ("DO_AMBIGUITYFILTER", "the ambiguity filter"),
        ("WRITE_GENOTYPES", "writing genotypes"),
        ("DO_ANALYTICS", "analytics"),
    ]:
        if p.get(key, "false").lower() == "true":
            print(f"warning: {key}=true ignored; {what} is not supported", file=sys.stderr)
    nom, mac, manifest = resolve_data_paths(args)
    # Hapl-o-Mat orders loci alphabetically (it keeps them in a std::map).
    loci_res = sorted(parse_loci(p["LOCI_AND_RESOLUTIONS"]))
    return run_estimate(
        input_path=p["FILENAME_INPUT"],
        fmt=args.format.lower(),
        loci_res=loci_res,
        nom_dir=nom,
        mac=mac,
        manifest=manifest,
        out_freqs=Path(p.get("FILENAME_HAPLOTYPEFREQUENCIES", "results/hfs.dat")),
        out_history=Path(p["FILENAME_EPSILON_LOGL"]) if "FILENAME_EPSILON_LOGL" in p else None,
        epsilon=float(p.get("EPSILON", "1e-6")),
        init=p.get("INITIALIZATION_HAPLOTYPEFREQUENCIES", "perturbation"),
        seed=int(p.get("SEED", "0")),
        cut=float(p.get("CUT_HAPLOTYPEFREQUENCIES", "1e-6")),
        renormalize=p.get("RENORMALIZE_HAPLOTYPEFREQUENCIES", "true").lower() == "true",
        min_genotype_frequency=float(p.get("MINIMAL_FREQUENCY_GENOTYPES", "1e-5")),
        max_completions=args.max_completions,
        threads=args.threads,
        quiet=args.quiet,
    )


def cmd_summary(args) -> int:
    if args.mac_codes:
        mac = Path(args.mac_codes)
    else:
        folder, manifest = find_data(args.data_dir)
        mac = folder / "mac" / "alpha.v3.txt" if folder and "mac" in manifest else None
    loci = args.loci.split(",") if args.loci else None
    ds, fmt = read_input(args.input, args.format, mac, loci)
    print(f"{ds!r} ({fmt.upper()})")
    missing, half = ds.missing_counts(), ds.one_allele_counts()
    for locus in ds.loci:
        print(f"  {locus:6} untyped: {missing[locus]:6}   one allele typed: {half[locus]:6}")
    for line, sid, reason in ds.issues[:20]:
        print(f"  line {line} (id {sid}): {reason}")
    if len(ds.issues) > 20:
        print(f"  ... and {len(ds.issues) - 20} more")
    return 0


# ---------------------------------------------------------------- argparse


def _tsv(path: Path, header: list[str], rows) -> None:
    with path.open("w") as f:
        f.write("\t".join(header) + "\n")
        for row in rows:
            f.write("\t".join("" if v is None else str(v) for v in row) + "\n")


SERIOUS = {"parent_excluded", "possible_duplicate", "not_full_sibling", "trio_inconsistent", "inheritance_mismatch"}


def cmd_families(args) -> int:
    """Family data workflow. Writes to the output folder:

    flags.tsv        pedigree problems: excluded parents, duplicated rows,
                     non-siblings, recombinants, single-locus mismatches
    mendel.tsv       children whose typing cannot come from their parents
    phases.tsv       each family's four parental haplotypes (LOCUS*? = unknown)
    children.tsv     which haplotype each child inherited from each parent
    hfs.dat          haplotype frequencies from founders (parents) only
    epsilon.dat      EM convergence trace
    diplotypes.tsv   most likely haplotype pair per founder
    comparison.json  EM on unphased founders vs. family phasing
    run.json         release, parameters and counts
    """
    from hlahaplotools._core import (
        audit_families,
        check_mendelian,
        compare_to_families,
        estimate_families,
        read_family_table,
    )

    say = (lambda *a: None) if args.quiet else print
    nom_dir, mac, manifest = resolve_data_paths(args)
    loci_res = parse_loci(args.loci)
    loci, res = [l for l, _ in loci_res], {l: r for l, r in loci_res}
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    t0 = time.perf_counter()
    data = read_family_table(args.input, mac_codes=str(mac) if mac and mac.exists() else None, loci=loci)
    nom = Nomenclature(str(nom_dir))
    fams = data.families()
    say(
        f"Read {len(data)} people in {len(fams)} families ({len(data.unrelated)} unrelated); "
        f"{len(data.issues)} rows or members left out"
    )
    for line, sid, reason in data.issues[:10]:
        say(f"  {sid or f'line {line}'}: {reason}")

    mendel = check_mendelian(data, nom, res)
    _tsv(
        out / "mendel.tsv",
        ["family", "sample", "locus", "message"],
        ([m["family"], m["sample"], m["locus"], m["message"]] for m in mendel),
    )
    say(f"Mendel check: {len(mendel)} problems (mendel.tsv)")

    flags = audit_families(data, nom, res, cross_family_duplicates=args.cross_family_duplicates)
    _tsv(
        out / "flags.tsv",
        ["family", "sample", "relative", "kind", "loci", "message"],
        ([f["family"], f["sample"], f["relative"], f["kind"], ",".join(f["loci"]), f["message"]] for f in flags),
    )
    by_kind: dict[str, int] = {}
    for f in flags:
        by_kind[f["kind"]] = by_kind.get(f["kind"], 0) + 1
    say(f"Pedigree checks: {len(flags)} flags (flags.tsv)" + (f": {by_kind}" if by_kind else ""))
    for f in [f for f in flags if f["kind"] in SERIOUS][:10]:
        say(f"  {f['kind']}: {f['sample'] or f['family']}: {f['message'][:160]}")

    em = dict(
        epsilon=args.epsilon,
        init=args.init,
        seed=args.seed,
        min_genotype_frequency=args.min_genotype_frequency,
        max_completions=args.max_completions,
        threads=args.threads,
    )
    est, phases = estimate_families(data, nom, res, use_family_phase=not args.no_family_phase, **em)
    statuses = {s: sum(p.status == s for p in phases) for s in ("complete", "partial", "failed")}
    say(f"Phasing: {statuses['complete']} complete, {statuses['partial']} partial, {statuses['failed']} failed")
    _tsv(
        out / "phases.tsv",
        [
            "family",
            "father",
            "mother",
            "status",
            "father_hap1",
            "father_hap2",
            "mother_hap1",
            "mother_hap2",
            "inconsistent_loci",
            "ambiguous_loci",
            "notes",
        ],
        (
            [
                p.family,
                p.father,
                p.mother,
                p.status,
                *p.haplotypes,
                ",".join(p.inconsistent_loci),
                ",".join(p.ambiguous_loci),
                "; ".join(p.notes),
            ]
            for p in phases
        ),
    )
    _tsv(
        out / "children.tsv",
        ["family", "child", "paternal_haplotype", "maternal_haplotype"],
        ([p.family, c, pat, mat] for p in phases for c, pat, mat in p.children),
    )

    est.write_frequencies(str(out / "hfs.dat"), cut=args.cut, renormalize=not args.no_renormalize)
    est.write_history(str(out / "epsilon.dat"))
    _tsv(out / "diplotypes.tsv", ["sample", "haplotype1", "haplotype2", "probability"], est.best_diplotypes())
    say(
        f"EM from {est.n_samples} founders: {len(est.haplotypes)} haplotypes, {est.iterations} iterations, "
        f"converged={est.converged}"
    )

    # Check the EM itself: founders as unrelated people vs. what the families show.
    unphased, _ = estimate_families(data, nom, res, use_family_phase=False, **em)
    cmp = compare_to_families(unphased, phases)
    (out / "comparison.json").write_text(json.dumps(cmp, indent=2) + "\n")
    if cmp["n_compared"]:
        say(
            f"EM vs families: {cmp['n_exact']}/{cmp['n_compared']} founders' haplotype pairs identical, "
            f"{cmp['mean_locus_agreement']:.1%} of loci agree; frequency difference (TVD) "
            f"{cmp['total_variation']:.3f}"
        )

    run_info = {
        "tool": "hlahaplotools families",
        "version": _version(),
        "input": str(args.input),
        "loci_and_resolutions": res,
        "imgthla_release": manifest.get("imgthla", {}).get("release") or nom.source,
        "mac_last_updated": manifest.get("mac", {}).get("last_updated"),
        "use_family_phase": not args.no_family_phase,
        "parameters": {**em, "cut": args.cut, "renormalize": not args.no_renormalize},
        "people": len(data),
        "families": len(fams),
        "phasing": statuses,
        "mendel_problems": len(mendel),
        "pedigree_flags": by_kind,
        "founders_used": est.n_samples,
        "iterations": est.iterations,
        "converged": est.converged,
        "log_likelihood": est.log_likelihood,
        "em_vs_family": {k: cmp[k] for k in ("n_compared", "n_exact", "mean_locus_agreement", "total_variation")},
        "seconds": round(time.perf_counter() - t0, 3),
        "left_out": [{"line": l, "id": i, "reason": r} for l, i, r in data.issues]
        + [{"id": i, "reason": r} for _, i, r in est.issues],
    }
    (out / "run.json").write_text(json.dumps(run_info, indent=2) + "\n")
    say(f"Wrote results to {out}/")
    return 0 if est.converged else 2


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="hlahaplotools",
        description="HLA haplotype frequency estimation (Hapl-o-Mat compatible, with missing-data handling).",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {_version()}")
    sub = parser.add_subparsers(dest="command", required=True)

    def data_options(p):
        p.add_argument("--data-dir", help="folder with manifest.json (default: $HLAHAPLOTOOLS_DATA or ./data)")
        p.add_argument("--nomenclature", help="lookup tables folder (default: from the manifest)")
        p.add_argument("--mac-codes", help="NMDP alpha.v3.txt (default: from the manifest)")

    def run_options(p):
        p.add_argument(
            "--max-completions",
            type=int,
            default=10_000,
            help="most haplotypes created to fill in one partly typed pattern (default 10000)",
        )
        p.add_argument("--threads", type=int, help="worker threads (default: all cores; results are identical)")
        p.add_argument("-q", "--quiet", action="store_true")

    def em_options(p):
        p.add_argument("--epsilon", type=float, default=1e-6, help="stop when no frequency changes more (default 1e-6)")
        p.add_argument("--init", choices=INITS, default="numberOccurrence")
        p.add_argument("--seed", type=int, default=0)
        p.add_argument("--cut", type=float, default=1e-6, help="drop haplotypes at or below this frequency")
        p.add_argument("--no-renormalize", action="store_true", help="do not rescale kept frequencies to sum to 1")
        p.add_argument(
            "--min-genotype-frequency",
            type=float,
            default=1e-5,
            help="leave out samples with more than 1/this possible genotypes (default 1e-5)",
        )

    e = sub.add_parser("estimate", help="estimate haplotype frequencies")
    e.add_argument("input", help="MAC (.dat/.mac) or GLSC (.glc/.glsc) file")
    e.add_argument("--loci", required=True, help="loci and resolutions, e.g. A:g,B:g,C:g,DRB1:g")
    e.add_argument("--format", choices=["mac", "glsc"], help="default: guess from the file extension")
    e.add_argument("-o", "--output", default="results", help="output folder (default: results)")
    em_options(e)
    data_options(e)
    run_options(e)
    e.set_defaults(func=cmd_estimate)

    f = sub.add_parser(
        "families",
        help="family data: Mendel check, phasing, founder-only estimate, EM-vs-family comparison",
        description=cmd_families.__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    f.add_argument("input", help="CSV/TSV with FAMILY_ID, Family_Member (F/M/C1..) and LOCUS_1/LOCUS_2 columns")
    f.add_argument("--loci", required=True, help="loci and resolutions, e.g. A:g,B:g,C:g,DRB1:g")
    f.add_argument("-o", "--output", default="results", help="output folder (default: results)")
    f.add_argument(
        "--no-family-phase",
        action="store_true",
        help="estimate from founders as unrelated people, ignoring the phase the children reveal",
    )
    f.add_argument(
        "--cross-family-duplicates",
        action="store_true",
        help="also flag people in different families with identical typing (noisy with common haplotypes)",
    )
    em_options(f)
    data_options(f)
    run_options(f)
    f.set_defaults(func=cmd_families)

    h = sub.add_parser("haplomat", help="run from a Hapl-o-Mat parameter file")
    h.add_argument("format", choices=["MAC", "GLSC", "mac", "glsc"], help="input format, as for Hapl-o-Mat")
    h.add_argument("parameters", help="Hapl-o-Mat parameter file, e.g. parametersMAC")
    data_options(h)
    run_options(h)
    h.set_defaults(func=cmd_haplomat)

    s = sub.add_parser("summary", help="summarise an input file without estimating")
    s.add_argument("input")
    s.add_argument("--format", choices=["mac", "glsc"])
    s.add_argument("--loci", help="comma-separated loci to keep")
    s.add_argument("--data-dir")
    s.add_argument("--mac-codes")
    s.set_defaults(func=cmd_summary)
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    sys.exit(args.func(args))
