"""Command-line interface: `hlahaplotools {estimate,families,stats,impute,haplomat,rename,summary}`."""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from hlahaplotools._core import (
    Nomenclature,
    hardy_weinberg,
    impute,
    linkage_disequilibrium,
    read_gls,
    read_glsc,
    read_hml,
    read_mac,
    read_ngsengine,
    sniff_xml,
)
from hlahaplotools.correction import adjust_fdr
from hlahaplotools.genotypes import read_genotypes, write_genotypes
from hlahaplotools.multistart import estimate_multistart

INITS = ["numberOccurrence", "perturbation", "equal", "random"]


def _version() -> str:
    try:
        return version("hlahaplotools")
    except PackageNotFoundError:
        return "unknown"


# ---------------------------------------------------------------- data lookup


def find_data(data_dir: str | None) -> tuple[Path | None, dict]:
    """Locate the data folder written by `hlahaplotools download-data`.

    Order: --data-dir, $HLAHAPLOTOOLS_DATA, ./data. Returns (folder, manifest).
    """
    for candidate in (data_dir, os.environ.get("HLAHAPLOTOOLS_DATA"), "data"):
        if candidate and (Path(candidate) / "manifest.json").exists():
            folder = Path(candidate)
            return folder, json.loads((folder / "manifest.json").read_text(encoding="utf-8"))
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
            "error: no nomenclature tables found. Run `hlahaplotools download-data` (once), "
            "or pass --nomenclature DIR (a data/derived/<release> folder or Hapl-o-Mat's data/)."
        )
    return nom, mac, manifest


# ---------------------------------------------------------------- running


FORMATS = ["mac", "glsc", "gls", "hml", "ngsengine", "genotypes"]


def load_nomenclature(nom_dir: Path, args=None) -> Nomenclature:
    """Lookup tables, reading serological antigens or pre-2010 names if asked (--serology, --v2-names)."""
    serology = bool(getattr(args, "serology", False))
    v2 = bool(getattr(args, "v2_names", False))
    try:
        return Nomenclature(str(nom_dir), serology=serology, v2_names=v2)
    except ValueError as e:
        sys.exit(f"error: {e}")


def xml_files(path: Path) -> list[Path]:
    """A folder of XML reports (e.g. one HML file per sample), or one file."""
    if path.is_dir():
        files = sorted(p for p in path.iterdir() if p.suffix.lower() in (".hml", ".xml") and p.is_file())
        if not files:
            sys.exit(f"error: no .hml or .xml files in {path}")
        return files
    return [path]


def guess_format(path: str) -> str:
    """Input format from the extension, or from the content of XML files."""
    p = Path(path)
    low = p.name.lower()
    if p.is_dir() or low.endswith((".hml", ".xml")):
        kinds = {sniff_xml(str(f)) or ("hml" if f.suffix.lower() == ".hml" else None) for f in xml_files(p)}
        if len(kinds) == 1 and None not in kinds:
            return kinds.pop()
        sys.exit(f"error: cannot tell the XML layout of {path} (found {sorted(map(str, kinds))}); pass --format")
    if low.endswith(".pull"):
        return "gls"
    if low.endswith((".glc", ".glsc")):
        return "glsc"
    return "mac"


def read_input(
    path: str,
    fmt: str | None,
    mac: Path | None,
    loci: list[str] | None,
    half_typed: str = "homozygous",
    drb345: bool = True,
    glid: str | None = None,
    include_unapproved: bool = False,
):
    fmt = fmt or guess_format(path)
    if fmt == "genotypes":
        ds, issues = read_genotypes(path, loci)
        for line, sid, reason in issues[:5]:
            print(f"warning: {path}:{line} {sid}: {reason}", file=sys.stderr)
        if len(issues) > 5:
            print(f"warning: {len(issues) - 5} more people left out of {path}", file=sys.stderr)
        return ds, fmt
    mac_codes = str(mac) if mac and mac.exists() else None
    common = dict(mac_codes=mac_codes, loci=loci, half_typed=half_typed, drb345=drb345)
    if fmt == "gls":
        glid = glid or str(Path(path).with_suffix(".glid"))
        if not Path(glid).exists():
            sys.exit(f"error: GLS input needs a GLID file; {glid} not found (pass --glid)")
        return read_gls(path, glid, **common), fmt
    if fmt in ("hml", "ngsengine"):
        files = [str(f) for f in xml_files(Path(path))]
        if fmt == "hml":
            return read_hml(files, **common), fmt
        return read_ngsengine(files, include_unapproved=include_unapproved, **common), fmt
    reader = read_glsc if fmt == "glsc" else read_mac
    return reader(path, **common), fmt


def describe_notes(notes: list, say, where: str | None = None) -> None:
    """Print loci the reader left untyped on purpose (e.g. results not approved)."""
    if not notes:
        return
    say(f"  {len(notes)} locus results left untyped while reading:")
    for _, sid, reason in notes[:10]:
        say(f"    {sid}: {reason}")
    if len(notes) > 10:
        say(f"    ... and {len(notes) - 10} more" + (f" (see {where})" if where else ""))


def describe_fills(fills: dict, say) -> None:
    """Print what reading filled in, per locus."""
    for locus, f in fills.items():
        parts = []
        if f["homozygous"]:
            parts.append(f"{f['homozygous']} single alleles copied (homozygous)")
        if f["absent_gene"]:
            parts.append(f"{f['absent_gene']} gene copies set absent from DRB1")
        if f["drb_completed"]:
            parts.append(f"{f['drb_completed']} genotypes completed from DRB1")
        if f["duplicate_read_as_one"]:
            parts.append(f"{f['duplicate_read_as_one']} repeated alleles read as one copy (DRB1 allows one)")
        if f["drb1_conflicts"]:
            parts.append(f"{f['drb1_conflicts']} typings with more copies than DRB1 allows (left as typed)")
        say(f"  filled at {locus}: " + "; ".join(parts))


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
    half_typed: str = "homozygous",
    drb345: bool = True,
    starts: int = 1,
    completion: str = "linked",
    glid: str | None = None,
    include_unapproved: bool = False,
    serology: bool = False,
    v2_names: bool = False,
    out_genotypes: Path | None = None,
) -> int:
    say = (lambda *a: None) if quiet else print
    loci = [l for l, _ in loci_res]
    t0 = time.perf_counter()
    ds, fmt = read_input(input_path, fmt, mac, loci, half_typed, drb345, glid, include_unapproved)
    nom = load_nomenclature(nom_dir, argparse.Namespace(serology=serology, v2_names=v2_names))
    t1 = time.perf_counter()
    say(f"Read {len(ds)} samples from {input_path} ({fmt.upper()}); {len(ds.issues)} lines unreadable")
    describe_fills(ds.fill_counts(), say)
    describe_notes(ds.notes, say, str(out_freqs.with_name("run.json")))
    missing, half = ds.missing_counts(), ds.one_allele_counts()
    for locus, res in loci_res:
        say(f"  {locus:6} -> {res:6} untyped: {missing[locus]:6}  one allele typed: {half[locus]:6}")

    est, start_summary = estimate_multistart(
        ds,
        nom,
        [r for _, r in loci_res],
        starts=starts,
        epsilon=epsilon,
        init=init,
        seed=seed,
        min_genotype_frequency=min_genotype_frequency,
        max_completions=max_completions,
        threads=threads,
        completion=completion,
        diplotypes=False,  # the command line does not write them; costly with many open patterns
    )
    t2 = time.perf_counter()
    if starts > 1:
        lls = [s[2] for s in start_summary]
        best = max(range(len(lls)), key=lambda k: lls[k])
        # Differences below 0.01 log units are convergence noise at the usual epsilon.
        lower = sum(l < lls[best] - 0.01 for l in lls)
        if lower == 0:
            say(f"Starts: {starts}; all reached the same fit (logL within 0.01 of {lls[best]:.6f})")
        else:
            say(
                f"Starts: {starts}; best logL {lls[best]:.6f} from start {best + 1} ({start_summary[best][0]}, "
                f"seed {start_summary[best][1]}); {lower} of {starts} starts ended lower, "
                f"the worst by {lls[best] - min(lls):.3f}"
            )
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
    if out_genotypes:
        out_genotypes.parent.mkdir(parents=True, exist_ok=True)
        write_genotypes(out_genotypes, ds, nom, [r for _, r in loci_res], min_genotype_frequency)
        say(f"Wrote each person's resolved genotypes to {out_genotypes}")
    if out_history:
        out_history.parent.mkdir(parents=True, exist_ok=True)
        est.write_history(str(out_history))
    kept = est.frequencies(cut=cut, renormalize=renormalize)

    run_info = {
        "tool": "hlahaplotools",
        "version": _version(),
        "input": str(input_path),
        **({"glid": glid or str(Path(input_path).with_suffix(".glid"))} if fmt == "gls" else {}),
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
            "starts": starts,
            "completion": completion,
            "serology": serology,
            "v2_names": v2_names,
        },
        "starts": [
            {"init": i, "seed": sd, "log_likelihood": ll, "iterations": it, "converged": cv}
            for i, sd, ll, it, cv in start_summary
        ],
        "filling": {"half_typed": half_typed, "drb345_from_drb1": drb345, "filled": ds.fill_counts()},
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
        "loci_left_untyped": [{"id": i, "reason": r} for _, i, r in ds.notes],
        "samples_left_out": [{"id": i, "reason": r} for _, i, r in est.issues],
    }
    run_path = out_freqs.with_name("run.json")
    run_path.write_text(json.dumps(run_info, indent=2) + "\n", encoding="utf-8", newline="\n")
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


def cmd_rename(args) -> int:
    """Rewrite a typing file with current IMGT/HLA allele names.

    Names that IMGT/HLA has deleted or renamed since the typing was reported
    (for example C*04:09N, now C*04:09L, or B*47:01:01:01, now B*47:01:01:03)
    are replaced with the current name given in Deleted_alleles.txt. Everything
    else is written unchanged. The changes are listed in <output>.renames.tsv.
    """
    from hlahaplotools.rename import load_renames, rename_file

    if args.imgt_dir:
        imgt = Path(args.imgt_dir)
    else:
        folder, manifest = find_data(args.data_dir)
        if folder is None or "imgthla" not in manifest:
            sys.exit("error: no IMGT/HLA files found; pass --imgt-dir or run `hlahaplotools download-data`")
        imgt = folder / "imgthla" / manifest["imgthla"]["release"]
    fmt = args.format
    if fmt is None:
        low = args.input.lower()
        fmt = "glsc" if low.endswith((".glc", ".glsc")) else "families" if low.endswith((".csv", ".tsv")) else "mac"
    renames = load_renames(imgt)
    out = Path(args.output)
    seen = rename_file(Path(args.input), out, fmt, renames)
    table = out.with_name(out.name + ".renames.tsv")
    with table.open("w", encoding="utf-8", newline="\n") as f:
        f.write("old\tnew\ttypings\n")
        for (old, new), n in sorted(seen.items(), key=lambda x: -x[1]):
            f.write(f"{old}\t{new}\t{n}\n")
    print(f"IMGT/HLA {imgt.name}: {len(renames)} retired names known")
    print(f"Wrote {out}: {sum(seen.values())} typings renamed ({len(seen)} names); changes in {table}")
    for (old, new), n in sorted(seen.items(), key=lambda x: -x[1]):
        print(f"  {old} -> {new}: {n}")
    return 0


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
        half_typed=args.half_typed,
        drb345=not args.no_drb345,
        starts=args.starts,
        completion=args.completion,
        glid=args.glid,
        include_unapproved=args.include_unapproved,
        serology=args.serology,
        v2_names=args.v2_names,
        out_genotypes=out_dir / "genotypes.dat" if args.write_genotypes else None,
    )


def read_parameters(path: str) -> dict[str, str]:
    params = {}
    for line in Path(path).read_text(encoding="utf-8-sig").splitlines():
        line = line.strip()
        if line and not line.startswith("#") and "=" in line:
            key, _, value = line.partition("=")
            params[key.strip()] = value.strip()
    return params


def cmd_haplomat(args) -> int:
    """Run from a Hapl-o-Mat parameter file (MAC, GLSC, GLS or READ)."""
    p = read_parameters(args.parameters)
    fmt = args.format.lower()
    if fmt == "read":
        # A genotype list written with WRITE_GENOTYPES: codes are already resolved.
        input_path, glid, fmt = p["FILENAME_INPUT"], None, "genotypes"
        first = next((l for l in Path(input_path).read_text(encoding="utf-8-sig").splitlines() if l.strip()), "")
        if not first:
            sys.exit(f"error: {input_path} is empty")
        p.setdefault(
            "LOCI_AND_RESOLUTIONS", ",".join(g.split("*", 1)[0] + ":asItIs" for g in first.split()[-1].split("^"))
        )
    elif fmt == "gls":
        input_path, glid = p["FILENAME_PULL"], p["FILENAME_GLID"]
        # Untyped loci (GL ID 0) are always kept as untyped, so this switch has no effect here.
        if p.get("RESOLVE_MISSING_GENOTYPES", "false").lower() == "true":
            print("note: RESOLVE_MISSING_GENOTYPES ignored; untyped loci are always kept", file=sys.stderr)
    else:
        input_path, glid = p["FILENAME_INPUT"], None
    for key, what in [
        ("DO_AMBIGUITYFILTER", "the ambiguity filter"),
        ("DO_ANALYTICS", "analytics"),
    ]:
        if p.get(key, "false").lower() == "true":
            print(f"warning: {key}=true ignored; {what} is not supported", file=sys.stderr)
    nom, mac, manifest = resolve_data_paths(args)
    # Hapl-o-Mat orders loci alphabetically (it keeps them in a std::map).
    # Older Hapl-o-Mat parameter files count digits: 2d, 4d, 6d, 8d = 1f, 2f, 3f, 4f.
    old = {"2d": "1f", "4d": "2f", "6d": "3f", "8d": "4f"}
    loci_res = sorted((l, old.get(r, r)) for l, r in parse_loci(p["LOCI_AND_RESOLUTIONS"]))
    write_geno = p.get("WRITE_GENOTYPES", "false").lower() == "true"
    return run_estimate(
        input_path=input_path,
        fmt=fmt,
        glid=glid,
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
        half_typed=args.half_typed,
        drb345=not args.no_drb345,
        serology=args.serology,
        v2_names=args.v2_names,
        out_genotypes=Path(p.get("FILENAME_GENOTYPES", "results/genotypes.dat")) if write_geno else None,
    )


def cmd_summary(args) -> int:
    if args.mac_codes:
        mac = Path(args.mac_codes)
    else:
        folder, manifest = find_data(args.data_dir)
        mac = folder / "mac" / "alpha.v3.txt" if folder and "mac" in manifest else None
    loci = args.loci.split(",") if args.loci else None
    ds, fmt = read_input(
        args.input, args.format, mac, loci, args.half_typed, not args.no_drb345, args.glid, args.include_unapproved
    )
    print(f"{ds!r} ({fmt.upper()})")
    describe_fills(ds.fill_counts(), print)
    describe_notes(ds.notes, print)
    missing, half = ds.missing_counts(), ds.one_allele_counts()
    for locus in ds.loci:
        print(f"  {locus:6} untyped: {missing[locus]:6}   one allele typed: {half[locus]:6}")
    for line, sid, reason in ds.issues[:20]:
        print(f"  line {line} (id {sid}): {reason}")
    if len(ds.issues) > 20:
        print(f"  ... and {len(ds.issues) - 20} more")
    return 0


# ---------------------------------------------------------------- argparse


def cmd_stats(args) -> int:
    """Hardy–Weinberg tests per locus and linkage disequilibrium per pair of loci.

    Hardy–Weinberg uses the people whose genotype is unambiguous at the target
    resolution (untyped, half-typed and ambiguous people are counted apart): observed
    and expected heterozygosity, F = 1 - Ho/He, and an exact test from random pairings
    of the same alleles. A blank second allele read as homozygous adds homozygotes; the
    number is reported, and --half-typed unknown leaves such people out instead.

    Linkage disequilibrium uses the two-locus haplotype frequencies estimated by EM from
    everyone typed at either locus: Hedrick's D', Wn (Cramér's V of the haplotype
    table), a likelihood-ratio test of linkage equilibrium (chi-square, approximate when
    many alleles are rare), and D' per allele pair. With --ld-permutations, linkage
    equilibrium is tested by permutation: one locus is shuffled across people and the
    two-locus likelihood re-estimated (the chi-square statistic is reported for reference;
    with many rare alleles its chi-square distribution does not hold).

    Both tables add q_value: the Benjamini–Hochberg false discovery rate adjusted P,
    across the loci (hwe.tsv) or the pairs (ld.tsv) of the run.
    """
    say = (lambda *a: None) if args.quiet else print
    nom_dir, mac, manifest = resolve_data_paths(args)
    loci_res = parse_loci(args.loci)
    loci = [l for l, _ in loci_res]
    t0 = time.perf_counter()
    ds, fmt = read_input(
        args.input, args.format, mac, loci, args.half_typed, not args.no_drb345, args.glid, args.include_unapproved
    )
    nom = load_nomenclature(nom_dir, args)
    res = [r for _, r in loci_res]
    say(f"Read {len(ds)} samples from {args.input} ({fmt.upper()}); {len(ds.issues)} lines unreadable")
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    hwe = hardy_weinberg(ds, nom, res, permutations=args.permutations, seed=args.seed, threads=args.threads)
    t1 = time.perf_counter()
    for r, q in zip(hwe, adjust_fdr([r["p_value"] for r in hwe]), strict=True):
        r["q_value"] = q
    cols = ["locus", "n", "alleles", "observed_heterozygosity", "expected_heterozygosity", "f", "p_value",
            "q_value", "untyped", "half_typed", "ambiguous", "homozygous_filled"]  # fmt: skip
    _tsv(out / "hwe.tsv", cols, ([r[c] for c in cols] for r in hwe))
    say(f"Hardy–Weinberg ({args.permutations:,} random tables per locus; hwe.tsv):")
    for r in hwe:
        p = "–" if r["p_value"] is None else f"{r['p_value']:.4f}  q={r['q_value']:.4f}"
        note = f"  ({r['homozygous_filled']} read as homozygous from a blank)" if r["homozygous_filled"] else ""
        ho, he = r["observed_heterozygosity"], r["expected_heterozygosity"]
        say(f"  {r['locus']:6} n={r['n']:6}  Ho={ho:.3f}  He={he:.3f}  F={r['f']:+.3f}  P={p}{note}")

    if args.pairs == "adjacent":
        pairs = list(zip(loci, loci[1:], strict=False))
    else:
        pairs = [(a, b) for i, a in enumerate(loci) for b in loci[i + 1 :]]
    ld = linkage_disequilibrium(
        ds, nom, res, pairs, epsilon=args.epsilon, seed=args.seed, starts=args.starts,
        min_genotype_frequency=args.min_genotype_frequency, max_completions=args.max_completions,
        completion=args.completion, threads=args.threads, permutations=args.ld_permutations,
    )  # fmt: skip
    t2 = time.perf_counter()
    for r, q in zip(ld, adjust_fdr([r["p_value"] for r in ld]), strict=True):
        r["q_value"] = q
    cols = ["locus_a", "locus_b", "n", "alleles_a", "alleles_b", "d_prime", "wn", "chi2", "df", "p_value", "q_value"]
    _tsv(out / "ld.tsv", cols, ([r[c] for c in cols] for r in ld))
    pcols = ["allele_a", "allele_b", "haplotype_frequency", "frequency_a", "frequency_b", "d", "d_prime"]
    _tsv(out / "ld_alleles.tsv", ["locus_a", "locus_b", *pcols],
         ([r["locus_a"], r["locus_b"], *[p[c] for c in pcols]] for r in ld for p in r["pairs"]))  # fmt: skip
    perm = f"; P from {args.ld_permutations:,} permutations" if args.ld_permutations else ""
    say(f"Linkage disequilibrium (ld.tsv; per allele pair in ld_alleles.tsv{perm}):")
    for r in ld:
        p = "–" if r["p_value"] is None else f"{r['p_value']:.4f}  q={r['q_value']:.4f}"
        say(f"  {r['locus_a']}~{r['locus_b']:6} n={r['n']:6}  D'={r['d_prime']:.3f}  Wn={r['wn']:.3f}"
            f"  P={p}")  # fmt: skip

    info = {
        "tool": "hlahaplotools",
        "version": _version(),
        "input": str(args.input),
        "format": fmt,
        "loci_and_resolutions": dict(loci_res),
        "imgthla_release": manifest.get("imgthla", {}).get("release") or nom.source,
        "mac_last_updated": manifest.get("mac", {}).get("last_updated"),
        "parameters": {
            "permutations": args.permutations,
            "ld_permutations": args.ld_permutations,
            "seed": args.seed,
            "pairs": args.pairs,
            "starts": args.starts,
            "epsilon": args.epsilon,
            "completion": args.completion,
            "half_typed": args.half_typed,
            "drb345_from_drb1": not args.no_drb345,
        },  # fmt: skip
        "samples_read": len(ds),
        "unreadable_lines": [{"line": l, "id": i, "reason": r} for l, i, r in ds.issues],
        "seconds": {"hardy_weinberg": round(t1 - t0, 3), "linkage_disequilibrium": round(t2 - t1, 3)},
    }
    (out / "run.json").write_text(json.dumps(info, indent=2) + "\n", encoding="utf-8", newline="\n")
    say(f"Wrote hwe.tsv, ld.tsv, ld_alleles.tsv and run.json to {out}")
    return 0


def read_frequency_table(path: Path) -> list[tuple[str, float]]:
    """hfs.dat: haplotype<TAB>frequency per line."""
    out = []
    for n, line in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
        if not line.strip():
            continue
        parts = line.split()
        try:
            out.append((parts[0], float(parts[1])))
        except IndexError, ValueError:
            sys.exit(f"error: {path}:{n}: expected 'haplotype<TAB>frequency', got {line!r}")
    return out


def cmd_impute(args) -> int:
    """Each person's most likely haplotype pairs and genotypes, from known haplotype frequencies.

    Every pair of haplotypes in the frequency table that fits a person's typing gets the
    probability f(h1) f(h2) (both orders when they differ), times the weight of the
    typings it fits, normalised over the person's pairs: the imputation of GRIMM and
    py-graph-imputation, here over a flat frequency table. Untyped or half-typed loci are
    imputed from the haplotypes that fit the rest. With several --frequencies NAME=FILE,
    each pair comes from one population, weighted by its --prior (default equal), and each
    person's population posteriors are written too.

    Loci and resolutions default to those in run.json next to the (first) frequency file,
    so the typing is resolved as it was for the estimate. People whose typing no pair in
    the table fits are listed in run.json.

    Output: pairs.tsv (id, rank, population, haplotype_1, haplotype_2, probability),
    genotypes.tsv (id, rank, genotype, probability), populations.tsv (with several
    populations) and run.json.
    """
    say = (lambda *a: None) if args.quiet else print
    pops = []
    for spec in args.frequencies:
        name, _, path = spec.rpartition("=")
        path = Path(path)
        if not path.exists():
            sys.exit(f"error: frequency file {path} not found")
        pops.append((name or path.parent.name or path.stem, path))
    if len({n for n, _ in pops}) != len(pops):
        sys.exit("error: give each --frequencies a distinct NAME=FILE")
    priors = {}
    for spec in args.prior or []:
        name, _, w = spec.partition("=")
        priors[name] = float(w)
    unknown = set(priors) - {n for n, _ in pops}
    if unknown:
        sys.exit(f"error: --prior for unknown population(s) {sorted(unknown)}")
    if args.loci:
        loci_res = parse_loci(args.loci)
    else:
        run = pops[0][1].with_name("run.json")
        if not run.exists():
            sys.exit(f"error: no --loci given and no {run} to take them from")
        loci_res = list(json.loads(run.read_text(encoding="utf-8"))["loci_and_resolutions"].items())
    loci = [l for l, _ in loci_res]
    nom_dir, mac, manifest = resolve_data_paths(args)
    t0 = time.perf_counter()
    ds, fmt = read_input(
        args.input, args.format, mac, loci, args.half_typed, not args.no_drb345, args.glid, args.include_unapproved
    )
    nom = load_nomenclature(nom_dir, args)
    tables = [(name, read_frequency_table(path), priors.get(name, 1.0)) for name, path in pops]
    say(f"Read {len(ds)} samples from {args.input} ({fmt.upper()}); "
        + ", ".join(f"{n}: {len(t):,} haplotypes" for n, t, _ in tables))  # fmt: skip
    try:
        rows, issues = impute(
            ds, nom, tables, [r for _, r in loci_res], top=args.top, max_candidates=args.max_candidates,
            threads=args.threads,
        )  # fmt: skip
    except ValueError as e:
        sys.exit(f"error: {e}")
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    _tsv(out / "pairs.tsv", ["id", "rank", "population", "haplotype_1", "haplotype_2", "probability"],
         ([r["id"], k, *p] for r in rows for k, p in enumerate(r["pairs"], 1)))  # fmt: skip
    _tsv(out / "genotypes.tsv", ["id", "rank", "genotype", "probability"],
         ([r["id"], k, *g] for r in rows for k, g in enumerate(r["genotypes"], 1)))  # fmt: skip
    if len(pops) > 1:
        _tsv(out / "populations.tsv", ["id", "population", "posterior"],
             ([r["id"], n, p] for r in rows for n, p in r["populations"]))  # fmt: skip
    approx = sum(r["approximate"] for r in rows)
    top1 = [r["pairs"][0][3] for r in rows]
    info = {
        "tool": "hlahaplotools",
        "version": _version(),
        "input": str(args.input),
        "format": fmt,
        "frequencies": {n: str(p) for n, p in pops},
        "priors": {n: priors.get(n, 1.0) for n, _ in pops},
        "loci_and_resolutions": dict(loci_res),
        "imgthla_release": manifest.get("imgthla", {}).get("release") or nom.source,
        "parameters": {
            "top": args.top,
            "max_candidates": args.max_candidates,
            "serology": args.serology,
            "v2_names": args.v2_names,
            "half_typed": args.half_typed,
        },  # fmt: skip
        "samples_read": len(ds),
        "samples_imputed": len(rows),
        "samples_approximate": approx,
        "unreadable_lines": [{"line": l, "id": i, "reason": r} for l, i, r in ds.issues],
        "samples_left_out": [{"id": i, "reason": r} for _, i, r in issues],
        "seconds": round(time.perf_counter() - t0, 3),
    }
    (out / "run.json").write_text(json.dumps(info, indent=2) + "\n", encoding="utf-8", newline="\n")
    if rows:
        top1.sort()
        say(f"Imputed {len(rows):,} of {len(ds):,} people; median probability of the best pair "
            f"{top1[len(top1) // 2]:.3f}; {len(issues)} left out (run.json)")  # fmt: skip
    else:
        say(f"No one could be imputed; {len(issues)} left out (run.json)")
    if approx:
        say(
            f"  {approx} people had more than {args.max_candidates:,} candidate haplotypes; the most frequent were used"
        )
    say(f"Wrote pairs.tsv and genotypes.tsv to {out}")
    return 0


def _tsv(path: Path, header: list[str], rows) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(header) + "\n")
        for row in rows:
            f.write("\t".join("" if v is None else str(v) for v in row) + "\n")


SERIOUS = {"parent_excluded", "possible_duplicate", "not_full_sibling", "trio_inconsistent", "inheritance_mismatch"}


def cmd_families(args) -> int:
    """Family data workflow. Writes to the output folder:

    flags.tsv        pedigree problems: excluded parents, duplicated rows,
                     non-siblings, recombinants, single-locus mismatches
    mendel.tsv       children whose typing cannot come from their parents
    phases.tsv       each family's parental haplotypes: A and B (father), C and D
                     (mother); A and C are the lowest-numbered child's. LOCUS*? = unknown
    children.tsv     each child's haplotypes: A or B from the father, C or D from the mother
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
    data = read_family_table(
        args.input,
        mac_codes=str(mac) if mac and mac.exists() else None,
        loci=loci,
        half_typed=args.half_typed,
        drb345=not args.no_drb345,
    )
    nom = load_nomenclature(nom_dir, args)
    fams = data.families()
    say(
        f"Read {len(data)} people in {len(fams)} families ({len(data.unrelated)} unrelated); "
        f"{len(data.issues)} rows or members left out"
    )
    for line, sid, reason in data.issues[:10]:
        say(f"  {sid or f'line {line}'}: {reason}")
    describe_fills(data.fill_counts(), say)

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
        completion=args.completion,
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
            "hap_A",
            "hap_B",
            "hap_C",
            "hap_D",
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
        ["family", "child", "paternal", "maternal", "haplotypes"],
        (
            [p.family, c, pat or "?", mat or "?", f"{pat or '?'}{mat or '?'}"]
            for p in phases
            for c, pat, mat in p.children
        ),
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
    (out / "comparison.json").write_text(json.dumps(cmp, indent=2) + "\n", encoding="utf-8", newline="\n")
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
        "filling": {
            "half_typed": args.half_typed,
            "drb345_from_drb1": not args.no_drb345,
            "filled": data.fill_counts(),
        },
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
    (out / "run.json").write_text(json.dumps(run_info, indent=2) + "\n", encoding="utf-8", newline="\n")
    say(f"Wrote results to {out}/")
    return 0 if est.converged else 2


def _prog() -> str:
    """The command name as typed: hlahaplo is the short form of hlahaplotools."""
    name = Path(sys.argv[0]).stem if sys.argv and sys.argv[0] else ""
    return "hlahaplo" if name == "hlahaplo" else "hlahaplotools"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=_prog(),
        description="Extended HLA and MHC haplotypes in families and populations: pedigree phasing, and EM "
        "frequency estimation that keeps untyped and half-typed loci (Hapl-o-Mat compatible).",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {_version()}")
    sub = parser.add_subparsers(dest="command", required=True)

    def input_options(p):
        p.add_argument("--format", choices=FORMATS, help="default: guess from the extension or XML content")
        p.add_argument("--glid", help="GLS input: the GLID file (default: the pull file with a .glid suffix)")
        p.add_argument(
            "--include-unapproved",
            action="store_true",
            help="NGSengine exports: also read locus results not approved in review (default: left untyped)",
        )

    def data_options(p):
        p.add_argument("--data-dir", help="folder with manifest.json (default: $HLAHAPLOTOOLS_DATA or ./data)")
        p.add_argument("--nomenclature", help="lookup tables folder (default: from the manifest)")
        p.add_argument("--mac-codes", help="NMDP alpha.v3.txt (default: from the manifest)")
        p.add_argument(
            "--serology",
            action="store_true",
            help="read typings that are serological antigens (B*62, DRB1*15) as their alleles (WMDA tables)",
        )
        p.add_argument(
            "--v2-names", action="store_true", help="read pre-2010 allele names (A*0201, A*010101) as current names"
        )
        fill_options(p)

    def fill_options(p):
        p.add_argument(
            "--half-typed",
            choices=["homozygous", "unknown"],
            default="homozygous",
            help="one allele typed, the other blank: homozygous copies it across (reporting convention, "
            "default); unknown keeps the second allele unknown",
        )
        p.add_argument(
            "--no-drb345",
            action="store_true",
            help="do not complete DRB3/4/5 from DRB1 (absent genes, hemizygous typings)",
        )

    def run_options(p):
        p.add_argument(
            "--max-completions",
            type=int,
            default=10_000,
            help="most haplotypes created to fill in one partly typed pattern (default 10000)",
        )
        p.add_argument(
            "--completion",
            choices=["linked", "all"],
            default="linked",
            help="how untyped loci are filled in: alleles seen with the most strongly linked typed locus "
            "(linked, default) or every combination of the alleles seen (all)",
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
    e.add_argument(
        "input",
        help="MAC (.dat/.mac), GLSC (.glc/.glsc), GLS pull file (.pull, with --glid), HML (.hml), "
        "NGSengine XML export, or a folder of HML/XML files",
    )
    e.add_argument("--loci", required=True, help="loci and resolutions, e.g. A:g,B:g,C:g,DRB1:g")
    input_options(e)
    e.add_argument("-o", "--output", default="results", help="output folder (default: results)")
    e.add_argument(
        "--write-genotypes",
        action="store_true",
        help="also write each person's resolved genotypes to genotypes.dat (Hapl-o-Mat's WRITE_GENOTYPES; "
        "read back with --format genotypes)",
    )
    em_options(e)
    e.add_argument(
        "--starts",
        type=int,
        default=10,
        help="EM start values: --init first, then perturbed starts (seed+1, ...); the best fit is kept "
        "(default 10; 1 reproduces a single Hapl-o-Mat run)",
    )
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

    st = sub.add_parser(
        "stats",
        help="Hardy–Weinberg tests per locus and linkage disequilibrium per pair of loci",
        description=cmd_stats.__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    st.add_argument("input", help="typing file, as for estimate")
    st.add_argument("--loci", required=True, help="loci and resolutions, e.g. A:2f,B:2f,C:2f,DRB1:2f")
    input_options(st)
    st.add_argument("-o", "--output", default="results/stats", help="output folder (default: results/stats)")
    st.add_argument("--permutations", type=int, default=10_000, help="random tables per Hardy–Weinberg test")
    st.add_argument(
        "--pairs",
        choices=["all", "adjacent"],
        default="all",
        help="locus pairs for linkage disequilibrium: every pair (default) or neighbours in --loci order",
    )
    st.add_argument("--starts", type=int, default=10, help="EM starts for each two-locus estimate (default 10)")
    st.add_argument(
        "--ld-permutations",
        type=int,
        default=0,
        help="permutation test of linkage equilibrium per pair: datasets with one locus shuffled "
        "(each needs an EM; default 0, no P-value)",
    )
    em_options(st)
    data_options(st)
    run_options(st)
    st.set_defaults(func=cmd_stats)

    im = sub.add_parser(
        "impute",
        help="each person's most likely haplotype pairs and genotypes, from known frequencies",
        description=cmd_impute.__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    im.add_argument("input", help="typing file, as for estimate")
    im.add_argument(
        "--frequencies",
        action="append",
        required=True,
        metavar="[NAME=]FILE",
        help="haplotype frequencies (hfs.dat); repeat with NAME=FILE for several populations",
    )
    im.add_argument("--prior", action="append", metavar="NAME=WEIGHT", help="population prior (default: equal)")
    im.add_argument("--loci", help="loci and resolutions (default: from run.json next to the frequency file)")
    input_options(im)
    im.add_argument("-o", "--output", default="results/impute", help="output folder (default: results/impute)")
    im.add_argument("--top", type=int, default=10, help="pairs and genotypes written per person (default 10)")
    im.add_argument(
        "--max-candidates",
        type=int,
        default=5000,
        help="most candidate haplotypes per person and population; above it the most frequent are kept (default 5000)",
    )
    data_options(im)
    im.add_argument("--threads", type=int, help="worker threads (default: all cores)")
    im.add_argument("-q", "--quiet", action="store_true")
    im.set_defaults(func=cmd_impute)

    h = sub.add_parser("haplomat", help="run from a Hapl-o-Mat parameter file")
    h.add_argument(
        "format",
        choices=["MAC", "GLSC", "GLS", "READ", "mac", "glsc", "gls", "read"],
        help="input format, as for Hapl-o-Mat",
    )
    h.add_argument("parameters", help="Hapl-o-Mat parameter file, e.g. parametersMAC")
    data_options(h)
    run_options(h)
    h.set_defaults(func=cmd_haplomat)

    r = sub.add_parser(
        "rename",
        help="rewrite a typing file with current IMGT/HLA names (deleted or renamed alleles updated)",
        description=cmd_rename.__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    r.add_argument("input", help="MAC, GLSC or family table")
    r.add_argument("output", help="file to write; a <output>.renames.tsv table of changes is written beside it")
    r.add_argument("--format", choices=["mac", "glsc", "families"], help="default: guess from the extension")
    r.add_argument("--data-dir")
    r.add_argument("--imgt-dir", help="IMGT/HLA folder with Allelelist.txt and Deleted_alleles.txt")
    r.set_defaults(func=cmd_rename)

    s = sub.add_parser("summary", help="summarise an input file without estimating")
    s.add_argument("input")
    input_options(s)
    s.add_argument("--loci", help="comma-separated loci to keep")
    s.add_argument("--data-dir")
    s.add_argument("--mac-codes")
    fill_options(s)
    s.set_defaults(func=cmd_summary)

    from hlahaplotools import reference_data

    dd = sub.add_parser(
        "download-data",
        help="download the IPD-IMGT/HLA files and NMDP MAC table, and build the lookup tables (run once)",
    )
    reference_data.add_arguments(dd)
    dd.set_defaults(func=reference_data.run)
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    sys.exit(args.func(args))
