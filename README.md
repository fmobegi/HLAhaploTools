# HLAhaploTools

> This repository previously held the HLAhaploTools **R package**; it now
> holds its successor, written in Rust with a Python API. The last R version
> is tagged [`r-package-v0.99.0`](https://github.com/fmobegi/HLAhaploTools/tree/r-package-v0.99.0).

The Python package and command are called `hlahaplotools`.

HLA haplotype frequency estimation by EM. A Rust port of
[Hapl-o-Mat](https://github.com/DKMS/Hapl-o-Mat) with a Python API and
command line, that also uses samples with **untyped or half-typed loci**
instead of discarding them.

- **Same results as Hapl-o-Mat** on fully typed data: within 1e-14 on
  Hapl-o-Mat's example files and on simulated data at every resolution
  (g, G, P, 1f–4f, asItIs). See [Compatibility](#compatibility).
- **Faster**: 16× on 20,000 ambiguous 6-locus samples (5.5 s vs 87 s, 8 cores).
- **Deterministic**: bit-identical output on any number of threads.
- **Traceable**: every run records the IMGT/HLA release and MAC table date.
- **Family data**: Mendel checks, exact phasing of parents from their
  children, founder-only estimation, and a check of the EM against the
  families. See [Families](#families).

## Setup (Linux / WSL)

Requires [uv](https://docs.astral.sh/uv/) and a Rust toolchain (`rustup`, Rust 1.85+).

```bash
uv run python scripts/build_nomenclature.py   # IMGT/HLA 3.66.0 + NMDP MAC table into data/ (~650 MB)
uv run hlahaplotools --help                   # first run compiles the Rust extension
```

The IMGT/HLA release is pinned in `PINNED_RELEASE` in
[scripts/build_nomenclature.py](scripts/build_nomenclature.py). Change it
deliberately (`--release 3.65.0`) and report the release with your results.
The NMDP MAC table has no release number; it is downloaded once and kept, and
its `LAST UPDATED` date is recorded in `data/manifest.json`
(`--refresh-mac` updates it).

`data/` then holds:

| Path | Contents |
|---|---|
| `imgthla/<release>/` | hla_nom_g.txt, hla_nom_p.txt, hla_nom.txt, Allelelist.txt, Deleted_alleles.txt |
| `derived/<release>/` | AllAllelesExpanded.txt, LargeG.txt, P.txt, Smallg.txt — the lookup tables used for resolution, byte-identical to Hapl-o-Mat's `prepareData` output |
| `mac/alpha.v3.txt` | NMDP multiple-allele codes |
| `manifest.json` | release, source URLs, checksums, MAC date |

## Command line

```bash
# Estimate at g resolution; writes results/hfs.dat, results/epsilon.dat, results/run.json
uv run hlahaplotools estimate population.dat --loci A:g,B:g,C:g,DQB1:g,DRB1:g

# Drop-in for a Hapl-o-Mat parameter file
uv run hlahaplotools haplomat MAC parametersMAC

# Look at an input file without estimating
uv run hlahaplotools summary population.glc
```

Input formats are Hapl-o-Mat's **MAC** (tab-separated, two columns per locus,
`24:DXTW`-style typings) and **GLSC** (an ID then one GL string per locus).
Untyped alleles are empty cells (or `-`, `.`, `NA`, `?`) in MAC files; in GLSC
files a locus is untyped if it is left out, and half-typed if written
`A*01:01+`.

Outputs, in Hapl-o-Mat's formats:

| File | Contents |
|---|---|
| `hfs.dat` | `haplotype<TAB>frequency`, most frequent first, cut at `--cut` (1e-6) and renormalised |
| `epsilon.dat` | per iteration: largest frequency change, log-likelihood |
| `run.json` | IMGT/HLA release, MAC date, all parameters, sample counts, timings, and every sample left out with the reason |

Main options: `--epsilon` (1e-6), `--init numberOccurrence|perturbation|equal|random`,
`--seed`, `--min-genotype-frequency` (1e-5), `--max-completions` (10,000),
`--threads`. Run `hlahaplotools estimate --help` for all of them.

## Python

```python
import hlahaplotools as h

ds = h.read_mac("population.dat", mac_codes="data/mac/alpha.v3.txt")  # or h.read_glsc(...)
ds.sample(0)  # {'A': 'A*01:01+A*02:01/A*02:02', 'B': None, 'C': 'C*07:01+'}
ds.missing_counts()  # untyped samples per locus

nom = h.Nomenclature("data/derived/3.66.0")
est = h.estimate(ds, nom, {"A": "g", "B": "g", "C": "g"}, epsilon=1e-6, threads=8)
est.frequencies(cut=1e-6)[:5]  # [('A*01:01g~B*08:01g~C*07:01g', 0.083...), ...]
est.issues  # samples left out, with reasons
est.write_frequencies("hfs.dat")
```

Passing `mac_codes` as a path loads only the codes the input file uses (a few
seconds, a few MB); `MacCodes(path)` loads the whole table (~650 MB) for reuse.

## Families

Relatives share haplotypes, so treating family members as unrelated samples
counts the same haplotypes several times. For family data, use a table with
one row per person:

```text
FAMILY_ID,Family_Member,A_1,A_2,B_1,B_2,DRB1_1,DRB1_2
fam01,F,A*01:01,A*02:01,B*08:01,B*07:02,DRB1*03:01,DRB1*15:01
fam01,M,03:01,11:01,35:01,44:02,04:01,13:01
fam01,C1,A*01:01,A*03:01,B*08:01,B*35:01,DRB1*03:01,DRB1*04:01
```

CSV, TSV or semicolon-separated. Columns may be named `A_1`, `A.1`, `A1`,
`A_ALLELE_1` or `HLA-A_1`; alleles may be written `A*01:01`, `*01:01` or
`01:01`, with `/` for alternatives. Roles: `F`/`Father`, `M`/`Mother`,
`C1`, `C2`, ... (also `Child2`, `Sib1`, ...). Rows without a family ID are
used as unrelated people.

```bash
uv run hlahaplotools families families.csv --loci A:g,B:g,C:g,DRB1:g,DQB1:g -o results/
```

| File | Contents |
|---|---|
| `flags.tsv` | Pedigree problems (see [Pedigree checks](#pedigree-checks)) |
| `mendel.tsv` | Children whose typing cannot come from their parents, per locus |
| `phases.tsv` | Each family's four parental haplotypes; `LOCUS*?` where unknown |
| `children.tsv` | Which haplotype each child inherited from each parent |
| `hfs.dat`, `epsilon.dat` | Frequencies estimated from founders (parents) only |
| `diplotypes.tsv` | Most likely haplotype pair per founder |
| `comparison.json` | EM on founders as unrelated people vs. the family phasing |

How it works:

- **Mendel check** at the target resolution, so a parent typed at 2 fields
  and a child at 3 are compared fairly.
- **Phasing** tries every way the children could have inherited the
  parents' haplotypes (the first 8 children) and keeps the patterns that fit
  the most loci. An untyped parent is rebuilt from the children as far as
  they show it. Loci no pattern fits (typing error or recombination) and loci
  the family cannot resolve are reported, never guessed.
- **Estimation** uses founders only. Loci the family phased enter the EM
  with known phase; the rest use the founder's own typing; a rebuilt parent's
  unknown alleles are treated as untyped. A family without typed parents
  that cannot be phased contributes its first child.

```python
fam = h.read_family_table("families.csv")
h.audit_families(fam, nom, "g")  # pedigree problems (see below)
h.check_mendelian(fam, nom, "g")  # per-locus Mendel failures
phases = h.phase_families(fam, nom, "g")  # FamilyPhase per family
est, phases = h.estimate_families(fam, nom, "g")  # founders only, with family phase
unphased, _ = h.estimate_families(fam, nom, "g", use_family_phase=False)
h.compare_to_families(unphased, phases)  # how well plain EM recovers the families
est.best_diplotypes()  # (sample, haplotype 1, haplotype 2, probability)
```

On simulated data (1,000 families of 1–4 children, 10% of parents untyped,
DQB1 20% untyped; total variation distance to the true frequencies, three
replicates), using the family phase gave 0.080 / 0.089 / 0.085, against
0.086 / 0.092 / 0.092 for founders without phase and 0.086 / 0.093 / 0.092
when everyone is treated as unrelated (which also overstates the sample size
about twofold). On 200 simulated families, every fully phased parent
matched the true haplotypes.

### Pedigree checks

`audit_families` (and `flags.tsv`) reports members whose typing does not fit
their stated relationships. Each flag names the member, the relative
involved, the loci, and a plain-language explanation.

| Kind | Meaning |
|---|---|
| `parent_excluded` | Child shares no allele with a parent at 2+ loci: not that parent's child, or a sample swap |
| `single_locus_mismatch` | One locus does not fit: usually a typing error, null allele or mutation |
| `trio_inconsistent` | Each parent fits alone, but no allele from each explains the child |
| `possible_duplicate` | Child identical to a parent and not fitting the other parent (a copied row), or two identical parents |
| `not_full_sibling` | With a parent untyped, a child carries alleles the siblings' parents could not also have passed on. If any of several children could be the odd one out, one family-level flag names them all |
| `recombination` | Inheritance from one parent switches once along chromosome 6; the interval is given |
| `inheritance_mismatch` | Inheritance switches several times, or a child lacks alleles an untyped parent must have |
| `structure` | Unknown member roles, a second father or mother, untypable rows |

Children are compared with haplotypes phased from the *rest* of the family,
leaving out members already shown not to belong. A child identical to a
parent is only flagged when it also fails to fit the other parent, because
parent–child matches also happen naturally. Identical people across families
are only reported with `cross_family_duplicates=True` (`--cross-family-duplicates`),
since unrelated people with common haplotypes often match.

Detection on simulated families ([scripts/validate_audit.py](scripts/validate_audit.py);
2,000 families of 2–4 children, 6 loci, 15% of parents untyped, two seeds):

| Planted problem | Both parents typed: flagged | A parent untyped: flagged |
|---|---|---|
| Child from another father | 88–94% | 53–69% |
| Child's row copied from the mother | 81–87% | 88–97% |
| Recombinant child | 79–81% | 36–47% |
| One-allele typo | 98–100% | 50–61% |
| *Clean children (false alarms)* | *1.1–1.6%* | *6.8–7.7%* |

On the example pedigrees in [examples/families](examples/families) (9
families, 40 people, 18 loci), the only flags are for `Family3/C1`: the father
is excluded at 6 loci and the row is identical to the mother's at 13 of 14
loci, i.e. the mother's typing copied into the child's row. To reproduce, and
to test phasing by leaving members out:

```bash
uv run python scripts/validate_pedigrees.py examples/families/*.csv
```

## Missing data

Hapl-o-Mat discards any sample that lacks a requested locus. hlahaplotools
keeps it: the probability of a sample is summed over all haplotype pairs that
agree with what *was* typed. An untyped locus can take any allele; a
half-typed locus `a+?` requires `a` on at least one haplotype. Written as a
signed sum over partially specified haplotypes, this gives exact E and M
steps, and reduces to Hapl-o-Mat's for fully typed samples (the derivation is
in [estimate.rs](crates/hlahaplotools-core/src/estimate.rs)).

If a partly typed sample matches no haplotype seen elsewhere, the open loci
are filled in with every allele observed at them, up to `--max-completions`
new haplotypes per pattern.

Validation on simulated data ([scripts/validate_missing.py](scripts/validate_missing.py),
5,000 samples, 5 loci, 300 true haplotypes, 3 replicates; error is the total
variation distance to the true frequencies, lower is better):

| C and DQB1 untyped | Samples left if dropped | With missing data | Dropping incomplete | Fully typed (best case) |
|---|---|---|---|---|
| 5% | 4504 / 5000 | 0.048 | 0.051 | 0.048 |
| 20% | 3204 / 5000 | 0.050 | 0.062 | 0.049 |
| 40% | 1786 / 5000 | 0.051 | 0.080 | 0.048 |

With 15% of the remaining C and DQB1 typings also half-typed, errors were
0.046 vs 0.063 (10% untyped) and 0.051 vs 0.088 (30% untyped). With B, C and
DQB1 each 30% untyped: 0.054 vs 0.082.

## Compatibility

Checked with [scripts/compare_haplomat.py](scripts/compare_haplomat.py), which
runs the original C++ Hapl-o-Mat and hlahaplotools on the same input and
tables:

| Data | Loci : resolution | Largest difference |
|---|---|---|
| Hapl-o-Mat example, MAC | A:g, B:g, DRB1:1f | 0 |
| Hapl-o-Mat example, GLSC | A:g, C:g, DPB1:g, DQB1:g, DRB1:P | 0 |
| Simulated, 3,000 MAC | all g | 0 |
| Simulated, 3,000 MAC | A:G, B:P, C:1f, DQB1:2f, DRB1:G | 1e-14 |
| Simulated, 3,000 MAC | A:3f, B:4f, C:asItIs, DQB1:P, DRB1:g | 0 |
| Simulated, 20,000 GLSC, ambiguous | 6 loci at g | 2e-6 |

Where results differ, it is on data with heavy typing ambiguity, where the EM
has many equally likely solutions. Both programs follow the same iterations
to 12+ digits for the first ~50 steps; after that, rounding decides which
near-zero haplotypes Hapl-o-Mat-style pruning (frequency < 1e-14) removes, and
the two can settle on different, nearly equally likely solutions.

Known differences:

- **Not supported**: Hapl-o-Mat's GLS (GLID) and READ formats, the ambiguity
  filter (`DO_AMBIGUITYFILTER`), `WRITE_GENOTYPES` and analytics. The
  `haplomat` subcommand warns and ignores these settings.
- **Initialisation**: `numberOccurrence` (the default here) matches exactly.
  `perturbation` and `random` use a different random number generator, so
  they start from different (seeded, reproducible) values.
- **Log-likelihood** values differ from Hapl-o-Mat's by a constant (it counts
  repeated phase combinations); the trajectories are the same.
- **`XX` codes** such as `A*01:XX` are read as "any A*01 allele". NMDP's MAC
  table has no `XX` entry, so Hapl-o-Mat discards such samples.
- **Loci order** follows `--loci`; the `haplomat` subcommand sorts loci
  alphabetically, as Hapl-o-Mat does.

## Layout

| Path | Contents |
|---|---|
| `crates/hlahaplotools-core/` | Pure Rust: parsers, MAC codes, nomenclature tables, resolution, EM, families, output |
| `src/lib.rs` | PyO3 bindings (`hlahaplotools._core`) |
| `src/hlahaplotools/` | Python package and command line |
| `scripts/build_nomenclature.py` | Downloads the pinned data and builds the lookup tables |
| `scripts/derive_tables.py` | The lookup-table builder (port of Hapl-o-Mat's prepareData) |
| `scripts/simulate.py` | Simulates populations or families (`--format families`) from known haplotype frequencies |
| `scripts/compare_haplomat.py` | Runs the C++ original and hlahaplotools side by side |
| `scripts/validate_missing.py` | Missing-data validation against simulated truth |
| `scripts/validate_pedigrees.py` | Family phasing on real pedigrees: trios, parentage, one parent missing (leave-members-out) |
| `scripts/validate_audit.py` | Pedigree checks: detection of planted errors in simulated families |
| `examples/families/` | Example pedigrees (from the R package): 9 families, 40 people, 18 loci |
| `reference/Hapl-o-Mat/` | Clone of the original (not in git), used by the comparison tests |

## Development

```bash
cargo test -p hlahaplotools-core      # Rust unit tests (no Python needed)
uv run pytest                          # Python API, CLI, and comparison tests
uv run ruff check . && cargo clippy --all-targets
```

The comparison tests need the original built locally, with its tables made
from the same IMGT/HLA release:

```bash
git clone https://github.com/DKMS/Hapl-o-Mat reference/Hapl-o-Mat
make -C reference/Hapl-o-Mat
mkdir -p reference/Hapl-o-Mat/data
cp data/derived/3.66.0/*.txt reference/Hapl-o-Mat/data/   # identical to its prepareData output
# its MAC table format: alpha.v3.txt without headers, '*' columns, CRs or leading tabs
tail -n +3 data/mac/alpha.v3.txt | tr -d '*\r' | sed 's/^\t//' | grep -v '^\s*$' \
    > reference/Hapl-o-Mat/data/MultipleAlleleCodes.txt
```

(Alternatively run its own `prepareData/BuildData.py`, which downloads the
latest release rather than the pinned one.)

For benchmarks, the release profile in `Cargo.toml` uses `lto = "fat"` and
`codegen-units = 1`; `uv run` builds in release mode. Projects on WSL's
`/mnt/` drives run noticeably slower than inside the Linux file system
(e.g. `~/hlahaplotools`).

## License

GPL-3.0-or-later (see [LICENSE](LICENSE)). The estimation and resolution
logic closely follows [Hapl-o-Mat](https://github.com/DKMS/Hapl-o-Mat)
(© DKMS gGmbH, GPL-3.0), so this project uses the same license.
