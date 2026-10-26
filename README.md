# HLAhaploTools

HLAhaploTools estimates HLA haplotype frequencies from population and family typing data with an expectation–maximisation (EM) algorithm. It has a Rust core, a Python API and a command line (`hlahaplotools`).

## Key strengths

- **Ambiguous typings**: NMDP multiple-allele codes (`24:DXTW`), allele lists (`02:01/02:02`), GL strings, and g, G and P groups, resolved against a pinned IPD-IMGT/HLA release.
- **Any target resolution** per locus: first to fourth field, g, G or P groups, or alleles exactly as typed.
- **Incomplete typing**: samples with an untyped locus, or with only one allele typed at a locus, are kept in the estimate instead of discarded.
- **Family data**: Mendel checks, phasing of parents from their children (including untyped parents), estimation from founders only, and checks that flag non-paternity, sample swaps, copied rows, non-siblings and recombinants.
- **Reproducible**: results are bit-identical on any number of threads, and every run records the IMGT/HLA release, MAC table date and all parameters.

## Contents

- [Installation](#installation)
- [Reference data](#reference-data)
- [Quick start](#quick-start)
- [Input formats](#input-formats)
- [Resolutions](#resolutions)
- [Command line](#command-line)
- [Output files](#output-files)
- [Python API](#python-api)
- [Methods](#methods)
- [Validation](#validation)
- [Limitations](#limitations)
- [Development](#development)
- [Citation](#citation)
- [License and acknowledgements](#license-and-acknowledgements)

## Installation

Requirements: Linux, macOS, or Windows with WSL; Python 3.14 or newer; [uv](https://docs.astral.sh/uv/); and a Rust toolchain 1.85 or newer ([rustup](https://rustup.rs/)).

```bash
git clone https://github.com/fmobegi/HLAhaploTools.git
cd HLAhaploTools
uv sync                       # builds the Rust extension and installs the package
uv run hlahaplotools --version
```

The release build uses link-time optimisation, so the first build takes a few minutes.

## Reference data

Typings are resolved against the IPD-IMGT/HLA nomenclature and the NMDP multiple-allele-code (MAC) table. Download them once:

```bash
uv run python scripts/build_nomenclature.py     # about 650 MB, mostly the MAC table
```

This writes `data/`:

| Path | Contents |
| --- | --- |
| `imgthla/<release>/` | `hla_nom_g.txt`, `hla_nom_p.txt`, `hla_nom.txt`, `Allelelist.txt`, `Deleted_alleles.txt` |
| `derived/<release>/` | Lookup tables used for resolution: `AllAllelesExpanded.txt`, `LargeG.txt`, `P.txt`, `Smallg.txt` |
| `mac/alpha.v3.txt` | NMDP multiple-allele codes |
| `manifest.json` | Release, source URLs, SHA-256 checksums, MAC table date |

The IMGT/HLA release is pinned (`PINNED_RELEASE` in [scripts/build_nomenclature.py](scripts/build_nomenclature.py), currently 3.66.0). To use another release, run with `--release 3.65.0`. The MAC table has no release number: it is downloaded once and kept, and its `LAST UPDATED` date is recorded; `--refresh-mac` downloads a new copy. Report both with any results.

Commands find `data/` in the current folder, in `$HLAHAPLOTOOLS_DATA`, or with `--data-dir`.

## Quick start

### Population data

[examples/population/population.dat](examples/population/population.dat) holds 1,000 simulated people typed at A, C, B, DRB1 and DQB1, drawn from the known frequencies in [truth.tsv](examples/population/truth.tsv). DQB1 is untyped in about 10% of them, and one C allele is missing in about 5%.

```bash
uv run hlahaplotools estimate examples/population/population.dat \
    --loci A:g,C:g,B:g,DRB1:g,DQB1:g -o results/population
```

Expected output:

```text
Read 1000 samples from examples/population/population.dat (MAC); 0 lines unreadable
  A      -> g      untyped:      0  one allele typed:      0
  C      -> g      untyped:      0  one allele typed:     43
  B      -> g      untyped:      0  one allele typed:      0
  DRB1   -> g      untyped:      0  one allele typed:      0
  DQB1   -> g      untyped:     93  one allele typed:      0
EM: 1000 samples (131 incomplete), 620 genotypes, 13588 candidate haplotypes (4986 filled in), 53 iterations, converged=True, logL=-6903.436695
Wrote 180 haplotypes to results/population/hfs.dat (details in results/population/run.json)
```

The most frequent haplotypes in `results/population/hfs.dat`, with the true frequencies used in the simulation:

| Haplotype (g resolution) | Estimated | True |
| --- | --- | --- |
| `A*02:1192~C*08:150~B*53:31~DRB1*08:81~DQB1*05:66` | 0.2265 | 0.2128 |
| `A*11:81~C*03:03g~B*35:61~DRB1*13:165~DQB1*04:24` | 0.0989 | 0.0993 |
| `A*32:136~C*15:05g~B*52:82~DRB1*04:231~DQB1*05:01g` | 0.0665 | 0.0636 |
| `A*02:313~C*07:1118~B*14:55~DRB1*14:95~DQB1*03:122` | 0.0430 | 0.0463 |
| `A*30:253~C*12:03g~B*27:05g~DRB1*14:129~DQB1*03:122` | 0.0295 | 0.0362 |

Differences of this size are sampling noise for 1,000 people.

### Family data

[examples/families/](examples/families) holds two family tables: 9 families, 40 people, typed at up to 18 loci.

```bash
uv run hlahaplotools families examples/families/family_typing_data.csv \
    --loci A:2f,C:2f,B:2f,DRB1:2f,DQA1:2f,DQB1:2f,DPA1:2f,DPB1:2f -o results/families
```

Expected output:

```text
Read 14 people in 3 families (0 unrelated); 0 rows or members left out
Mendel check: 6 problems (mendel.tsv)
Pedigree checks: 2 flags (flags.tsv): {'parent_excluded': 1, 'trio_inconsistent': 1}
  parent_excluded: Family3/C1: shares no allele with the father at 5 of 8 loci compared: not the biological father, or a sample swap (...)
  trio_inconsistent: Family3/C1: each parent fits alone but not both together at 1 loci: typing error, null allele or mutation (...)
Phasing: 2 complete, 1 partial, 0 failed
EM from 6 founders: 264 haplotypes, 1 iterations, converged=True
...
```

`Family3/C1` does not fit the father at 5 loci. Its typing matches the mother's at every locus except A: the mother's row appears to have been copied into the child's. Families 1 and 2 are phased completely. Family3 is phased only partly, because the flagged child is still included; correct or remove flagged members and rerun for clean phasing and estimates.

The other example table, `validated.csv`, has 6 families and no flags.

## Input formats

Missing typings are always allowed: empty cells, or `-`, `.`, `NA`, `N/A`, `?`. A locus with neither allele typed is *untyped*; a locus with one allele typed is *half-typed*. Both are used in the estimate (see [Methods](#methods)).

### MAC table (`.dat`, `.mac`)

Tab-separated, with a header naming each locus twice. Typings have no locus prefix and may be alleles, groups or NMDP multiple-allele codes:

```text
id A A B B DRB1 DRB1
s1 01:01 24:DXTW 08:01 07:02 03:01 15:01
s2 02:01 02:01 44:02  04:01 07:01
```

`s2` is half-typed at B. Codes such as `01:XX` mean "any allele in family A*01". MAC codes need the MAC table, which is found automatically through `data/manifest.json` or given with `--mac-codes`.

### GL string table (`.glc`, `.glsc`)

One person per line: an ID, then one GL string per locus in any order, separated by whitespace. `|` separates alternative genotypes, `+` the two alleles, `/` alternative alleles:

```text
1 A*01:01+A*02:01/A*02:02 B*08:01+B*07:02 DRB1*13:01+DRB1*04:01|DRB1*13:02+DRB1*04:04
2 A*03:01+ DRB1*15:01+DRB1*15:01
```

A locus left out of a line is untyped (`2` has no B); `A*03:01+` is half-typed.

### Family table (`.csv`, `.tsv`)

One row per person, two columns per locus, plus family and member columns:

```text
FAMILY_ID,Family_Member,A_1,A_2,B_1,B_2,DRB1_1,DRB1_2
fam01,F,A*01:01,A*02:01,B*08:01,B*07:02,DRB1*03:01,DRB1*15:01
fam01,M,03:01,11:01,35:01,44:02,04:01,13:01
fam01,C1,A*01:01,A*03:01,B*08:01,B*35:01,DRB1*03:01,DRB1*04:01
```

- Separators: comma, tab or semicolon; quoted fields are allowed.
- Locus columns: `A_1`, `A.1`, `A1`, `A_ALLELE_1` or `HLA-A_1`.
- Alleles: `A*01:01`, `*01:01` or `01:01`; `/` for alternatives; MAC codes.
- Roles: `F`, `Father`, `Dad` / `M`, `Mother`, `Mum` / `C1`, `C2`, `Child3`, `Sib1`, ... Rows without a family ID are used as unrelated people.
- A `Sample_ID` column, if present, names people; otherwise `family/member` is used (e.g. `Family3/C1`).

## Resolutions

Each locus is resolved to a target resolution, given as `LOCUS:RESOLUTION` (e.g. `--loci A:g,B:g,DRB1:2f`):

| Code | Target | Example for `A*02:01:01:01` |
| --- | --- | --- |
| `1f` | First field | `A*02` |
| `2f` | Two fields | `A*02:01` |
| `3f` | Three fields | `A*02:01:01` |
| `4f` | Four fields | `A*02:01:01:01` |
| `g` | g group (identical antigen-recognition domains, nulls included) | `A*02:01g` |
| `G` | G group (identical exons 2 and 3, or 2 for class II) | `A*02:01:01G` |
| `P` | P group (identical protein in the antigen-recognition domains) | `A*02:01P` |
| `asItIs` | Typings used exactly as written, without lookup | `A*02:01:01:01` |

A typing that covers several alleles at the target resolution (e.g. `A*02` at `2f`, or a MAC code) contributes to each with equal weight. Expression suffixes (N, L, S, Q) are kept. Loci outside the IMGT/HLA tables (e.g. MICA, MICB) can be used with `asItIs`.

## Command line

```bash
uv run hlahaplotools estimate  INPUT --loci A:g,B:g,... [-o results]   # population data
uv run hlahaplotools families  INPUT --loci A:g,B:g,... [-o results]   # family data
uv run hlahaplotools haplomat  MAC|GLSC PARAMETER_FILE                 # run from a parameter file
uv run hlahaplotools summary   INPUT                                   # describe an input file
```

Options for `estimate` and `families`:

| Option | Default | Meaning |
| --- | --- | --- |
| `--loci` | required | Loci and resolutions, e.g. `A:g,C:g,B:g,DRB1:g` |
| `-o`, `--output` | `results` | Output folder |
| `--format` | from extension | `mac` or `glsc` (`estimate` only) |
| `--epsilon` | `1e-6` | Stop when no haplotype frequency changes by more than this |
| `--init` | `numberOccurrence` | Start values: `numberOccurrence`, `perturbation`, `equal`, `random` |
| `--seed` | `0` | Seed for `perturbation` and `random` |
| `--cut` | `1e-6` | Leave out haplotypes at or below this frequency in `hfs.dat` |
| `--no-renormalize` | off | Do not rescale the kept frequencies to sum to 1 |
| `--min-genotype-frequency` | `1e-5` | Leave out a sample whose typing allows more than 1/this genotypes |
| `--max-completions` | `10000` | Most haplotypes created to fill in one partly typed pattern |
| `--threads` | all cores | Worker threads; results do not depend on it |
| `--no-family-phase` | off | `families` only: estimate from founders as unrelated people |
| `--cross-family-duplicates` | off | `families` only: also flag identical people in different families |
| `--data-dir`, `--nomenclature`, `--mac-codes` | from `data/manifest.json` | Reference data locations |

`haplomat` reads a parameter file of `KEY=value` lines (`FILENAME_INPUT`, `FILENAME_HAPLOTYPEFREQUENCIES`, `FILENAME_EPSILON_LOGL`, `LOCI_AND_RESOLUTIONS`, `EPSILON`, `CUT_HAPLOTYPEFREQUENCIES`, `RENORMALIZE_HAPLOTYPEFREQUENCIES`, `INITIALIZATION_HAPLOTYPEFREQUENCIES`, `SEED`, `MINIMAL_FREQUENCY_GENOTYPES`). Loci are then ordered alphabetically.

The exit status is 0 on success and 2 if the EM did not converge.

## Output files

### `estimate`

| File | Contents |
| --- | --- |
| `hfs.dat` | `haplotype<TAB>frequency`, most frequent first. Alleles are joined with `~` in the order given by `--loci`. |
| `epsilon.dat` | One line per iteration: largest frequency change, log-likelihood |
| `run.json` | Tool version, input, loci and resolutions, IMGT/HLA release, MAC table date, all parameters, numbers of samples read, used and incomplete, iterations, log-likelihood, timings, and every unreadable line and left-out sample with the reason |

### `families`

| File | Contents |
| --- | --- |
| `flags.tsv` | Pedigree problems: family, member, relative, kind, loci, explanation (see [Pedigree checks](#pedigree-checks)) |
| `mendel.tsv` | Per child and locus: typing that cannot come from the parents; families that cannot be phased |
| `phases.tsv` | Per family: phasing status, the father's haplotypes `hap_A` and `hap_B` and the mother's `hap_C` and `hap_D` (`LOCUS*?` where unknown), inconsistent and ambiguous loci, notes |
| `children.tsv` | Per child, in child-number order: `paternal` (A or B), `maternal` (C or D) and both together (e.g. `BC`); `?` where the data do not decide |
| `hfs.dat`, `epsilon.dat` | Haplotype frequencies estimated from founders only |
| `diplotypes.tsv` | Most likely haplotype pair per founder, with its posterior probability |
| `comparison.json` | EM on founders as unrelated people compared with the family phasing (meaningful for larger cohorts) |
| `run.json` | As for `estimate`, plus family counts, phasing status counts and flag counts |

## Python API

```python
import hlahaplotools as h

nom = h.Nomenclature("data/derived/3.66.0")

# Population data
ds = h.read_mac("examples/population/population.dat", mac_codes="data/mac/alpha.v3.txt")
ds.sample(0)  # {'A': 'A*02:1192+A*32:136', 'C': 'C*08:150+', ...}  (C half-typed)
ds.missing_counts()  # untyped samples per locus
est = h.estimate(ds, nom, {"A": "g", "C": "g", "B": "g", "DRB1": "g", "DQB1": "g"}, epsilon=1e-6)
est.frequencies(cut=1e-6)[:3]  # [(haplotype, frequency), ...]
est.best_diplotypes()[:3]  # [(sample, haplotype 1, haplotype 2, probability), ...]
est.issues  # samples left out, with reasons
est.write_frequencies("hfs.dat")

# Family data (loci covered by the IMGT/HLA tables)
loci = ["A", "C", "B", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
fam = h.read_family_table("examples/families/validated.csv", loci=loci)
h.audit_families(fam, nom, "2f")  # pedigree problems
h.check_mendelian(fam, nom, "2f")  # per-locus Mendel failures
phases = h.phase_families(fam, nom, "2f")  # one FamilyPhase per family
phases[0].labelled_haplotypes  # {"A": ..., "B": ..., "C": ..., "D": ...}
phases[0].children  # [("25-FS-5/C1", "A", "C"), ("25-FS-5/C2", "A", "D"), ...]
est, phases = h.estimate_families(fam, nom, "2f")
```

`resolutions` may be one code for all loci (`"g"`), a list in locus order, or a dict. `read_glsc` reads GL string tables. Passing `mac_codes` as a path loads only the codes the input uses; `h.MacCodes(path)` loads the whole table for reuse. Type hints are in [`_core.pyi`](src/hlahaplotools/_core.pyi).

## Methods

### Typing resolution

Each typed allele, group or MAC code is expanded to full IMGT/HLA allele names and mapped to the target resolution. When a typing covers several target codes, each receives an equal share of its weight; when a side of a genotype lists several candidates (`02:01/02:02`), each candidate receives an equal share; each alternative of a GL-string genotype list (`|`) receives weight 1/n. A sample is expanded into all combinations of its per-locus genotypes, each weighted by the product of these shares. A sample with more than `1/min_genotype_frequency` combinations is left out and reported.

### EM algorithm

Haplotype frequencies are estimated by maximum likelihood with the EM algorithm. The likelihood of a multi-locus genotype is the sum, over ordered haplotype pairs (h, h′) consistent with it, of f(h)·f(h′). The E step computes the expected number of copies of each haplotype; the M step divides by twice the number of samples. Iteration stops when no frequency changes by more than `epsilon`; frequencies below 1e-14 are set to zero. Start values are either haplotype occurrence counts over all phase combinations (`numberOccurrence`, the default), those counts with up to 10% seeded noise (`perturbation`), equal frequencies, or seeded random frequencies.

### Untyped and half-typed loci

A sample untyped at a locus is consistent with any allele there; a sample typed `a+?` requires `a` on at least one of its two haplotypes. Writing "a on one of the two" as [h = a] + [h′ = a] − [h = a and h′ = a] and multiplying over loci expresses the likelihood as a signed sum of products F(q)·F(q′), where q and q′ are haplotype *patterns* with open loci and F(q) is the total frequency of the haplotypes matching q. This gives exact E and M steps, which reduce to the standard ones for fully typed samples. Candidate haplotypes are those fully specified by the data; a pattern that matches none of them is filled in with every allele observed at its open loci (at most `max_completions` per pattern). The derivation is documented in [estimate.rs](crates/hlahaplotools-core/src/estimate.rs).

### Family phasing

The father's two haplotypes are labelled A and B, the mother's C and D. Children are taken in order of their member number (`C`, `C1`, `C2`, …, `C10`), whatever the row order in the file, and the lowest-numbered typed child defines the labels: it inherits A from the father and C from the mother. Every other child is then described as A or B plus C or D (e.g. `BC`). Inheritance is decided from all loci together, so a parent homozygous at some loci is resolved by the loci where its haplotypes differ. A parent whose two haplotypes are identical at every typed locus passes the same haplotype to every child; its children are then all reported as receiving A (father) or C (mother), with a note. A side stays `?` only when the data cannot decide it, for example when both parents carry the same two haplotypes and a child has one of each.

For each family, every way the children could have inherited one haplotype from each parent is tried (up to the 8 lowest-numbered children), and the inheritance patterns that explain the most loci are kept. At each locus the parental alleles consistent with every kept pattern are collected: a parent's alleles are assigned to haplotypes where all patterns agree, and left unknown (`LOCUS*?`) otherwise. An untyped parent is rebuilt from the children as far as they show it; when neither parent is typed, which rebuilt parent is the father is arbitrary. Loci no pattern explains are reported as inconsistent (typing error or recombination). All comparisons are made at the target resolution.

### Estimation from families

Relatives share haplotypes, so only founders (parents, typed or rebuilt) enter the estimate. Loci the family phased enter with known phase; other loci use the founder's own typing; unknown alleles of a rebuilt parent are treated as untyped. A family with no usable founders contributes its first child.

### Pedigree checks

| Kind | Meaning |
| --- | --- |
| `parent_excluded` | The child shares no allele with a parent at 2 or more loci: not that parent's child, or a sample swap |
| `single_locus_mismatch` | One locus does not fit: usually a typing error, null allele or mutation |
| `trio_inconsistent` | Each parent fits alone, but no allele from each explains the child |
| `possible_duplicate` | A child identical to a parent at 90% or more of loci that also fails to fit the other parent (a copied row), or two identical parents |
| `not_full_sibling` | With a parent untyped, a child carries alleles the other children's parents could not also have passed on. If any of several children could be the odd one out, one family-level flag names them all |
| `recombination` | The child's inheritance from one parent switches once along chromosome 6 (e.g. haplotype A at A and B, haplotype B at DRB1); the interval is given |
| `inheritance_mismatch` | The inheritance switches several times, or the child lacks alleles an untyped parent must have |
| `structure` | Unknown member roles, a second father or mother, typings that cannot be resolved |

Each child is compared with haplotypes phased from the rest of its family, leaving out members already shown not to belong. Parent–child genotype matches also happen naturally (when the other parent carries the same haplotype), so a match is only reported when the child also fails to fit the other parent. Identical genotypes between unrelated people with common haplotypes are not rare, so cross-family duplicates are only reported on request.

### Determinism

Work is split into fixed-size chunks that do not depend on the thread count, and partial sums are combined in a fixed order, so results are bit-identical on any number of threads.

## Validation

All results below are reproducible with the scripts in [scripts/](scripts). Simulated populations draw real allele names from IMGT/HLA 3.66.0 and give haplotypes skewed (Zipf-like) frequencies. Error is the total variation distance between estimated and true haplotype frequencies (0 = identical).

### Incomplete typing

`uv run python scripts/validate_missing.py` — 5,000 people, 5 loci, 300 true haplotypes, 3 replicates (mean ± standard deviation). C and DQB1 are left untyped at the given rate; "drop incomplete" estimates from fully typed people only.

| C and DQB1 untyped | People left if dropped | Error, all people | Error, drop incomplete | Error, before blanking |
| --- | --- | --- | --- | --- |
| 5% | 4,516 | 0.049 ± 0.001 | 0.052 ± 0.001 | 0.049 ± 0.001 |
| 20% | 3,219 | 0.051 ± 0.004 | 0.060 ± 0.004 | 0.050 ± 0.003 |
| 40% | 1,820 | 0.050 ± 0.003 | 0.081 ± 0.002 | 0.046 ± 0.003 |

With 15% of the remaining C and DQB1 typings also half-typed, errors were 0.049 against 0.065 when dropping (10% untyped) and 0.051 against 0.084 (30% untyped). With B, C and DQB1 each 30% untyped: 0.053 against 0.089.

### Estimation from families

`uv run python scripts/validate_family_estimates.py` — 1,000 families of 1–4 children per replicate, 10% of parents untyped, DQB1 20% untyped.

| Replicate | Founders, family phase | Founders, unphased | Everyone as unrelated |
| --- | --- | --- | --- |
| 1 | 0.076 | 0.080 | 0.084 |
| 2 | 0.078 | 0.083 | 0.085 |
| 3 | 0.078 | 0.082 | 0.082 |

Treating every family member as unrelated also counts about 4,300 people for about 2,000 independent founders, which overstates the precision.

### Phasing real pedigrees

`uv run python scripts/validate_pedigrees.py examples/families/*.csv` — the 9 example families, 18 loci at 2-field resolution. Each family is phased with all members as the reference, then again with members left out; called parental alleles are compared with the reference. Family3 is excluded from these figures because of its copied row.

| Scenario | Correct | Wrong | Alleles called |
| --- | --- | --- | --- |
| Both parents and one child (21 trios) | 999 | 0 | 94.7% |
| One parent left out, 2 children (rebuilt parent) | 708 | 0 | 73.9% |
| One parent left out, 3 children (rebuilt parent) | 276 | 0 | 96.2% |

Every recovered allele of a rebuilt parent (1,026) is part of that parent's actual typing. In the parentage test, all 171 pairings of a child with another family's parents fail the Mendel check (2 to 12 inconsistent loci), and 20 of 21 true children pass; the exception is `Family3/C1`. With 2 children, the untransmitted haplotype of a missing parent is often never observed, which explains the lower recovery.

### Pedigree checks

`uv run python scripts/validate_audit.py --seed 1` (and `--seed 2`) — 2,000 families of 2–4 children, 6 loci, 15% of parents untyped, with problems planted in some children.

| Planted problem | Both parents typed | A parent untyped |
| --- | --- | --- |
| Child from another father | 86–93% flagged | 51–67% |
| Child's row copied from the mother | 84–89% | 86–93% |
| Recombinant child | 80–84% | 41–49% |
| One allele changed (typo) | 97–99% | 44–59% |
| Clean children flagged (false alarms) | 1.2–1.4% | 6.7–7.9% |

## Limitations

- Phasing does not model recombination; a recombinant child shows up as inconsistent loci (or a `recombination` flag) and those loci stay unphased. Only the first 8 children of a family are used to phase it.
- An empty DRB3, DRB4 or DRB5 cell is treated as untyped, not as an absent gene.
- Loci outside the IMGT/HLA tables (H, J, MICA, MICB, ...) must be used with `asItIs`, and their alleles must then be written consistently.
- Excel files are not read; save them as CSV. GL-string ID (GLID) files are not read.
- On highly ambiguous data, the EM likelihood can have several nearly equal maxima; results then depend slightly on the start values.
- With small cohorts, haplotype frequencies and the EM-versus-family comparison are not informative.

## Development

```bash
cargo test -p hlahaplotools-core       # Rust unit tests
uv run pytest                           # Python API, command line and data tests
uv run ruff check . && cargo clippy --all-targets
```

Tests that need the reference data are skipped until `scripts/build_nomenclature.py` has been run.

| Path | Contents |
| --- | --- |
| `crates/hlahaplotools-core/` | Rust core: parsers, MAC codes, nomenclature tables, resolution, EM, families, pedigree checks, output |
| `src/lib.rs` | Python bindings (`hlahaplotools._core`) |
| `src/hlahaplotools/` | Python package and command line |
| `scripts/build_nomenclature.py` | Downloads the pinned reference data and builds the lookup tables |
| `scripts/derive_tables.py` | Builds the lookup tables from an IMGT/HLA release |
| `scripts/simulate.py` | Simulates populations or families (`--format families`, optionally with planted problems) |
| `scripts/validate_*.py` | The validations above |
| `examples/` | Example population and family data |
| `tests/` | Test suite |

## Citation

If you use HLAhaploTools, please cite this repository, and report the IPD-IMGT/HLA release and MAC table date recorded in `run.json`.

## License and acknowledgements

HLAhaploTools is released under the GNU General Public License, version 3 or later ([LICENSE](LICENSE)).

Some functions and scripts are ported from [Hapl-o-Mat](https://github.com/DKMS/Hapl-o-Mat) by DKMS gGmbH, under the GPL-3.0 license: typing resolution and the lookup-table build, the EM algorithm, and the MAC, GL-string and parameter-file input formats.

Nomenclature data come from the [IPD-IMGT/HLA Database](https://www.ebi.ac.uk/ipd/imgt/hla/), and multiple-allele codes from the [NMDP](https://hml.nmdp.org/mac/).

This repository previously held an R package of the same name; its last version is tagged [`r-package-v0.99.0`](https://github.com/fmobegi/HLAhaploTools/tree/r-package-v0.99.0).
