# Getting started

The user guide of HLAhaploTools: how to run it, what it reads and writes, and how it works. Install it and download the reference data first ([README](../README.md#installation)). From a clone, put `uv run` before each command. For scripting, see the [Python API](python-api.md).

## Contents

- [Quick start](#quick-start)
- [Input formats](#input-formats)
- [Resolutions](#resolutions)
- [Command line](#command-line)
- [Output files](#output-files)
- [Good practice](#good-practice)
- [How it works](#how-it-works)

## Quick start

### Population data

```bash
hlahaplotools estimate INPUT --loci A:g,C:g,B:g,DRB1:g,DQB1:g -o results/population
```

`INPUT` is a typing file in any of the [input formats](#input-formats). The screen reports how many people were read, what was filled in at each locus, how many people are untyped or half-typed at each locus, whether the EM starts agreed, and where the results were written:

```text
results/population/
├── hfs.dat       haplotype frequencies, most frequent first
├── epsilon.dat   EM convergence, one line per iteration
└── run.json      what was run: input, loci, reference data versions, parameters, people used and left out
```

Each line of `hfs.dat` is one haplotype, its alleles joined by `~` in the order of `--loci`, then a tab and its frequency:

```text
A*01:01~C*07:01~B*08:01~DRB1*03:01~DQB1*02:01 0.0712
```

### Family data

```bash
hlahaplotools families FAMILY_TABLE --loci A:2f,C:2f,B:2f,DRB1:2f,DQB1:2f -o results/families
```

`FAMILY_TABLE` has one row per person, with family and member columns (see [Family table](#family-table-csv-tsv)). The results:

```text
results/families/
├── flags.tsv        possible pedigree problems (non-paternity, sample swaps, copied rows, recombinants)
├── mendel.tsv       typings a child cannot have inherited from its parents
├── phases.tsv       each family's four parental haplotypes (A, B from the father; C, D from the mother)
├── children.tsv     which haplotypes each child inherited
├── hfs.dat          haplotype frequencies estimated from the founders (parents) only
├── diplotypes.tsv   most likely haplotype pair of each founder
├── comparison.json  population-style estimate compared with the family phasing
└── run.json         as for estimate
```

Check `flags.tsv` first: correct or remove flagged people and rerun. Every file is described under [Output files](#output-files). Example data are in [examples/](../examples).

## Input formats

| Format          | Files                  | Notes                                                                                                                                                                                                                                                |
| --------------- | ---------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| MAC table       | `.dat`, `.mac`     | Tab-separated, a header naming each locus twice; typings without locus prefix (`01:01`, `24:DXTW`, `01:XX`). MAC codes need the NMDP table (found through `data/manifest.json`, or `--mac-codes`).                                         |
| GL string table | `.glc`, `.glsc`    | One person per line: an ID, then a [GL String](https://pmc.ncbi.nlm.nih.gov/articles/PMC10524834/) (`^` loci, `\|` alternative genotypes, `+` the two copies, `/` alternative alleles). Hapl-o-Mat's whitespace-separated layout is also read. |
| GLS             | `.pull` + `.glid`  | GL identifier files: the GLID file numbers each locus GL string, the pull file lists one GL ID per locus per person (`0` = untyped).                                                                                                               |
| HML             | `.hml`, or a folder  | NMDP HML 1.0 and 1.0.1, from `<glstring>` or `<haploid>`. GL strings given only by URI are left untyped and listed.                                                                                                                              |
| NGSengine       | `.xml`, or a folder  | GenDx project export; each locus's `ConclusiveTypingResult`. Results not approved in review are left untyped unless `--include-unapproved`.                                                                                                      |
| Genotype list   | `--format genotypes` | Resolved genotypes written by `estimate --write-genotypes` (Hapl-o-Mat's `WRITE_GENOTYPES`); estimate with `asItIs`.                                                                                                                           |
| Family table    | `.csv`, `.tsv`     | One row per person, family and member columns, two columns per locus (below). Used by `families`.                                                                                                                                                  |

The format is taken from the extension (and, for `.xml`/`.hml`, the content); `--format` overrides it. A few rules apply to every format:

- **Missing typings** (empty, `-`, `.`, `NA`, `N/A`, `?`) are allowed. A locus with neither allele typed is *untyped* and is used as such in the estimate.
- **A blank second allele** is read as homozygous, the usual reporting convention; `--half-typed unknown` keeps it unknown instead (*half-typed*).
- **DRB3, DRB4 and DRB5** are filled from the DRB1 typing: a gene copy that the DRB1 group cannot carry becomes the absent-gene code `DRBx*NNNN`, and a single allele is read as one or two copies as DRB1 implies. `--no-drb345` turns this off. Every fill is counted per locus on screen and in `run.json`.
- **Ambiguity** is kept in full: allele lists, genotype lists, MAC codes and groups are weighed by the EM, never resolved by hand. Phased GL strings (`~`) are read without their phase.
- **Serological antigens and pre-2010 names** (`B*62`, `A*0201`) leave the person out by default, with a message; `--serology` and `--v2-names` read them.
- **Conflicts**: if the same person and locus appear twice with different results, the locus is left untyped and listed.

### Family table (`.csv`, `.tsv`)

```text
FAMILY_ID,Family_Member,A_1,A_2,B_1,B_2,DRB1_1,DRB1_2
fam01,F,A*01:01,A*02:01,B*08:01,B*07:02,DRB1*03:01,DRB1*15:01
fam01,M,03:01,11:01,35:01,44:02,04:01,13:01
fam01,C1,A*01:01,A*03:01,B*08:01,B*35:01,DRB1*03:01,DRB1*04:01
```

Separators: comma, tab or semicolon. Locus columns: `A_1`, `A.1`, `A1`, `A_ALLELE_1` or `HLA-A_1`. Alleles: `A*01:01`, `*01:01` or `01:01`, `/` for alternatives, MAC codes. Roles: `F`/`Father`/`Dad`, `M`/`Mother`/`Mum`, `C1`, `C2`, `Child3`, `Sib1`, ... Rows without a family ID are used as unrelated people. A `Sample_ID` column, if present, names people; otherwise `family/member` (e.g. `fam01/C1`).

## Resolutions

Each locus is resolved to a target resolution, given as `LOCUS:RESOLUTION` (e.g. `--loci A:g,B:g,DRB1:2f`):

| Code | Target | Example for `A*02:01:01:01` |
| ---------- | -------------------------------------------------------------------------------------------------- | ----------------------------- |
| `1f` | First field | `A*02` |
| `2f` | Two fields | `A*02:01` |
| `3f` | Three fields | `A*02:01:01` |
| `4f` | Four fields | `A*02:01:01:01` |
| `g` | g group (identical antigen-recognition domains, nulls included) | `A*02:01g` |
| `G` | G group (identical exons 2 and 3, or 2 for class II) | `A*02:01:01G` |
| `P` | P group (identical protein in the antigen-recognition domains) | `A*02:01P` |
| `lgx` | g group named by its two fields, without the `g` (py-ard's `lgx`; `lg` is the same as `g`) | `A*02:01` |
| `asItIs` | Typings used exactly as written, without lookup | `A*02:01:01:01` |

A typing that covers several alleles at the target resolution (e.g. `A*02` at `2f`, or a MAC code) contributes to each with equal weight. Expression suffixes (N, L, S, Q) are kept. Every locus in the IPD-IMGT/HLA allele list can be resolved, including MICA, MICB and the non-classical and pseudogene loci (E, F, G, H, J, K, L, DQA2, TAP1, ...). IMGT/HLA defines G and P groups only for the classical and E/F/G-type loci; for the others, use `1f` to `4f` or `g` (the 2-field name). Loci not in the list at all can still be used with `asItIs`.

## Command line

```bash
hlahaplotools estimate  INPUT --loci A:g,B:g,... [-o results]          # population data
hlahaplotools families  INPUT --loci A:g,B:g,... [-o results]          # family data
hlahaplotools stats     INPUT --loci A:2f,B:2f,... [-o results/stats]  # Hardy–Weinberg and LD
hlahaplotools impute    INPUT --frequencies hfs.dat [-o results/impute] # each person's haplotype pairs
hlahaplotools haplomat  MAC|GLSC|GLS|READ PARAMETER_FILE               # run from a Hapl-o-Mat parameter file
hlahaplotools summary   INPUT                                          # describe an input file
hlahaplotools rename    INPUT OUTPUT                                   # update retired allele names
hlahaplotools download-data [--skip-mac]                               # reference data, once
```

`hlahaplotools COMMAND --help` lists every option of a command. Options of `estimate` and `families` (most are shared by `stats` and `impute`):

| Option                                              | Default                     | Meaning                                                                                                                               |
| --------------------------------------------------- | --------------------------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| `--loci`                                          | required                    | Loci and resolutions, e.g. `A:g,C:g,B:g,DRB1:g`                                                                                     |
| `-o`, `--output`                                | `results`                 | Output folder                                                                                                                         |
| `--format`                                        | from extension or content   | `mac`, `glsc`, `gls`, `hml`, `ngsengine` or `genotypes` (`estimate`, `stats`, `impute`, `summary`)                |
| `--glid`                                          | pull file with `.glid`    | GLS input: the GLID file                                                                                                              |
| `--include-unapproved`                            | off                         | NGSengine input: also read locus results not approved in review                                                                       |
| `--epsilon`                                       | `1e-6`                    | Stop when no haplotype frequency changes by more than this                                                                            |
| `--init`                                          | `numberOccurrence`        | Start values: `numberOccurrence`, `perturbation`, `equal`, `random`                                                           |
| `--seed`                                          | `0`                       | Seed for `perturbation` and `random`                                                                                              |
| `--starts`                                        | `10`                      | EM start values: `--init` first, then perturbed starts; the best fit is kept (`estimate`, `stats`)                              |
| `--cut`                                           | `1e-6`                    | Leave out haplotypes at or below this frequency in `hfs.dat`                                                                        |
| `--no-renormalize`                                | off                         | Do not rescale the kept frequencies to sum to 1                                                                                       |
| `--min-genotype-frequency`                        | `1e-5`                    | Leave out a sample whose typing allows more than 1/this genotypes                                                                     |
| `--max-completions`                               | `10000`                   | Most haplotypes created to fill in one partly typed pattern                                                                           |
| `--completion`                                    | `linked`                  | How untyped loci are filled in: `linked` (alleles that go with the most strongly linked typed locus) or `all` (every combination) |
| `--threads`                                       | all cores                   | Worker threads; results do not depend on it                                                                                           |
| `--half-typed`                                    | `homozygous`              | One allele typed, the other blank: `homozygous` (copy it) or `unknown`                                                            |
| `--no-drb345`                                     | off                         | Do not fill DRB3/4/5 from DRB1                                                                                                        |
| `--serology`                                      | off                         | Read typings that are serological antigens (`B*62`) as their alleles                                                                |
| `--v2-names`                                      | off                         | Read pre-2010 allele names (`A*0201`) as current names                                                                              |
| `--write-genotypes`                               | off                         | `estimate` only: also write each person's resolved genotypes to `genotypes.dat`                                                   |
| `--no-family-phase`                               | off                         | `families` only: estimate from founders as unrelated people                                                                         |
| `--cross-family-duplicates`                       | off                         | `families` only: also flag identical people in different families                                                                   |
| `--data-dir`, `--nomenclature`, `--mac-codes` | from `data/manifest.json` | Reference data locations                                                                                                              |
| `-q`, `--quiet`                                 | off                         | Less screen output                                                                                                                    |

`stats` adds `--permutations` (random tables per Hardy–Weinberg test, default 10,000), `--pairs all|adjacent` (every pair of loci, or only neighbours in `--loci` order) and `--ld-permutations` (shuffled datasets per linkage-disequilibrium test; default 0, no P-value). A permutation P cannot be smaller than 1/(N + 1), so choose N large enough for the threshold after multiple-testing correction.

`impute` gives each person's most likely haplotype pairs and multilocus genotypes from a frequency table (`hfs.dat` from `estimate`, or any `haplotype<TAB>frequency` file), as GRIMM and py-graph-imputation do. Loci and resolutions are taken from the `run.json` beside the frequency file unless `--loci` is given, so the typing is resolved as it was for the estimate. Repeat `--frequencies NAME=FILE` for several populations (with `--prior NAME=WEIGHT`, default equal) to also get each person's population posteriors. `--top` (default 10) sets how many pairs and genotypes are written per person; `--max-candidates` (default 5,000) caps the haplotypes considered per person, keeping the most frequent. Only haplotypes in the table can be proposed: impute from an uncut table (`estimate --cut 0`) to keep people who carry a rare haplotype.

`haplomat` runs from an existing Hapl-o-Mat parameter file, with a single EM start, so earlier analyses can be repeated unchanged. `WRITE_GENOTYPES=true` writes the genotype list to `FILENAME_GENOTYPES`, and `READ` estimates from such a list. The ambiguity filter (`DO_AMBIGUITYFILTER`) and the haplotype analytics (`DO_ANALYTICS`) are not implemented, and are ignored with a warning: ambiguity is kept in full and weighed by the EM instead.

`rename` rewrites a typing file (MAC, GLSC or family table) with current IMGT/HLA names. Alleles deleted or renamed since the typing was reported (e.g. `C*04:09N`, now `C*04:09L`, or `B*47:01:01:01`, now `B*47:01:01:03`) are replaced with the name `Deleted_alleles.txt` gives for them, following chains of renames, and everything else is written unchanged. The changes are listed in `OUTPUT.renames.tsv`.

`summary` describes an input file without estimating (`--loci` keeps only the loci listed).

The exit status is 0 on success and 2 if the EM did not converge.

## Output files

### `estimate`

| File              | Contents                                                                                                                                                                                                                                                                                                                                                                                  |
| ----------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `hfs.dat`       | `haplotype<TAB>frequency`, most frequent first. Alleles are joined with `~` in the order given by `--loci`.                                                                                                                                                                                                                                                                         |
| `epsilon.dat`   | One line per iteration: largest frequency change, log-likelihood                                                                                                                                                                                                                                                                                                                          |
| `genotypes.dat` | With `--write-genotypes`: each person's resolved genotypes (see [Input formats](#input-formats))                                                                                                                                                                                                                                                                                         |
| `run.json`      | Tool version, input, loci and resolutions, IMGT/HLA release, MAC table date, all parameters, how blank alleles were filled (rules and counts per locus), numbers of samples read, used and incomplete, iterations, log-likelihood of every start, timings, and every unreadable line, locus result left untyped while reading (`loci_left_untyped`) and left-out sample with the reason |

### `families`

| File                         | Contents                                                                                                                                                                                 |
| ---------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `flags.tsv`                | Pedigree problems: family, member, relative, kind, loci, explanation (see [Pedigree checks](#pedigree-checks))                                                                            |
| `mendel.tsv`               | Per child and locus: typing that cannot come from the parents; families that cannot be phased                                                                                            |
| `phases.tsv`               | Per family: phasing status, the father's haplotypes `hap_A` and `hap_B` and the mother's `hap_C` and `hap_D` (`LOCUS*?` where unknown), inconsistent and ambiguous loci, notes |
| `children.tsv`             | Per child, in child-number order: `paternal` (A or B), `maternal` (C or D) and both together (e.g. `BC`); `?` where the data do not decide                                       |
| `hfs.dat`, `epsilon.dat` | Haplotype frequencies estimated from founders only                                                                                                                                       |
| `diplotypes.tsv`           | Most likely haplotype pair per founder, with its posterior probability                                                                                                                   |
| `comparison.json`          | EM on founders as unrelated people compared with the family phasing (informative only for larger cohorts)                                                                                |
| `run.json`                 | As for `estimate`, plus family counts, phasing status counts and flag counts                                                                                                           |

### `impute`

| File                | Contents                                                                                                                                                                |
| ------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `pairs.tsv`       | Per person: rank, population, the two haplotypes and the probability of the pair, most likely first                                                                     |
| `genotypes.tsv`   | Per person: rank, multilocus unphased genotype (`A*01:01+A*02:01^B*...`) and its probability, summed over the pairs that give it                                      |
| `populations.tsv` | With several populations: each person's posterior for each                                                                                                              |
| `run.json`        | Inputs, frequency files and priors, loci, parameters, people imputed, and the people left out with the reason (unreadable typing, or no pair in the table that fits it) |

### `stats`

| File               | Contents                                                                                                                                                                                                              |
| ------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `hwe.tsv`        | Per locus: people tested, alleles, observed and expected heterozygosity, F, exact-test P, its false-discovery-rate q-value, and people not tested (untyped, half-typed, ambiguous) or read as homozygous from a blank |
| `ld.tsv`         | Per pair of loci: people used, alleles at each locus, D′, Wn, likelihood-ratio statistic, nominal degrees of freedom, and the permutation P-value and its q-value (with `--ld-permutations`)                       |
| `ld_alleles.tsv` | Per pair of loci and allele pair: haplotype frequency, allele frequencies, D and D′, most frequent first                                                                                                             |
| `run.json`       | Input, loci and resolutions, reference data versions and parameters                                                                                                                                                   |

## Good practice

1. **Look at the input first** with `hlahaplotools summary INPUT`: people read, unreadable lines, untyped and half-typed people per locus, and what was filled in.
2. **Update retired allele names** with `hlahaplotools rename INPUT OUTPUT` (e.g. `C*04:09N`, now `C*04:09L`); otherwise people with a retired name are left out.
3. **Choose the loci and resolution.** `2f` suits most work; `g` or `G` groups mixed typing depths. Every extra locus multiplies the possible haplotypes, and a typing coarser than the target is spread over all the alleles it covers, so start with the loci you need.
4. **Check the reporting conventions**: add `--half-typed unknown` if a blank second allele means "not typed" in your laboratory.
5. **Keep the defaults** (ten EM starts, linked completion); the screen says whether the starts agreed.
6. **Run `stats`** to check Hardy–Weinberg proportions; a locus far from them often points to typing problems or a mixed population.
7. **Report** the IPD-IMGT/HLA release and MAC table date from `run.json`, with the loci, resolution and number of people used. Treat rare haplotypes (below about 1 in 1,000) as uncertain.

Memory grows with the number of candidate haplotypes, which grows with the loci, the typing depth and the people untyped at several loci: five classical loci take seconds and well under 1 GB, eleven loci with incomplete typing can take minutes and several GB. `--threads` changes the speed, never the result.

## How it works

### Typing resolution

Each typed allele, group or MAC code is expanded to full IMGT/HLA allele names and mapped to the target resolution. When a typing covers several target codes, each receives an equal share of its weight; when a side of a genotype lists several candidates (`02:01/02:02`), each candidate receives an equal share; each alternative of a GL-string genotype list (`|`) receives weight 1/n. A sample is expanded into all combinations of its per-locus genotypes, each weighted by the product of these shares. A sample with more than `1/min_genotype_frequency` combinations is left out and reported.

### EM algorithm and starts

Haplotype frequencies are estimated by maximum likelihood with the EM algorithm. The likelihood of a multi-locus genotype is the sum, over ordered haplotype pairs (h, h′) consistent with it, of f(h)·f(h′). Iteration stops when no frequency changes by more than `epsilon`; frequencies below 1e-14 are set to zero. Start values are haplotype occurrence counts over all phase combinations (`numberOccurrence`, the default), those counts with up to 10% seeded noise (`perturbation`), equal frequencies, or seeded random frequencies.

The EM converges to a stationary point near its start, which need not be the maximum. `estimate` therefore runs 10 starts by default, all on one model built once: the `--init` start, then 9 perturbed starts (seeds `seed + 1` to `seed + 9`), and keeps the highest log-likelihood. Every start's log-likelihood is recorded in `run.json`, and the screen output says whether the starts disagreed. `--starts 1` reproduces a single Hapl-o-Mat run.

### Blank alleles

Blanks are filled when the file is read, before resolution, with the rules in [Input formats](#input-formats); `DRBx*NNNN` is an allele like any other and stays `DRBx*NNNN` at every resolution. In unrelated people the homozygous reading cannot be checked. In families it can: after resolution, each family and locus where no inheritance pattern fits is checked again with copied alleles read as unknown, and the copies of the fewest members that make the locus fit are undone (noted in `phases.tsv`). The rules are in [fill.rs](../crates/hlahaplotools-core/src/fill.rs).

### Untyped and half-typed loci

A sample untyped at a locus is consistent with any allele there; a sample typed `a+?` requires `a` on at least one of its two haplotypes. The likelihood is written as a signed sum over haplotype *patterns* with open loci, which gives exact E and M steps that reduce to the standard ones for fully typed samples. Candidate haplotypes are those fully specified by the data; a pattern that matches none of them is filled in at its open loci. With `--completion linked` (the default) an open locus takes only the alleles that go with the pattern's allele at its most strongly linked typed locus: the known linked pairs first (DQA1 with DQB1, DPA1 with DPB1, DRB3/4/5 with DRB1, C with B, MICA and MICB with B), then the nearest typed locus along the MHC, with the allele associations estimated by a two-locus EM. `--completion all` fills in every combination of the alleles seen (at most `max_completions` per pattern), which can create millions of candidates when several loci are untyped together. The derivation is documented in [estimate.rs](../crates/hlahaplotools-core/src/estimate.rs).

### Family phasing

The father's two haplotypes are labelled A and B, the mother's C and D. Children are taken in order of their member number (`C`, `C1`, `C2`, …, `C10`), whatever the row order in the file, and the lowest-numbered typed child defines the labels: it inherits A from the father and C from the mother. Every other child is then described as A or B plus C or D (e.g. `BC`). Inheritance is decided from all loci together, so a parent homozygous at some loci is resolved by the loci where its haplotypes differ. A side stays `?` only when the data cannot decide it.

For each family, every way the children could have inherited one haplotype from each parent is tried (up to the 8 lowest-numbered children), and the inheritance patterns that explain the most loci are kept. A parent's alleles are assigned to haplotypes where all kept patterns agree, and left unknown (`LOCUS*?`) otherwise. An untyped parent is rebuilt from the children as far as they show it; when neither parent is typed, which rebuilt parent is the father is arbitrary. Loci no pattern explains are reported as inconsistent (typing error or recombination). All comparisons are made at the target resolution.

### Estimation from families

Relatives share haplotypes, so only founders (parents, typed or rebuilt) enter the estimate. Loci the family phased enter with known phase; other loci use the founder's own typing; unknown alleles of a rebuilt parent are treated as untyped. A family with no usable founders contributes its first child. Treating every family member as an unrelated person would count shared haplotypes more than once and overstate the precision.

### Pedigree checks

| Kind                      | Meaning                                                                                                                                                                                              |
| ------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `parent_excluded`       | The child shares no allele with a parent at 2 or more loci: not that parent's child, or a sample swap                                                                                                |
| `single_locus_mismatch` | One locus does not fit: usually a typing error, null allele or mutation                                                                                                                              |
| `trio_inconsistent`     | Each parent fits alone, but no allele from each explains the child                                                                                                                                   |
| `possible_duplicate`    | A child identical to a parent at 90% or more of loci that also fails to fit the other parent (a copied row), or two identical parents                                                                |
| `not_full_sibling`      | With a parent untyped, a child carries alleles the other children's parents could not also have passed on. If any of several children could be the odd one out, one family-level flag names them all |
| `recombination`         | The child's inheritance from one parent switches once along chromosome 6 (e.g. haplotype A at A and B, haplotype B at DRB1); the interval is given                                                   |
| `inheritance_mismatch`  | The inheritance switches several times, or the child lacks alleles an untyped parent must have                                                                                                       |
| `structure`             | Unknown member roles, a second father or mother, typings that cannot be resolved                                                                                                                     |

Each child is compared with haplotypes phased from the rest of its family, leaving out members already shown not to belong. A parent–child genotype match is only reported when the child also fails to fit the other parent. Identical genotypes between unrelated people with common haplotypes are not rare, so cross-family duplicates are only reported on request.

### Hardy–Weinberg and linkage disequilibrium

**Hardy–Weinberg.** Only people whose genotype at the locus is unambiguous at the target resolution are tested; untyped, half-typed and ambiguous people are counted apart. Expected heterozygosity is Nei's unbiased estimate, 2n/(2n−1)(1 − Σp²), and F = 1 − Ho/He (positive when heterozygotes are too few). The exact test compares the probability of the observed genotype table, given the allele counts, with that of tables formed by pairing the same alleles at random. With large samples even small departures are significant, so read F alongside P.

**Linkage disequilibrium.** For each pair of loci, haplotype frequencies are estimated by EM from everyone typed at either locus. D′ for an allele pair divides D by its largest possible value with the same sign; the overall D′ (Hedrick) is Σ p_a q_b |D′_ab|; and Wn (Cohen) is Cramér's V of the haplotype table. With `--ld-permutations N`, linkage equilibrium is tested by permutation: one locus is shuffled across people and the two-locus likelihood re-estimated. The likelihood-ratio statistic is also written, with its nominal degrees of freedom but no χ² P-value, because the χ² approximation fails with many rare alleles.

**Multiple testing.** `hwe.tsv` and `ld.tsv` give `q_value`, the Benjamini–Hochberg false-discovery-rate adjusted P across the loci or pairs of the run; report it rather than the raw P.

### Imputation

For a person with possible genotypes g (with weights w, from their typing), every pair of haplotypes h₁, h₂ in the frequency table gets P(h₁, h₂) ∝ f(h₁) f(h₂) Σ_orders Π_loci Σ_g w · [h₁, h₂ fit g], normalised over the person's pairs. Untyped and half-typed loci constrain only the typed allele, so the pairs impute the rest. With several populations, both haplotypes of a pair come from one population, weighted by its prior, and the population posteriors are the shares of the pairs from each.

### Determinism

Work is split into fixed-size chunks that do not depend on the thread count, and partial sums are combined in a fixed order, so results are bit-identical on any number of threads.
