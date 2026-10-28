# HLAhaploTools

HLAhaploTools estimates HLA haplotype frequencies from population and family typing data with an expectation–maximisation (EM) algorithm.

> [!CAUTION]
> **Research use only.** HLAhaploTools is research software. It is not a medical device and has not been validated for clinical use: do not use it, or its haplotype assignments, pedigree checks or frequencies, to make decisions about patient care, donor selection or transplantation. Results are statistical estimates and can be wrong. It is provided without warranty of any kind (see [LICENSE](LICENSE)).

## Key strengths

- **Laboratory reports as input**: MAC and GL string tables, GL identifier (GLS) pull and GLID files, NMDP HML files, and GenDx NGSengine XML exports, read with all their ambiguity.
- **Ambiguous typings**: NMDP multiple-allele codes (`24:DXTW`, `01:XX`), allele lists (`02:01/02:02`), GL strings (phased ones too), and g, G and P groups, resolved against a pinned IPD-IMGT/HLA release. Serological antigens (`B62`) and pre-2010 allele names (`A*0201`) can be read as well.
- **Any target resolution**: per locus: first to fourth field, g, G or P groups, or alleles exactly as typed.
- **Incomplete typing**: samples with an untyped locus, or with only one allele typed at a locus, are kept in the estimate instead of discarded.
- **HLA reporting conventions**: a blank second allele is read as homozygous, and blank DRB3, DRB4 and DRB5 cells are filled from the DRB1 typing (absent genes included). In families, a homozygous reading that breaks Mendelian inheritance is undone.
- **Population statistics**: Hardy–Weinberg exact tests per locus, and linkage disequilibrium (D′, Wn and a permutation test) per pair of loci, estimated from incompletely typed people too, with false-discovery-rate adjusted P-values.
- **Imputation**: each person's most likely haplotype pairs and genotypes from known haplotype frequencies (one or several populations), including the loci they were not typed at.
- **Family data**: Mendel checks, phasing of parents from their children (including untyped parents), estimation from founders only, and checks that flag non-paternity, sample swaps, copied rows, non-siblings and recombinants.
- **Reproducible**: results are bit-identical on any number of threads, and every run records the IMGT/HLA release, MAC table date and all parameters.

## Contents

- [Installation](#installation)
- [Reference data](#reference-data)
- [Quick start](#quick-start)
- [Recommended workflow](#recommended-workflow)
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
- [Use of AI](#use-of-ai)

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
| ---------------------- | -------------------------------------------------------------------------------------------------------- |
| `imgthla/<release>/` | `hla_nom_g.txt`, `hla_nom_p.txt`, `hla_nom.txt`, `Allelelist.txt`, `Deleted_alleles.txt` |
| `derived/<release>/` | Lookup tables used for resolution: `AllAllelesExpanded.txt`, `LargeG.txt`, `P.txt`, `Smallg.txt` |
| `mac/alpha.v3.txt` | NMDP multiple-allele codes |
| `manifest.json` | Release, source URLs, SHA-256 checksums, MAC table date |

The IMGT/HLA release is pinned (`PINNED_RELEASE` in [scripts/build_nomenclature.py](scripts/build_nomenclature.py), currently 3.66.0). To use another release, run with `--release 3.65.0`. The MAC table has no release number: it is downloaded once and kept, and its `LAST UPDATED` date is recorded; `--refresh-mac` downloads a new copy. Report both with any results.

Commands find `data/` in the current folder, in `$HLAHAPLOTOOLS_DATA`, or with `--data-dir`.

## Quick start

### Population data

```bash
uv run hlahaplotools estimate INPUT --loci A:g,C:g,B:g,DRB1:g,DQB1:g -o results/population
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
uv run hlahaplotools families FAMILY_TABLE --loci A:2f,C:2f,B:2f,DRB1:2f,DQB1:2f -o results/families
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

Check `flags.tsv` first: correct or remove flagged people and rerun. Every file is described under [Output files](#output-files).

## Recommended workflow

1. **Look at the input first.** `hlahaplotools summary INPUT` reads the file without estimating anything: how many people were read, which lines could not be read and why, how many people are untyped or half-typed at each locus, and what the reporting conventions filled in. Fix unreadable lines before going further; anyone left out of the estimate is listed in `run.json`.
2. **Update retired allele names.** Typings reported years ago may use names IPD-IMGT/HLA has since changed (e.g. `C*04:09N`, now `C*04:09L`). `hlahaplotools rename INPUT OUTPUT` writes a corrected copy and a table of every change; without it, people with a retired name are left out.
3. **Choose the loci and resolution.** Two-field resolution (`2f`) suits most work; use `g` or `G` when typings mix depths and you want them grouped. A typing coarser than the target (e.g. `A*02` at `2f`) is spread over all the alleles it covers, which costs time and adds uncertainty. Every extra locus multiplies the number of possible haplotypes, so start with the loci you need.
4. **Check the reporting conventions.** By default a blank second allele means homozygous and blank DRB3/4/5 cells are filled from DRB1. If your laboratory leaves a second allele blank because it was not typed, add `--half-typed unknown`.
5. **Run with the defaults.** Ten EM starts and linked completion of untyped loci are the defaults for good reason; the screen says whether the starts agreed. If they did not, the best fit is kept and the spread is recorded in `run.json`.
6. **Check Hardy–Weinberg and linkage disequilibrium.** `hlahaplotools stats INPUT --loci ...` tests each locus for Hardy–Weinberg proportions and measures linkage disequilibrium between loci. A locus far from Hardy–Weinberg often points to typing problems (missed heterozygotes, or blanks read as homozygous) or to a mixed population.
7. **Report the run.** Quote the IPD-IMGT/HLA release and MAC table date from `run.json`, together with the loci, resolution and number of people used. Rare haplotypes, below about 1 in 1,000, are uncertain at any cohort size typical of HLA studies (see [Validation](#validation)).

**Run time and memory.** Measured on a laptop (4 cores, 16 GB) for a cohort of 15,633 people typed at mixed depth with DQA1 and DPA1 untyped in about 30%, at two-field resolution with ten starts:

| Loci                    | People incomplete | Candidate haplotypes | Time    | Peak memory |
| ----------------------- | ----------------- | -------------------- | ------- | ----------- |
| 5 (A, B, C, DRB1, DQB1) | 217               | 0.26 million         | 13 s    | 0.3 GB      |
| 6 (+ DPB1)              | 403               | 0.7 million          | 46 s    | 0.6 GB      |
| 8 (+ DQA1, DPA1)        | 4,760             | 4.7 million          | 2.5 min | 1.8 GB      |
| 11 (+ DRB3, DRB4, DRB5) | 4,831             | 16 million           | 8 min   | 6.1 GB      |

Memory grows with the number of candidate haplotypes, which grows with the loci and with how many people are untyped at several of them. Family analyses of a few hundred families take seconds to minutes.

## Input formats

Missing typings are always allowed: empty cells, or `-`, `.`, `NA`, `N/A`, `?`. A locus with neither allele typed is *untyped* and is used in the estimate as such (see [Methods](#methods)). A locus with one allele typed and the other blank is read following HLA reporting practice:

- **Classical loci** (A, B, C, DRB1, DQA1, DQB1, DPA1, DPB1, and any other locus except DRB3/4/5): the person is homozygous, and the typed allele is copied across. With `--half-typed unknown` the blank is kept as an unknown allele instead (*half-typed*, also used in the estimate).
- **DRB3, DRB4, DRB5**: these genes are carried with particular DRB1 groups (DR52: DRB1\*03, 11, 12, 13, 14 carry DRB3; DR53: DRB1\*04, 07, 09 carry DRB4; DR51: DRB1\*15, 16 carry DRB5; DRB1\*01, 08 and 10 carry none). From the person's DRB1 typing, each gene copy that cannot be present is set to the absent-gene code `DRBx*NNNN`, a single typed allele where two copies are expected is read as homozygous, the same allele written twice where only one copy is expected is read as one copy (a common way of reporting a single copy), and an expected copy that was not typed stays unknown. `NNNN` may also be written in the input. Nothing is filled when DRB1 is untyped or its candidates span several groups. Typings with more copies than DRB1 allows are left as typed and counted as conflicts. `--no-drb345` turns these rules off; a blank DRB3/4/5 cell is then never copied.

Every fill is counted per locus on screen and in `run.json`.

### MAC table (`.dat`, `.mac`)

Tab-separated, with a header naming each locus twice. Typings have no locus prefix and may be alleles, groups or NMDP multiple-allele codes:

```text
id A A B B DRB1 DRB1
s1 01:01 24:DXTW 08:01 07:02 03:01 15:01
s2 02:01 02:01 44:02  04:01 07:01
```

`s2` has one B allele typed, read as `B*44:02+B*44:02`. Codes such as `01:XX` mean "any allele in family A*01". MAC codes need the MAC table, which is found automatically through `data/manifest.json` or given with `--mac-codes`.

### GL string table (`.glc`, `.glsc`)

One person per line: an ID, then the person's [GL String](https://pmc.ncbi.nlm.nih.gov/articles/PMC10524834/) (Genotype List String 1.1). In the grammar, `^` separates loci, `|` alternative genotypes, `+` the two copies at a locus, and `/` alternative alleles:

```text
1 A*01:01+A*02:01/A*02:02^B*08:01+B*07:02^DRB1*13:01+DRB1*04:01|DRB1*13:02+DRB1*04:04
2 A*03:01+A*03:01^DRB1*15:01+DRB1*15:01
```

Person 1 has `A*01:01` and either `A*02:01` or `A*02:02` (`/` binds more tightly than `+`), and at DRB1 either `13:01+04:01` or `13:02+04:04` (two alternative genotypes, `|`). A locus left out is untyped (`2` has no B). The `HLA-` prefix is optional, and the loci may also be separated by whitespace instead of `^` (the layout of Hapl-o-Mat's GLSC files), e.g. `2 A*03:01+A*03:01 DRB1*15:01+DRB1*15:01`. Two shorthands that are not part of the GL String grammar are accepted for typing reports: `A*03:01+`, with the second copy empty, is one allele typed and is read as homozygous; and a locus written without `+` (`DRB4*01:01:01G`, or `DRB3*01:01|DRB3*01:02`) is one allele typed, as laboratories often report DRB3/4/5. Phased haplotypes (`~`) and possible gene locations (`?`) are not read.

### GLS: pull and GLID files (`.pull`, `.glid`)

Used by registries and laboratories that report through GL identifiers. The GLID file numbers each locus GL string; the pull file lists, per person, one GL ID per locus:

```text
populationData.glid          populationData.pull
1;A*01:01+A*02:01            s1;1:2
2;B*08:01+B*07:02/B*07:03    s2;3:0
3;A*03:01+A*11:01
```

GL ID `0` means the locus was not typed; the person is kept with that locus untyped. The GLID file is found beside the pull file (same name, `.glid`), or given with `--glid`.

### HML (`.hml`)

NMDP Histoimmunogenetics Markup Language, versions 1.0 and 1.0.1. Each `<sample id=...>` gives one person; each `<allele-assignment>` gives one locus, from its `<glstring>` or, when there is none, from its `<haploid>` alleles. A folder of HML files (often one per person) can be given as the input. Genotype lists (`|`) and allele lists (`/`) are kept in full and weighed by the EM, so no ambiguity is resolved by hand. A GL string given only by its URI (`<glstring uri="..."/>`) cannot be read offline; the locus is left untyped and listed, so resolve such URIs with the GL service and give the GL string as text.

In any GL string input (GL string tables, HML, NGSengine), a phased GL string, with `~` joining the loci of one haplotype (`A*01:01~B*08:01+A*02:01~B*07:02`), is read without its phase: each locus becomes an unphased genotype, and the person is listed with a note. The estimate then infers the phase like any other.

### NGSengine XML export (`.xml`)

The project export of GenDx NGSengine (`<ProjectXml>`, the TARR layout). For each `<Locus>`, the reviewed result in `<ConclusiveTypingResult>` is read for its `<SampleName>`. Locus results whose review status is not `Approved` (for example `Rejected`) and loci with no conclusive result are left untyped, and each is listed on screen and in `run.json`; `--include-unapproved` reads them anyway. Several exports, or a folder of them, can be read together. The file is read as a stream, so a 150 MB export with sequences takes a few seconds and little memory.

The layout of `.xml` and `.hml` files is recognised from their content; `--format hml` or `--format ngsengine` overrides it. If the same person and locus appear twice with different results, the locus is left untyped and listed.

### Genotype list (`--format genotypes`)

The resolved genotypes written by `estimate --write-genotypes` (Hapl-o-Mat's `WRITE_GENOTYPES` output and `READ` input): one line per possible multilocus genotype of a person, `id<TAB>type<TAB>weight<TAB>A*01:01g+A*02:01g^B*08:01g+B*07:02g`, a person's weights summing to 1 and an unknown copy written `A*?`. The codes are already resolved, so estimate it with `asItIs` at every locus; the result is the same as from the original input. A person's lines must be all combinations of per-locus genotypes with weights that multiply (as written by either tool); others are left out with a warning.

### Serological antigens and pre-2010 names

A typing that is a serological antigen entered as an allele (`B*62`, `DRB1*15`, or `62` in a B column) or a pre-2010 name (`A*0201`, `A*010101`) is not an IMGT/HLA allele name. By default such a person is left out, with a message saying what the typing looks like; a one-field typing that is neither (`DRB1*92`) is reported as such too. `--serology` reads antigens as the alleles that carry them, from the WMDA tables (`rel_dna_ser.txt`, with broad antigens covering their splits and associated antigens from `rel_ser_ser.txt`), each allele with an equal share; `--v2-names` reads pre-2010 names as their current names (`Nomenclature_2009.txt`). An allele group that exists (`A*02`) is always read as alleles, never as an antigen. Antigen typings carry little information at the allele level, so expect broad frequencies for the people typed this way.

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
- A `Sample_ID` column, if present, names people; otherwise `family/member` is used (e.g. `fam01/C1`).

## Resolutions

Each locus is resolved to a target resolution, given as `LOCUS:RESOLUTION` (e.g. `--loci A:g,B:g,DRB1:2f`):

| Code | Target | Example for `A*02:01:01:01` |
| ---------- | --------------------------------------------------------------- | ----------------------------- |
| `1f` | First field | `A*02` |
| `2f` | Two fields | `A*02:01` |
| `3f` | Three fields | `A*02:01:01` |
| `4f` | Four fields | `A*02:01:01:01` |
| `g` | g group (identical antigen-recognition domains, nulls included) | `A*02:01g` |
| `G` | G group (identical exons 2 and 3, or 2 for class II) | `A*02:01:01G` |
| `P` | P group (identical protein in the antigen-recognition domains) | `A*02:01P` |
| `lgx` | g group named by its two fields, without the `g` (py-ard's `lgx`; `lg` is the same as `g`) | `A*02:01` |
| `asItIs` | Typings used exactly as written, without lookup | `A*02:01:01:01` |

A typing that covers several alleles at the target resolution (e.g. `A*02` at `2f`, or a MAC code) contributes to each with equal weight. Expression suffixes (N, L, S, Q) are kept. Every locus in the IPD-IMGT/HLA allele list can be resolved, including MICA, MICB and the non-classical and pseudogene loci (E, F, G, H, J, K, L, DQA2, TAP1, ...), so haplotypes can extend beyond the classical loci. IMGT/HLA defines G and P groups only for the classical and E/F/G-type loci; for the others, use `1f` to `4f` or `g` (the 2-field name). Loci not in the list at all can still be used with `asItIs`.

## Command line

```bash
uv run hlahaplotools estimate  INPUT --loci A:g,B:g,... [-o results]   # population data
uv run hlahaplotools families  INPUT --loci A:g,B:g,... [-o results]   # family data
uv run hlahaplotools haplomat  MAC|GLSC|GLS|READ PARAMETER_FILE        # run from a parameter file
uv run hlahaplotools stats     INPUT --loci A:2f,B:2f,... [-o results/stats]  # Hardy–Weinberg and LD
uv run hlahaplotools impute    INPUT --frequencies hfs.dat [-o results/impute]  # each person's haplotype pairs
uv run hlahaplotools summary   INPUT                                   # describe an input file
uv run hlahaplotools rename    INPUT OUTPUT                            # update retired allele names
```

Options for `estimate` and `families`:

| Option                                              | Default                     | Meaning                                                                                                                               |
| --------------------------------------------------- | --------------------------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| `--loci`                                          | required                    | Loci and resolutions, e.g. `A:g,C:g,B:g,DRB1:g`                                                                                     |
| `-o`, `--output`                                | `results`                 | Output folder                                                                                                                         |
| `--format`                                        | from extension or content   | `mac`, `glsc`, `gls`, `hml`, `ngsengine` or `genotypes` (`estimate` and `summary`)                                      |
| `--glid`                                          | pull file with `.glid`    | GLS input: the GLID file                                                                                                              |
| `--include-unapproved`                            | off                         | NGSengine input: also read locus results not approved in review                                                                       |
| `--epsilon`                                       | `1e-6`                    | Stop when no haplotype frequency changes by more than this                                                                            |
| `--init`                                          | `numberOccurrence`        | Start values: `numberOccurrence`, `perturbation`, `equal`, `random`                                                           |
| `--seed`                                          | `0`                       | Seed for `perturbation` and `random`                                                                                              |
| `--starts`                                        | `10`                      | EM start values: `--init` first, then perturbed starts; the best fit is kept (`estimate` only)                                    |
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

`haplomat` runs from an existing Hapl-o-Mat parameter file, with a single EM start, so earlier analyses can be repeated unchanged. `WRITE_GENOTYPES=true` writes the genotype list to `FILENAME_GENOTYPES`, and `READ` estimates from such a list. The ambiguity filter (`DO_AMBIGUITYFILTER`) and the haplotype analytics (`DO_ANALYTICS`) are not implemented, and are ignored with a warning: ambiguity is kept in full and weighed by the EM instead.

`impute` gives each person's most likely haplotype pairs and multilocus genotypes from a frequency table (`hfs.dat` from `estimate`, or any `haplotype<TAB>frequency` file), as GRIMM and py-graph-imputation do. Loci and resolutions are taken from the `run.json` beside the frequency file unless `--loci` is given, so the typing is resolved as it was for the estimate. Repeat `--frequencies NAME=FILE` for several populations (with `--prior NAME=WEIGHT`, default equal) to also get each person's population posteriors. `--top` (default 10) sets how many pairs and genotypes are written per person; `--max-candidates` (default 5,000) caps the haplotypes considered per person, keeping the most frequent.

`rename` rewrites a typing file (MAC, GLSC or family table) with current IMGT/HLA names. Alleles deleted or renamed since the typing was reported (e.g. `C*04:09N`, now `C*04:09L`, or `B*47:01:01:01`, now `B*47:01:01:03`) are replaced with the name `Deleted_alleles.txt` gives for them, following chains of renames, and everything else is written unchanged. The changes are listed in `OUTPUT.renames.tsv`. Without this step, a person with a retired name is left out of the estimate and reported in `run.json`.

`stats` takes the same input options as `estimate`, plus `--permutations` (random tables per Hardy–Weinberg test, default 10,000), `--pairs all|adjacent` (every pair of loci, or only neighbours in `--loci` order) and `--ld-permutations` (shuffled datasets per linkage-disequilibrium test; default 0, no P-value).

The exit status is 0 on success and 2 if the EM did not converge.

## Output files

### `estimate`

| File            | Contents                                                                                                                                                                                                                                                                                                                                                                   |
| --------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `hfs.dat`     | `haplotype<TAB>frequency`, most frequent first. Alleles are joined with `~` in the order given by `--loci`.                                                                                                                                                                                                                                                          |
| `epsilon.dat` | One line per iteration: largest frequency change, log-likelihood                                                                                                                                                                                                                                                                                                           |
| `run.json`    | Tool version, input, loci and resolutions, IMGT/HLA release, MAC table date, all parameters, how blank alleles were filled (rules and counts per locus), numbers of samples read, used and incomplete, iterations, log-likelihood, timings, and every unreadable line, locus result left untyped while reading (`loci_left_untyped`) and left-out sample with the reason |

### `families`

| File                         | Contents                                                                                                                                                                                 |
| ---------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `flags.tsv`                | Pedigree problems: family, member, relative, kind, loci, explanation (see [Pedigree checks](#pedigree-checks))                                                                            |
| `mendel.tsv`               | Per child and locus: typing that cannot come from the parents; families that cannot be phased                                                                                            |
| `phases.tsv`               | Per family: phasing status, the father's haplotypes `hap_A` and `hap_B` and the mother's `hap_C` and `hap_D` (`LOCUS*?` where unknown), inconsistent and ambiguous loci, notes |
| `children.tsv`             | Per child, in child-number order: `paternal` (A or B), `maternal` (C or D) and both together (e.g. `BC`); `?` where the data do not decide                                       |
| `hfs.dat`, `epsilon.dat` | Haplotype frequencies estimated from founders only                                                                                                                                       |
| `diplotypes.tsv`           | Most likely haplotype pair per founder, with its posterior probability                                                                                                                   |
| `comparison.json`          | EM on founders as unrelated people compared with the family phasing (meaningful for larger cohorts)                                                                                      |
| `run.json`                 | As for `estimate`, plus family counts, phasing status counts and flag counts                                                                                                           |

### `impute`

| File               | Contents |
| ------------------ | -------- |
| `pairs.tsv`        | Per person: rank, population, the two haplotypes and the probability of the pair, most likely first |
| `genotypes.tsv`    | Per person: rank, multilocus unphased genotype (`A*01:01+A*02:01^B*...`) and its probability, summed over the pairs that give it |
| `populations.tsv`  | With several populations: each person's posterior for each |
| `run.json`         | Inputs, frequency files and priors, loci, parameters, people imputed, and the people left out with the reason (unreadable typing, or no pair in the table that fits it) |

### `stats`

| File              | Contents |
| ----------------- | -------- |
| `hwe.tsv`         | Per locus: people tested, alleles, observed and expected heterozygosity, F, exact-test P, its false-discovery-rate q-value, and people not tested (untyped, half-typed, ambiguous) or read as homozygous from a blank |
| `ld.tsv`          | Per pair of loci: people used, alleles at each locus, D′, Wn, likelihood-ratio statistic, nominal degrees of freedom, and the permutation P-value and its q-value (with `--ld-permutations`) |
| `ld_alleles.tsv`  | Per pair of loci and allele pair: haplotype frequency, allele frequencies, D and D′, most frequent first |
| `run.json`        | Input, loci and resolutions, reference data versions and parameters |

## Python API

```python
import hlahaplotools as h

nom = h.Nomenclature("data/derived/3.66.0")

# Population data
ds = h.read_mac("population.dat", mac_codes="data/mac/alpha.v3.txt")
ds.sample(0)  # one person: {locus: GL string, or None if untyped}
ds.fill_counts()  # what was filled in at each locus: {locus: {"homozygous": n, "absent_gene": n, ...}}
ds.missing_counts()  # untyped people per locus
h.read_mac("population.dat", half_typed="unknown")  # keep a blank second allele unknown
est = h.estimate(ds, nom, {"A": "g", "C": "g", "B": "g", "DRB1": "g", "DQB1": "g"}, epsilon=1e-6)
est, starts = h.estimate_multistart(ds, nom, "g", starts=10)  # best of 10 starts, as the command line does
est.frequencies(cut=1e-6)[:3]  # [(haplotype, frequency), ...]
est.best_diplotypes()[:3]  # [(sample, haplotype 1, haplotype 2, probability), ...]
est.issues  # samples left out, with reasons
est.write_frequencies("hfs.dat")

# Family data (loci covered by the IMGT/HLA tables)
loci = ["A", "C", "B", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
fam = h.read_family_table("families.csv", loci=loci)
h.audit_families(fam, nom, "2f")  # pedigree problems
h.check_mendelian(fam, nom, "2f")  # per-locus Mendel failures
phases = h.phase_families(fam, nom, "2f")  # one FamilyPhase per family
phases[0].labelled_haplotypes  # the four parental haplotypes: {"A": ..., "B": ..., "C": ..., "D": ...}
phases[0].children  # [(child, paternal haplotype "A"/"B", maternal haplotype "C"/"D"), ...]
est, phases = h.estimate_families(fam, nom, "2f")
```

`resolutions` may be one code for all loci (`"g"`), a list in locus order, or a dict. `h.hardy_weinberg(ds, nom, "2f")` and `h.linkage_disequilibrium(ds, nom, "2f", pairs=[("A", "B")])` return the statistics of `stats` as lists of dicts. `read_glsc` reads GL string tables, `read_gls(pull, glid)` GLS input, `read_hml(paths)` HML files and `read_ngsengine(paths, include_unapproved=False)` NGSengine exports (one path or a list); `ds.notes` lists locus results left untyped while reading, and `h.sniff_xml(path)` tells HML from NGSengine. All readers take `half_typed="homozygous"|"unknown"` and `drb345=True|False`. Passing `mac_codes` as a path loads only the codes the input uses; `h.MacCodes(path)` loads the whole table for reuse. Type hints are in [`_core.pyi`](src/hlahaplotools/_core.pyi).

## Methods

### Typing resolution

Each typed allele, group or MAC code is expanded to full IMGT/HLA allele names and mapped to the target resolution. When a typing covers several target codes, each receives an equal share of its weight; when a side of a genotype lists several candidates (`02:01/02:02`), each candidate receives an equal share; each alternative of a GL-string genotype list (`|`) receives weight 1/n. A sample is expanded into all combinations of its per-locus genotypes, each weighted by the product of these shares. A sample with more than `1/min_genotype_frequency` combinations is left out and reported.

### EM algorithm

Haplotype frequencies are estimated by maximum likelihood with the EM algorithm. The likelihood of a multi-locus genotype is the sum, over ordered haplotype pairs (h, h′) consistent with it, of f(h)·f(h′). The E step computes the expected number of copies of each haplotype; the M step divides by twice the number of samples. Iteration stops when no frequency changes by more than `epsilon`; frequencies below 1e-14 are set to zero. Start values are either haplotype occurrence counts over all phase combinations (`numberOccurrence`, the default), those counts with up to 10% seeded noise (`perturbation`), equal frequencies, or seeded random frequencies.

The EM converges to a stationary point near its start, which need not be the maximum. On data with several loci or many rare haplotypes, the occurrence-count start often stops on a symmetric point far below it: when the possible haplotypes of a person are carried by nobody else, an even start stays even. `estimate` therefore runs 10 starts by default, all on one model built once: the `--init` start, then 9 perturbed starts (occurrence counts with seeded noise, seeds `seed + 1` to `seed + 9`), and keeps the highest log-likelihood. Every start's log-likelihood is recorded in `run.json`, and the screen output says whether the starts disagreed. `--starts 1` reproduces a single run.

### Blank alleles

Blanks are filled when the file is read, before resolution, with the rules in [Input formats](#input-formats); `DRBx*NNNN` is an allele like any other and stays `DRBx*NNNN` at every resolution. In unrelated people the homozygous reading cannot be checked. In families it can: after resolution, each family and locus where no inheritance pattern fits is checked again with copied alleles read as unknown. The copies of the fewest members that make the locus fit are undone; when several choices fit equally, all their members are undone. A note names the locus and the members (in `phases.tsv`). When undoing copies does not help, they are kept and the locus is reported as inconsistent as usual. The rules are in [fill.rs](crates/hlahaplotools-core/src/fill.rs).

### Untyped and half-typed loci

A sample untyped at a locus is consistent with any allele there; a sample typed `a+?` requires `a` on at least one of its two haplotypes. Writing "a on one of the two" as [h = a] + [h′ = a] − [h = a and h′ = a] and multiplying over loci expresses the likelihood as a signed sum of products F(q)·F(q′), where q and q′ are haplotype *patterns* with open loci and F(q) is the total frequency of the haplotypes matching q. This gives exact E and M steps, which reduce to the standard ones for fully typed samples. Candidate haplotypes are those fully specified by the data. A pattern that matches none of them is filled in at its open loci. By default (`--completion linked`) an open locus takes only the alleles that go with the pattern's allele at its most strongly linked typed locus: the known linked pairs first (DQA1 with DQB1, DPA1 with DPB1, DRB3/4/5 with DRB1, C with B, MICA and MICB with B), then the nearest typed locus along the MHC. Which alleles go together is estimated by a two-locus EM of the pair over everyone typed at both, keeping alleles with a conditional probability of at least 0.1%; when nothing links, every allele seen at the locus is used. `--completion all` fills in every combination of the alleles seen (at most `max_completions` per pattern), which can create millions of haplotypes when several loci are untyped together. Completions that only one person's data could support cannot be told apart by the EM, so the exhaustive set adds candidates without information: in simulations with DQA1 and DPA1 untyped in 30% of people, both settings gave the same accuracy, with linked completion using a third of the candidates. The derivation is documented in [estimate.rs](crates/hlahaplotools-core/src/estimate.rs).

### Family phasing

The father's two haplotypes are labelled A and B, the mother's C and D. Children are taken in order of their member number (`C`, `C1`, `C2`, …, `C10`), whatever the row order in the file, and the lowest-numbered typed child defines the labels: it inherits A from the father and C from the mother. Every other child is then described as A or B plus C or D (e.g. `BC`). Inheritance is decided from all loci together, so a parent homozygous at some loci is resolved by the loci where its haplotypes differ. A parent whose two haplotypes are identical at every typed locus passes the same haplotype to every child; its children are then all reported as receiving A (father) or C (mother), with a note. A side stays `?` only when the data cannot decide it, for example when both parents carry the same two haplotypes and a child has one of each.

For each family, every way the children could have inherited one haplotype from each parent is tried (up to the 8 lowest-numbered children), and the inheritance patterns that explain the most loci are kept. At each locus the parental alleles consistent with every kept pattern are collected: a parent's alleles are assigned to haplotypes where all patterns agree, and left unknown (`LOCUS*?`) otherwise. An untyped parent is rebuilt from the children as far as they show it; when neither parent is typed, which rebuilt parent is the father is arbitrary. Loci no pattern explains are reported as inconsistent (typing error or recombination). All comparisons are made at the target resolution.

### Estimation from families

Relatives share haplotypes, so only founders (parents, typed or rebuilt) enter the estimate. Loci the family phased enter with known phase; other loci use the founder's own typing; unknown alleles of a rebuilt parent are treated as untyped. A family with no usable founders contributes its first child.

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

Each child is compared with haplotypes phased from the rest of its family, leaving out members already shown not to belong. Parent–child genotype matches also happen naturally (when the other parent carries the same haplotype), so a match is only reported when the child also fails to fit the other parent. Identical genotypes between unrelated people with common haplotypes are not rare, so cross-family duplicates are only reported on request.

### Hardy–Weinberg and linkage disequilibrium

**Hardy–Weinberg.** Only people whose genotype at the locus is unambiguous at the target resolution are tested; untyped, half-typed and ambiguous people are counted apart. Expected heterozygosity is Nei's unbiased estimate, 2n/(2n−1)(1 − Σp²), and F = 1 − Ho/He (positive when heterozygotes are too few). The exact test compares the probability of the observed genotype table, given the allele counts, with that of tables formed by pairing the same alleles at random; P is the share of random tables at most as probable as the observed one. With large samples even small departures are significant, so read F alongside P. A blank second allele read as homozygous adds homozygotes; their number is reported, and `--half-typed unknown` leaves those people out.

**Linkage disequilibrium.** For each pair of loci, haplotype frequencies are estimated by EM from everyone typed at either locus. With allele frequencies p and q and D = p_ab − p_a q_b: D′ for an allele pair divides D by its largest possible value with the same sign; the overall D′ (Hedrick) is Σ p_a q_b |D′_ab|; and Wn (Cohen) is √(Σ D²/(p_a q_b) / (min(k, l) − 1)), Cramér's V of the haplotype table, 0 under equilibrium and 1 when each allele at one locus goes with one at the other. With `--ld-permutations N`, linkage equilibrium is tested by permutation: one locus is shuffled across people, which keeps each locus's own data but breaks any association, the two-locus likelihood is re-estimated, and P is the share of the N shuffled datasets that fit at least as well as the real one (each needs an EM, so this takes time). The likelihood-ratio statistic 2(logL₂ − logL_A − logL_B) is also written, with its nominal (k − 1)(l − 1) degrees of freedom, but no χ² P-value: HLA loci at two fields have hundreds of mostly rare alleles, and the χ² approximation then fails.

**Multiple testing.** A run tests many loci and pairs at once, so `hwe.tsv` and `ld.tsv` also give `q_value`, the Benjamini–Hochberg false-discovery-rate adjusted P across the loci or pairs of that run; report it rather than the raw P. A permutation P cannot be smaller than 1/(N + 1), so choose N large enough for the corrected threshold: with 55 pairs, 999 permutations give P ≥ 0.001, above a Bonferroni threshold of 0.05/55 ≈ 0.0009, while 1,999 permutations give P ≥ 0.0005.

### Imputation

For a person with possible genotypes g (with weights w, from their typing), every pair of haplotypes h₁, h₂ in the frequency table gets P(h₁, h₂) ∝ f(h₁) f(h₂) Σ_orders Π_loci Σ_g w · [h₁, h₂ fit g], over both orders when h₁ ≠ h₂, normalised over the person's pairs. Untyped and half-typed loci constrain only the typed allele, so the pairs impute the rest. With several populations, both haplotypes of a pair come from one population, weighted by its prior, and the population posteriors are the shares of the pairs from each. Only haplotypes in the table can be proposed, so a person whose typing no pair fits gets no result. Candidates are found from the most selective typed locus and checked at every locus, so about 15,000 people take a few seconds with four loci and under half a minute with eleven. A frequency table cut at a small frequency (`--cut`) lacks the rarest haplotypes, so people who carry one are left out; impute from an uncut table (`--cut 0`) to keep them.

### Determinism

Work is split into fixed-size chunks that do not depend on the thread count, and partial sums are combined in a fixed order, so results are bit-identical on any number of threads.

## Validation

All results below are reproducible with the scripts in [scripts/](scripts). Simulated populations draw real allele names from IMGT/HLA 3.66.0 and give haplotypes skewed (Zipf-like) frequencies. Error is the total variation distance between estimated and true haplotype frequencies (0 = identical).

### Agreement with an established implementation

From the same start, frequencies agree with Hapl-o-Mat's to within 10⁻¹⁰ on its example populations in every input format it reads (`scripts/compare_haplomat.py`, run by the test suite).

### Incomplete typing

`uv run python scripts/validate_missing.py` — 5,000 people, 5 loci, 300 true haplotypes, 3 replicates (mean ± standard deviation). C and DQB1 are left untyped at the given rate; "drop incomplete" estimates from fully typed people only.

| C and DQB1 untyped | People left if dropped | Error, all people | Error, drop incomplete | Error, before blanking |
| ------------------ | ---------------------- | ----------------- | ---------------------- | ---------------------- |
| 5%                 | 4,502                  | 0.049 ± 0.002    | 0.052 ± 0.003         | 0.049 ± 0.002         |
| 20%                | 3,188                  | 0.049 ± 0.001    | 0.062 ± 0.003         | 0.049 ± 0.001         |
| 40%                | 1,820                  | 0.057 ± 0.005    | 0.086 ± 0.008         | 0.053 ± 0.005         |

With 15% of the remaining C and DQB1 typings also half-typed, errors were 0.048 against 0.062 when dropping (10% untyped) and 0.054 against 0.082 (30% untyped). With B, C and DQB1 each 30% untyped: 0.054 against 0.087.

### Accuracy at the scale of a real cohort

`uv run python scripts/plasmode.py TRUTH_HFS COHORT -o OUT --reps 3` (plasmode simulation). The frequencies estimated from a real cohort of 15,633 people are taken as the truth; new cohorts of the same size are drawn from them, with each simulated person untyped at the loci one real person was untyped at (about 30% at DQA1 and DPA1), and re-estimated from everyone and from the fully typed people only. Three replicates, mean (range); "common" haplotypes are those with a true frequency of at least 0.1%.

| Loci | True haplotypes | Error, everyone (15,633) | Error, fully typed only (10,873) | Common haplotypes: mean relative error, everyone / fully typed |
| ---- | --------------- | ------------------------ | -------------------------------- | -------------------------------------------------------------- |
| 8    | 38,645          | 0.320 (0.317–0.322)     | 0.356 (0.355–0.357)             | 16% / 18%                                                      |
| 11   | 31,980          | 0.307 (0.306–0.308)     | 0.344 (0.339–0.348)             | 17% / 18%                                                      |

The total error is much larger than in the simulations above because real haplotype frequencies have a long tail: the 78 haplotypes at 0.1% or more carry only about a quarter of the frequency, and the rest is spread over tens of thousands of rare haplotypes that 31,000 chromosomes cannot resolve. Expect common haplotypes to be estimated to within about 15–20% of their frequency at this cohort size, and treat rare ones as uncertain. Using everyone lowered the error by 10–11%; in these simulations the untyped loci were assigned at random, so in real data, where typing is often not missing at random, the gain is larger.

### Hardy–Weinberg and linkage disequilibrium against PyPop

`hlahaplotools stats` and PyPop 1.4.2 (Guo–Thompson exact test; pairwise LD from emhaplofreq) on the same 5,000 people from a real cohort, typed at A, B, C, DRB1 and DQB1 at two-field resolution (86, 139, 67, 73 and 36 alleles):

| Measure | Agreement over the 10 locus pairs or 5 loci |
| ------- | ------------------------------------------- |
| D′ | within 0.0004 for every pair |
| Likelihood-ratio statistic 2(logL₂ − logL_A − logL_B) | within 2 units (of 3,300–37,000) for every pair |
| Wn | identical to three decimals for 8 pairs; A–B and B–DRB1 differ by 0.008, where PyPop's EM reached a likelihood about 0.5 units higher |
| Hardy–Weinberg | both reject equilibrium at every locus (PyPop P ≤ 9 × 10⁻⁵, HLAhaploTools P < 10⁻⁵) |

PyPop took 26 minutes, HLAhaploTools 56 seconds with 100,000 random tables per Hardy–Weinberg test.

### Estimation from families

`uv run python scripts/validate_family_estimates.py` : 1,000 families of 1–4 children per replicate, 10% of parents untyped, DQB1 20% untyped.

| Replicate | Founders, family phase | Founders, unphased | Everyone as unrelated |
| --------- | ---------------------- | ------------------ | --------------------- |
| 1         | 0.076                  | 0.080              | 0.084                 |
| 2         | 0.078                  | 0.083              | 0.085                 |
| 3         | 0.078                  | 0.082              | 0.082                 |

Treating every family member as unrelated also counts about 4,300 people for about 2,000 independent founders, which overstates the precision.

### Phasing real pedigrees

`uv run python scripts/validate_pedigrees.py examples/families/*.csv` : the 9 example families, 18 loci at 2-field resolution. Each family is phased with all members as the reference, then again with members left out; called parental alleles are compared with the reference. Blank alleles are filled as described in [Input formats](#input-formats), so DRB3/4/5 absent genes count as alleles. One family with a flagged child (a copied row) is excluded from these figures.

| Scenario                                         | Correct | Wrong | Alleles called |
| ------------------------------------------------ | ------- | ----- | -------------- |
| Both parents and one child (21 trios)            | 1,323   | 0     | 97.5%          |
| One parent left out, 2 children (rebuilt parent) | 910     | 0     | 76.1%          |
| One parent left out, 3 children (rebuilt parent) | 340     | 0     | 98.6%          |

Of the 1,273 alleles recovered for rebuilt parents, 1,266 are part of that parent's own typing, 4 are at loci where the parent is untyped, and 3 fill a blank second allele that the children show is not a homozygous copy. In the parentage test, all 171 pairings of a child with another family's parents fail the Mendel check (2 to 14 inconsistent loci), and 20 of 21 true children pass; the exception is the flagged child. With 2 children, the untransmitted haplotype of a missing parent is often never observed, which explains the lower recovery.

### Pedigree checks

`uv run python scripts/validate_audit.py --seed 1` (and `--seed 2`) : 2,000 families of 2–4 children, 6 loci, 15% of parents untyped, with problems planted in some children.

| Planted problem                       | Both parents typed | A parent untyped |
| ------------------------------------- | ------------------ | ---------------- |
| Child from another father             | 86–93% flagged    | 51–67%          |
| Child's row copied from the mother    | 84–89%            | 86–93%          |
| Recombinant child                     | 80–84%            | 41–49%          |
| One allele changed (typo)             | 97–99%            | 44–59%          |
| Clean children flagged (false alarms) | 1.2–1.4%          | 6.7–7.9%        |

## Limitations

- Phasing does not model recombination; a recombinant child shows up as inconsistent loci (or a `recombination` flag) and those loci stay unphased. Only the first 8 children of a family are used to phase it.
- A blank second allele is read as homozygous. If a laboratory leaves alleles blank because they were not typed, use `--half-typed unknown`. The DRB3/4/5 rules follow the usual DRB1 groups and do not cover rare haplotypes with other gene content.
- G and P resolutions are not available for MICA, MICB, H, J and the other loci without IMGT/HLA G and P groups.
- Each extra locus multiplies the number of phases per person, and loci left untyped in many people add many filled-in candidate haplotypes: 11 loci with incomplete typing need several GB of memory (see [Run time and memory](#recommended-workflow)).
- First-field typings (e.g. `A*11`) expand to hundreds of alleles at two fields. Estimation stays fast, but in a family a locus whose parent is untyped and whose child is typed this way can be left unphased (with a note), and the best haplotype pair of such a person is scored in one pass, so its probability is a lower bound.
- Excel files are not read; save them as CSV. Phased GL strings (`~`) and GL strings given only by URI in HML are not read.
- NGSengine exports are read from the `ConclusiveTypingResult` of each locus; other export layouts need converting to GLSC first.
- On highly ambiguous data, the EM likelihood can have several nearly equal maxima; results then depend slightly on the start values.
- With small cohorts, haplotype frequencies and the EM-versus-family comparison are not informative.

## Development

```bash
cargo test -p hlahaplotools-core       # Rust unit tests
uv run pytest                           # Python API, command line and data tests
uv run ruff check . && cargo clippy --all-targets
```

Tests that need the reference data are skipped until `scripts/build_nomenclature.py` has been run.

| Path                              | Contents                                                                                              |
| --------------------------------- | ----------------------------------------------------------------------------------------------------- |
| `crates/hlahaplotools-core/`    | Rust core: parsers, MAC codes, nomenclature tables, resolution, EM, families, pedigree checks, output |
| `src/lib.rs`                    | Python bindings (`hlahaplotools._core`)                                                             |
| `src/hlahaplotools/`            | Python package and command line                                                                       |
| `scripts/build_nomenclature.py` | Downloads the pinned reference data and builds the lookup tables                                      |
| `scripts/derive_tables.py`      | Builds the lookup tables from an IMGT/HLA release                                                     |
| `scripts/simulate.py`           | Simulates populations or families (`--format families`, optionally with planted problems)           |
| `scripts/validate_*.py`         | The validations above                                                                                 |
| `examples/`                     | Example population and family data                                                                    |
| `tests/`                        | Test suite                                                                                            |

## Citation

If you use HLAhaploTools, please cite this [repository](https://github.com/fmobegi/HLAhaploTools), and report the IPD-IMGT/HLA release and MAC table date recorded in `run.json`.
Please, also cite the original [Hapl-o-Mat](https://pmc.ncbi.nlm.nih.gov/articles/PMC5450239/) paper for the EM algorithm and other data preprocessing methods used in HLAhaploTools.

## License and acknowledgements

HLAhaploTools is released under the GNU General Public License, version 3 or later ([LICENSE](LICENSE)).

HLAhaploTools reimplements, in new Rust and Python code, some methods of [Hapl-o-Mat](https://github.com/DKMS/Hapl-o-Mat) (DKMS gGmbH, GPL-3.0): typing resolution and the lookup-table build, the EM algorithm, and its MAC, GL string, GLS and parameter-file input formats. Their rules follow Hapl-o-Mat's so that the logic and results agree.

Nomenclature data come from the [IPD-IMGT/HLA Database](https://www.ebi.ac.uk/ipd/imgt/hla/), and multiple-allele codes from the [NMDP](https://hml.nmdp.org/mac/).

This repository previously held an R package of the same name; its last version is tagged [`r-package-v0.99.0`](https://github.com/fmobegi/HLAhaploTools/tree/r-package-v0.99.0).

## Use of AI

> This software was developed with the assistance of Artificial Intelligence tools: GitHub Copilot, and Anthropic’s Claude Sonnet 5.5 and Claude Opus 5.5. They were used to port over existing code, generate code snippets, revise the authors' source code and documentation, and compare analysis test results between different tools, under the direction of the developer(s). All AI-generated material was extensively reviewed, tested and, where needed, corrected by the author to ensure accuracy, functionality, and compliance with project standards. Correctness was checked by an automated test suite, by agreement with Hapl-o-Mat on shared inputs, and by simulations with known true frequencies.
