# HLAhaploTools

HLAhaploTools assigns and estimates extended HLA and MHC haplotypes. In families it phases haplotypes from inheritance; in unrelated people it estimates haplotype frequencies with an expectation–maximisation (EM) algorithm, keeping people who are untyped or half-typed at some loci.

> [!CAUTION]
> **Research use only.** HLAhaploTools is research software. It is not a medical device and has not been validated for clinical use: do not use it, or its haplotype assignments, pedigree checks or frequencies, to make decisions about patient care, donor selection or transplantation. Results are statistical estimates and can be wrong. It is provided without warranty of any kind (see [LICENSE](LICENSE)).

## What it is and what it is not

HLAhaploTools **is**:

- a haplotype frequency estimator for unrelated people, reading laboratory reports (MAC and GL string tables, GL identifier files, NMDP HML, GenDx NGSengine XML) with all their ambiguity;
- usable at any IPD-IMGT/HLA locus, including MICA, MICB and the non-classical and pseudogene loci, at a resolution chosen per locus;
- a pedigree tool: Mendel checks, phasing of parents from their children (including untyped parents), estimation from founders only, and checks that flag non-paternity, sample swaps, copied rows, non-siblings and recombinants;
- a calculator for Hardy–Weinberg tests, linkage disequilibrium and imputation of haplotype pairs from known frequencies;
- compatible with [Hapl-o-Mat](https://github.com/DKMS/Hapl-o-Mat): it reads Hapl-o-Mat input formats and parameter files and, from the same start, gives the same frequencies.

HLAhaploTools **is not**:

- a clinical, donor-search or matching tool (see the caution above);
- an HLA typing or allele-calling program: it starts from typings that have already been assigned;
- a model of recombination: a recombinant child shows up as inconsistent loci (or a `recombination` flag), and those loci stay unphased;
- a way to make rare haplotypes certain: frequencies of rare haplotypes remain uncertain at the cohort sizes typical of HLA studies, whatever the tool.

## Contents

- [What it is and what it is not](#what-it-is-and-what-it-is-not)
- [Installation](#installation)
- [Usage](#usage)
- [Key strengths](#key-strengths)
- [Limitations](#limitations)
- [Contributing](#contributing)
- [Credits](#credits)
- [Citation](#citation)
- [License](#license)
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

`hlahaplo` is a short form of the command: `hlahaplo estimate ...` is the same as `hlahaplotools estimate ...`.

### Reference data

Typings are resolved against the IPD-IMGT/HLA nomenclature and the NMDP multiple-allele-code (MAC) table. Download them once:

```bash
uv run python scripts/build_nomenclature.py     # about 650 MB, mostly the MAC table
```

This writes `data/`:

| Path                   | Contents                                                                                                 |
| ---------------------- | -------------------------------------------------------------------------------------------------------- |
| `imgthla/<release>/` | `hla_nom_g.txt`, `hla_nom_p.txt`, `hla_nom.txt`, `Allelelist.txt`, `Deleted_alleles.txt`       |
| `derived/<release>/` | Lookup tables used for resolution: `AllAllelesExpanded.txt`, `LargeG.txt`, `P.txt`, `Smallg.txt` |
| `mac/alpha.v3.txt`   | NMDP multiple-allele codes                                                                               |
| `manifest.json`      | Release, source URLs, SHA-256 checksums, MAC table date                                                  |

The IMGT/HLA release is pinned (`PINNED_RELEASE` in [scripts/build_nomenclature.py](scripts/build_nomenclature.py), currently 3.66.0). To use another release, run with `--release 3.65.0`. The MAC table has no release number: it is downloaded once and kept, and its `LAST UPDATED` date is recorded; `--refresh-mac` downloads a new copy. Report both with any results.

Commands find `data/` in the current folder, in `$HLAHAPLOTOOLS_DATA`, or with `--data-dir`.

## Usage

```bash
uv run hlahaplotools estimate INPUT --loci A:2f,C:2f,B:2f,DRB1:2f,DQB1:2f -o results/population   # unrelated people
uv run hlahaplotools families FAMILY_TABLE --loci A:2f,C:2f,B:2f,DRB1:2f,DQB1:2f -o results/families  # families
```

Each command writes one folder:

```text
results/
├── population/          estimate
│   ├── hfs.dat          haplotype frequencies, most frequent first
│   ├── epsilon.dat      EM convergence
│   └── run.json         inputs, reference data versions, parameters, people used and left out
├── families/            families
│   ├── flags.tsv        pedigree problems: check these first
│   ├── mendel.tsv       typings a child cannot have inherited
│   ├── phases.tsv       each family's haplotypes: A, B (father), C, D (mother)
│   ├── children.tsv     which haplotypes each child inherited
│   ├── hfs.dat          frequencies from the founders
│   ├── epsilon.dat      EM convergence
│   ├── diplotypes.tsv   most likely haplotype pair per founder
│   ├── comparison.json  EM on the founders against the family phasing
│   └── run.json
├── stats/               stats: hwe.tsv, ld.tsv, ld_alleles.tsv, run.json
└── impute/              impute: pairs.tsv, genotypes.tsv, populations.tsv (several populations), run.json
```

`hlahaplotools --help` lists every command (`estimate`, `families`, `stats`, `impute`, `haplomat`, `summary`, `rename`). Example data are in [examples/](examples).

- **[Getting started](docs/getting-started.md)**: the user guide: quick start, input formats, resolutions, every command-line option, output files, good practice, and how it works.
- **[Python API](docs/python-api.md)**: if HLAhaploTools is one cog in a bigger machine, a Snakemake or Nextflow pipeline, a Jupyter notebook, or a loop over fifty locus sets at 2 a.m., call it from Python: copy the worked example and edit to taste.
- **[Reproducing the analyses](docs/reproducing-the-analyses.md)**: the benchmark and validation scripts of the manuscript.

## Key strengths

- **Laboratory reports in, ambiguity kept**: MAC, GL strings, GLS, HML and NGSengine, with MAC codes, allele lists and groups weighed by the EM.
- **Nobody left behind**: people untyped or half-typed at some loci stay in the estimate.
- **The whole MHC**: any IPD-IMGT/HLA locus, including MICA, MICB and the non-classical genes, at a resolution chosen per locus.
- **Families and unrelated people**: pedigree phasing and checks, founder-only estimates, and population estimates in one tool.
- **Reporting conventions built in**: blank second alleles and DRB3/4/5 filled as laboratories report them, and every fill counted.
- **Fast and reproducible**: a Rust core, ten EM starts by default, results identical on any number of threads, and every run recorded in `run.json`.
- **Hapl-o-Mat compatible**: its formats and parameter files are read, and from the same start the results agree.

## Limitations

- **Local maxima**: the EM can stop short of the best fit; ten starts reduce the risk but do not remove it, so check whether the starts agreed.
- **Rare haplotypes** remain uncertain at typical cohort sizes, and candidate haplotypes (and memory) grow fast with many loci and much untyped data.
- **Linked completion** is an approximation; `--completion all` is the exhaustive alternative.
- **Reporting conventions**: a blank second allele is read as homozygous unless `--half-typed unknown`, and the DRB3/4/5 rules do not cover rare DRB1 haplotypes with other gene content.
- **Families**: recombination is not modelled (recombinant loci stay unphased), and only the first 8 children of a family are used.
- **Imputation** proposes only haplotypes in the frequency table; it is not a graph search or donor matching.
- **Input**: no Excel files (save as CSV); phased GL strings are read without their phase; GL strings given only by URI are not read; NGSengine is read from the project export only; G and P groups are not available for MICA, MICB, H, J and other loci without them in IMGT/HLA.

## Contributing

Bug reports, questions and pull requests are welcome through [GitHub issues](https://github.com/fmobegi/HLAhaploTools/issues). For a bug, include the command, the screen output, `run.json`, and a small input that shows the problem. **Do not attach or commit typing data from your own cohorts**: reduce the problem to simulated data (`scripts/simulate.py`) or to the files in `examples/`.

Development setup and checks:

```bash
uv sync                                  # build the extension in development mode
cargo test -p hlahaplotools-core         # Rust unit tests
uv run pytest                            # Python API, command line and data tests
uv run ruff check . && uv run ruff format --check .
cargo clippy --all-targets && cargo fmt --check
```

Tests that need the reference data are skipped until `scripts/build_nomenclature.py` has been run. A change to estimation, resolution or file reading should come with a test, and should keep agreement with Hapl-o-Mat (`scripts/compare_haplomat.py`) unless the difference is intended and documented.

| Path                              | Contents                                                                                                 |
| --------------------------------- | -------------------------------------------------------------------------------------------------------- |
| `crates/hlahaplotools-core/`    | Rust core: parsers, MAC codes, nomenclature tables, resolution, EM, families, pedigree checks, output    |
| `src/lib.rs`                    | Python bindings (`hlahaplotools._core`)                                                                |
| `src/hlahaplotools/`            | Python package and command line                                                                          |
| `scripts/build_nomenclature.py` | Downloads the pinned reference data and builds the lookup tables                                         |
| `scripts/derive_tables.py`      | Builds the lookup tables from an IMGT/HLA release                                                        |
| `scripts/`                      | Benchmarks, validations and simulation (see [Reproducing the analyses](docs/reproducing-the-analyses.md)) |
| `examples/`                     | Example data: a simulated population, and family typings from the study described in the manuscript      |
| `tests/`                        | Test suite                                                                                               |

## Credits

HLAhaploTools was originally written by Fredrick M. Mobegi: [@fmobegi](https://github.com/fmobegi) at the Department of Clinical Immunology, [PathWest Laboratory Medicine WA](https://pathwest.health.wa.gov.au/), and [The University of Western Australia](https://www.uwa.edu.au/).

We thank the following people for their assistance in the development and testing of HLAhaploTools:

- Dianne De Santis
- Lloyd D'Orsogna
- [Hoiley Sham](https://github.com/h-sham)
- Samuel Bruce
- Sarah Nortmann
- Naser El-Lagta

HLAhaploTools reimplements, in new Rust and Python code, some methods of [Hapl-o-Mat](https://github.com/DKMS/Hapl-o-Mat) (DKMS gGmbH, GPL-3.0): typing resolution and the lookup-table build, the EM algorithm, and its MAC, GL string, GLS and parameter-file input formats. Their rules follow Hapl-o-Mat's so that the logic and results agree.

Nomenclature data come from the [IPD-IMGT/HLA Database](https://www.ebi.ac.uk/ipd/imgt/hla/), and multiple-allele codes from the [NMDP](https://hml.nmdp.org/mac/).

HLAhaploTools began in this repository as an R package for family segregation analysis; its last version is tagged [`r-package-v0.99.0`](https://github.com/fmobegi/HLAhaploTools/tree/r-package-v0.99.0). A copy of that prototype, with a related Shiny interface, is kept separately as a student's project at [h-sham/HLAhaploTools](https://github.com/h-sham/HLAhaploTools) and [h-sham/HLAhaploTools_Shiny](https://github.com/h-sham/HLAhaploTools_Shiny). These projects, while sharing the name, are not maintained with this repository, and do not share code with the current version.

Maintenance and future developments will be led by Fredrick M. Mobegi. Feel free to raise an [issue](https://github.com/fmobegi/HLAhaploTools/issues) or message the author if you need any support getting this tool running, or with suggestions for improvement.

<p float="center">
  <img src="docs/images/pathwest_logo.png" alt="PathWest Laboratory Medicine WA" width="400" height="90">
  <img src="docs/images/uwa_logo.png" alt="The University of Western Australia" width="400" height="90">
</p>

## Citation

A manuscript describing HLAhaploTools is in preparation. Until it is published, please cite this [repository](https://github.com/fmobegi/HLAhaploTools) with the version used, and report the IPD-IMGT/HLA release and MAC table date recorded in `run.json`. Please also cite the original [Hapl-o-Mat](https://pmc.ncbi.nlm.nih.gov/articles/PMC5450239/) paper for the EM algorithm and the data preprocessing methods that HLAhaploTools reimplements.

## License

HLAhaploTools is released under the GNU General Public License, version 3 or later ([LICENSE](LICENSE)).

## Use of AI

> This software was developed with the assistance of Artificial Intelligence tools: GitHub Copilot, and Anthropic’s Claude Sonnet 5.5 and Claude Opus 5.5. They were used to port over existing code, generate code snippets, revise the authors' source code and documentation, and compare analysis test results between different tools, under the direction of the developer(s). All AI-generated material was extensively reviewed, tested and, where needed, corrected by the author to ensure accuracy, functionality, and compliance with project standards. Correctness was checked by an automated test suite, by agreement with Hapl-o-Mat on shared inputs, and by simulations with known true frequencies.
