# Reproducing the analyses

The analyses of the accompanying manuscript are scripts in [scripts/](../scripts); the results are reported in the manuscript, not here. Every benchmark runs from one command, into one folder, on a laptop or a cluster node. The benchmarks need Linux (on Windows, use WSL): they cap and measure memory with Linux tools, and Hapl-o-Mat is built with `make`.

```bash
uv run python scripts/benchmark.py all --out runs/run1 --real COHORT.dat --threads 1,8 --rscript Rscript
uv run python scripts/benchmark.py pack --out runs/run1 --real COHORT.dat   # runs/run1_metrics.tar.gz
```

`all` runs each step in turn (`run`, `run-extra`, `extended`, `multistart`, `bootstrap`, `full-cohort`, `stats`, `plasmode`, `completion`, `missing`, `haplostats`, `report`; each is also a subcommand) and records exit codes and times in `steps.tsv` and the machine in `environment.json`. Use `--steps` or `--skip` to run part of it, `--haplomat`, `--pypop` and `--rscript` to point at the other tools, and `--plasmode-reps` for more plasmode replicates. `pack` keeps the metrics (JSON, TSV, logs, estimated frequencies), leaves out every typing file, and replaces the cohort's sample IDs (from `--real`) wherever logs or `run.json` name an individual. The comparison tools are installed separately: Hapl-o-Mat (`make -C reference/Hapl-o-Mat`, or `--haplomat PATH`), PyPop (`--pypop PATH`) and the R package haplo.stats (`--rscript`).

| Script                           | What it checks                                                                        |
| -------------------------------- | ------------------------------------------------------------------------------------- |
| `compare_haplomat.py`          | Agreement with Hapl-o-Mat on its example populations (run by the test suite)          |
| `validate_missing.py`          | Estimation with untyped and half-typed loci, against simulated truth                  |
| `plasmode.py`                  | Accuracy at the size and missingness of a real cohort                                 |
| `completion_truth.py`          | Linked against exhaustive completion of untyped loci                                  |
| `validate_family_estimates.py` | Estimation from founders, with and without family phase                               |
| `validate_pedigrees.py`        | Phasing of real pedigrees with members left out                                       |
| `validate_family_em.py`        | The EM on family members treated as unrelated, scored against the pedigree haplotypes |
| `validate_audit.py`            | Pedigree checks on simulated families with planted problems                           |
| `simulate.py`                  | Simulated populations or families with known haplotype frequencies                    |
