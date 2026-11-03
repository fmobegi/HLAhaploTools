# Python API

Everything the `hlahaplotools` command does is available from Python: the command line is a thin layer over the same Rust core. Use the API when HLAhaploTools is one step of something bigger, such as a Snakemake or Nextflow pipeline, a Jupyter notebook, a loop over many locus sets or cohort subsets, or a new method built on its parts.

```python
import hlahaplotools as h
```

Type hints for every function and object are in [`_core.pyi`](../src/hlahaplotools/_core.pyi). The package is at version 0.1, so the API may still change before the first release.

> [!IMPORTANT]
> `h.estimate()` and `h.estimate_families()` run **one** EM start by default; the command line runs ten. Pass `starts=10` (or use `h.estimate_multistart()`) to get what the command line gives.

## A worked example

```python
import json
import hlahaplotools as h

release = json.load(open("data/manifest.json"))["imgthla"]["release"]
nom = h.Nomenclature(f"data/derived/{release}")
mac = h.MacCodes("data/mac/alpha.v3.txt")            # load once, reuse for every file

# Unrelated people
ds = h.read_mac("examples/population/population.dat", mac_codes=mac)
print(len(ds), ds.missing_counts())                  # people read, untyped per locus
print(ds.issues)                                     # lines that could not be read
res = {"A": "2f", "C": "2f", "B": "2f", "DRB1": "2f", "DQB1": "2f"}
est = h.estimate(ds, nom, res, starts=10)            # ten starts, as the command line does
print(est.converged, est.log_likelihood, est.n_samples)
for haplotype, frequency in est.frequencies(cut=1e-4)[:5]:
    print(f"{frequency:.4f}  {haplotype}")
pairs = est.best_diplotypes()                        # [(id, haplotype 1, haplotype 2, probability), ...]
est.write_frequencies("hfs.dat")

# Hardy–Weinberg and linkage disequilibrium
hwe = h.hardy_weinberg(ds, nom, res)                 # one dict per locus
ld = h.linkage_disequilibrium(ds, nom, res, pairs=[("A", "B"), ("DRB1", "DQB1")])
q = h.adjust_fdr([row["p_value"] for row in hwe])    # Benjamini–Hochberg q-values

# Families
classical = ["A", "C", "B", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
fam = h.read_family_table("examples/families/validated.csv", mac_codes=mac, loci=classical)
flags = h.audit_families(fam, nom, "2f")             # check these first
phases = h.phase_families(fam, nom, "2f")
for p in phases:
    print(p.family, p.status, p.children)            # children: (child, "A"/"B", "C"/"D")
fam_est, _ = h.estimate_families(fam, nom, "2f", starts=10)

# Impute haplotype pairs from known frequencies
rows, left_out = h.impute(ds, nom, [("cohort", est.frequencies(cut=0.0), 1.0)], res, top=3)
```

## Objects

| Object | What it holds |
| --- | --- |
| `Nomenclature(path)` | The IMGT/HLA lookup tables (`data/derived/<release>`). `translate(code, resolution)` resolves one typing; `serology=True` and `v2_names=True` read antigens and pre-2010 names. |
| `MacCodes(path)` | The NMDP multiple-allele-code table, loaded once. `expand(locus, code)` lists the alleles a code covers. |
| `Dataset` | Typings of unrelated people: `ids`, `loci`, `sample(i)` (`{locus: GL string or None}`), `missing_counts()`, `one_allele_counts()`, `fill_counts()`, and `issues` and `notes` (lines not read, loci left untyped). |
| `Estimate` | An EM result: `frequencies(cut)`, `best_diplotypes()`, `log_likelihood`, `converged`, `iterations`, `starts` (every start's fit), `n_samples`, `n_samples_incomplete`, `n_candidates`, `timings`, `issues` (people left out), `write_frequencies()`, `write_history()`. |
| `FamilyData` | A family table: `dataset` (the typings as a `Dataset`), `families()`, `unrelated`, `issues`. |
| `FamilyPhase` | One family's phasing: `haplotypes` (`[A, B, C, D]`), `labelled_haplotypes`, `children`, `status`, `inconsistent_loci`, `ambiguous_loci`, `notes`, `parent_id(0 or 1)`. |

## Functions

| Task | Functions |
| --- | --- |
| Read typings | `read_mac`, `read_glsc`, `read_gls(pull, glid)`, `read_hml(paths)`, `read_ngsengine(paths, include_unapproved=False)`, `read_genotypes`. All take `mac_codes=`, `loci=`, `half_typed="homozygous"\|"unknown"` and `drb345=True\|False`. `sniff_xml(path)` tells HML from NGSengine. |
| Estimate | `estimate(ds, nom, res, starts=, completion=, threads=, ...)`; `estimate_multistart(ds, nom, res, starts=10)` returns the best fit and a summary of every start. |
| Families | `read_family_table`, `check_mendelian`, `phase_families`, `audit_families` (the pedigree flags of `flags.tsv`), `estimate_families(use_family_phase=True)`, `compare_to_families(estimate, phases)`. |
| Statistics | `hardy_weinberg(ds, nom, res, permutations=10_000)`, `linkage_disequilibrium(ds, nom, res, pairs=None, permutations=0)`, `adjust_fdr(p_values)`. |
| Imputation | `impute(ds, nom, populations, res, top=10)`, with `populations` a list of `(name, [(haplotype, frequency), ...], prior)`. |
| Genotype lists | `resolved_genotypes`, `write_genotypes`, `read_genotypes` (Hapl-o-Mat's `WRITE_GENOTYPES` format). |
| One typing | `classify_typing("B*62")`: does it look like a serological antigen or a pre-2010 name? |

`res` (resolutions) is one code for every locus (`"2f"`), a list in locus order, or a dict (`{"A": "2f", "MICA": "asItIs"}`); the codes are those of the [README](getting-started.md#resolutions).

## Tips

- **Check `.issues`.** Readers and estimates do not raise on a bad line or an unusable person: they record it in `issues` on the `Dataset`, `FamilyData` or `Estimate`, as the command line records it in `run.json`.
- **Reuse the reference data.** Load `Nomenclature` and `MacCodes` once; a `mac_codes` given as a path is read again (only the codes the file needs) at every call.
- **Threads change speed, not results.** The work runs in Rust, in parallel and outside Python's global lock; results are bit-identical for any `threads=`.
- **Record what you ran.** The API writes no `run.json`; keep the IMGT/HLA release (`nom.source`), the MAC table date (`data/manifest.json`) and your parameters with the results.
- **More examples.** The scripts in [scripts/](../scripts), such as `validate_family_em.py` and `plasmode.py`, are working programs built on this API.
