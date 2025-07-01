# HLAhaploTools

<!-- badges: start -->
[![R-CMD-check](https://github.com/fmobegi/HLAhaploTools/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fmobegi/HLAhaploTools/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/fmobegi/HLAhaploTools/graph/badge.svg)](https://app.codecov.io/gh/fmobegi/HLAhaploTools)
<!-- badges: end -->

The goal of HLAhaploTools is to provide a comprehensive, user-friendly toolkit for:

- Standardizing and cleaning HLA genotype data from various sources
- Estimating HLA haplotype and allele frequencies using robust EM algorithms
- Visualizing HLA allele and haplotype distributions across populations
- Supporting both classical (A, B, C, DRB1, DQB1, DQA1, DPB1, DPA1, DRB3, DRB4, DRB5) and non-classical HLA genes (E, F, G, H, J, K, L, MICA, MICB, HFE, DMA, DMB, DOA, DOB)
- Facilitating downstream immunogenetics research and reporting

## Supported HLA Gene Names

HLAhaploTools recognizes the following gene names (case-insensitive):

- **Classical HLA genes:**  
  A, B, C, DRB1, DRB3, DRB4, DRB5, DQA1, DQB1, DPA1, DPB1
- **Non-classical HLA genes and related loci:**  
  E, F, G, H, J, K, L, MICA, MICB, HFE, DMA, DMB, DOA, DOB

Example of allowed non-classical gene names:  
E, F, G, H, J, K, L, MICA, MICB, HFE, DMA, DMB, DOA, DOB

## Installation

You can install **HLAhaploTools** as a stable release from Bioconductor with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("HLAhaploTools")
```

Development versions can be accessed directly from GitHub using one of the following methods:

```r

# Using pak

# install.packages("pak")

pak::pak("fmobegi/HLAhaploTools")

# Using devtools

# install.packages("devtools")

devtools::install_github("fmobegi/HLAhaploTools")

# Using remotes

# install.packages("remotes")

remotes::install_github("fmobegi/HLAhaploTools")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(HLAhaploTools)

# Standardize column names
cleaned <- standardize_colnames(raw_hla_data)

# Extract loci
loci <- extract_loci(cleaned)

# Collapse genotypes
geno <- collapse_genotypes(cleaned, loci)

# Estimate haplotype frequencies
freqs <- run_em_algorithm(all_diplotypes)

# Plot allele frequencies (supports classical and non-classical genes)
Plot_HLA_allele_frequency(freqs)

# Plot diversity for a non-classical gene
Plot_HLA_diversity(freqs, gene = "E")
```

## Samples Input Formatting

The input data should be a tabular file (Coma-separated file CSV, TAB-separates TSV/TXT, Excel worksheet XLSX/XLS, or similar) with the following structure and column names. Each row corresponds to a family member with their HLA typing data:

| FAMILY_ID | Family_Member | A_1   | A_2   | B_1   | B_2   | C_1   | C_2   | DRB1_1 | DRB1_2 | DRB3_1 | DRB3_2 | DRB4_1 | DRB4_2 | DRB5_1 | DRB5_2 | DQA1_1 | ... |
|-----------|---------------|-------|-------|-------|-------|-------|-------|--------|--------|--------|--------|--------|--------|--------|--------|--------|-----|
| family01  | F             | A*30:…| A*25:…| B*18:…| B*18:…| C*05:…| C*12:…| DRB1*… | DRB1*… | DRB3*… | DRB3*… | NA     | NA     | NA     | DRB5*… | DQA1*… | ... |
| family01  | M             | A*26:…| A*02:…| B*44:…| B*13:…| C*16:…| C*06:…| DRB1*… | DRB1*… | NA     | NA     | DRB4*… | DRB4*… | NA     | NA     | DQA1*… | ... |
| family01  | C1            | A*30:…| A*26:…| B*18:…| B*44:…| C*05:…| C*16:…| DRB1*… | DRB1*… | DRB3*… | DRB3*… | NA     | DRB4*… | NA     | NA     | DQA1*… | ... |
| ...       | ...           | ...   | ...   | ...   | ...   | ...   | ...   | ...    | ...    | ...    | ...    | ...    | ...    | ...    | ...    | ...    | ... |

*Note: Allele values can be presented in full gene format (e.g., `A*01:01`), or in allele-only format (e.g., `*01:01` or `01:01`).
Family members (`Family_Member` column) must be designated as one of `F` for father, `M` for mother, or `C1, C2, Cn`, ... for children.
If your data is not `family typing data`, you can exclude this column and provide a samples identifier column with unique values.*

### Accepted column names

HLAhaploTools accepts a flexible naming convention for input columns. The package will automatically recognize and standardize various column name formats:

- **Standard naming format**: `GENE_1` and `GENE_2` (e.g., `A_1`, `A_2`)
- **Alternative formats**:
  - `GENE.1` and `GENE.2` (e.g., `A.1`, `A.2`)
  - `GENE1` and `GENE2` (e.g., `A1`, `A2`)
  - `GENE-1` and `GENE-2` (e.g., `A-1`, `A-2`)
  - `GENE_ALLELE1` and `GENE_ALLELE2` (e.g., `A_ALLELE1`, `A_ALLELE2`)

The standardization process will convert all recognized column names to the preferred format (`GENE_1`, `GENE_2`). The package supports all classical and non-classical HLA genes as listed below:



<!-- | Description    | Column Name in file   |
|----------------|-----------------------|
| FAMILY_ID      | Family ID             |
| Family_Member  | Family_Member         |
| HLA*A          | A_1, A_2              |
| HLA*B          | B_1, B_2              |
| HLA*C          | C_1, C_2              |
| HLA*DRB1       | DRB1_1, DRB1_2        |
| HLA*DRB3       | DRB3_1, DRB3_2        |
| HLA*DRB4       | DRB4_1, DRB4_2        |
| HLA*DRB5       | DRB5_1, DRB5_2        |
| HLA*DQA1       | DQA1_1, DQA1_2        |
| HLA*DQB1       | DQB1_1, DQB1_2        |
| HLA*DPA1       | DPA1_1, DPA1_2        |
| HLA*DPB1       | DPB1_1, DPB1_2        |
| HLA*F          | F_1, F_2              |
| HLA*G          | G_1, G_2              |
| HLA*H          | H_1, H_2              |
| HLA*J          | J_1, J_2              |
| HLA*E          | E_1, E_2              |
| HLA*MICA       | MICA_1, MICA_2        |
| HLA*MICB       | MICB_1, MICB_2        | -->

<table style="border-collapse: collapse; line-height: 1.1; font-size: 12px;">
  <thead>
    <tr>
      <th style="border: 0.5px solid #ddd; padding: 4px 6px; text-align: left;">Description</th>
      <th style="border: 0.5px solid #ddd; padding: 4px 6px; text-align: left;">Column Name in file</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">FAMILY_ID</td>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">Family ID</td>
    </tr>
    <tr>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">Family_Member</td>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">Family_Member</td>
    </tr>
    <tr>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">HLA*A</td>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">A_1, A_2</td>
    </tr>
    <tr>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">HLA*B</td>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">B_1, B_2</td>
    </tr>
    <tr>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">HLA*C</td>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">C_1, C_2</td>
    </tr>
    <tr>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">HLA*DRB1</td>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">DRB1_1, DRB1_2</td>
    </tr>
    <tr>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">HLA*DRB3</td>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">DRB3_1, DRB3_2</td>
    </tr>
    <tr>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">HLA*DRB4</td>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">DRB4_1, DRB4_2</td>
    </tr>
    <tr>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">HLA*DRB5</td>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">DRB5_1, DRB5_2</td>
    </tr>
    <tr>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">HLA*DQA1</td>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">DQA1_1, DQA1_2</td>
    </tr>
    <tr>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">HLA*DQB1</td>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">DQB1_1, DQB1_2</td>
    </tr>
    <tr>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">HLA*DPA1</td>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">DPA1_1, DPA1_2</td>
    </tr>
    <tr>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">HLA*DPB1</td>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">DPB1_1, DPB1_2</td>
    </tr>
    <tr>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">HLA*F</td>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">F_1, F_2</td>
    </tr>
    <tr>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">HLA*G</td>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">G_1, G_2</td>
    </tr>
    <tr>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">HLA*H</td>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">H_1, H_2</td>
    </tr>
    <tr>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">HLA*J</td>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">J_1, J_2</td>
    </tr>
    <tr>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">HLA*E</td>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">E_1, E_2</td>
    </tr>
    <tr>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">HLA*MICA</td>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">MICA_1, MICA_2</td>
    </tr>
    <tr>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">HLA*MICB</td>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">MICB_1, MICB_2</td>
    </tr>
        <tr>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">HLA*HFE</td>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">HFE_1, HFE_2</td>
    </tr>
    <tr>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">HLA*DMA</td>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">DMA_1, DMA_2</td>
    </tr>
    <tr>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">HLA*DMB</td>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">DMB_1, DMB_2</td>
    </tr>
    <tr>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">HLA*DOA</td>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">DOA_1, DOA_2</td>
    </tr>
    <tr>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">HLA*DOB</td>
      <td style="border: 0.5px solid #ddd; padding: 4px 6px;">DOB_1, DOB_2</td>
    </tr>
  </tbody>
</table>

#
