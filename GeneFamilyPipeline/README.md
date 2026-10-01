# Asimina triloba gene family analysis
# Gene Family Evolution Pipeline
### Pawpaw-Focused Comparative Genomics

This repository contains scripts used to analyze **gene family evolution across Annonaceae and related magnoliid taxa**, with a primary focus on common pawpaw (*Asimina triloba*).

The pipeline integrates genome and proteome quality assessment, chromosome-anchored protein filtering, orthogroup inference, functional annotation, gene family evolution modeling, visualization of gene family expansions and contractions, and downstream functional investigation of rapidly evolving gene families.

The final stage focuses on identifying rapidly evolving *Asimina triloba* gene families associated with biological processes potentially relevant to **cold tolerance and adaptation to temperate environments**.

The workflow uses the following major tools:

- **BUSCO** — proteome completeness assessment
- **OrthoFinder** — orthogroup inference and species-tree estimation
- **InterProScan** — protein functional annotation and GO assignment
- **CAFE5** — gene family expansion and contraction modeling
- **R / R Markdown** — visualization and downstream functional analyses

---

# Dataset Scope

The final comparative analysis includes six Annonaceae and related magnoliid species:

| Species | Family | Role |
|---|---|---|
| *Annona cherimola* | Annonaceae | Annonaceae comparison |
| *Annona montana* | Annonaceae | Annonaceae comparison |
| ***Asimina triloba*** | **Annonaceae** | **Focal species** |
| *Lindera megaphylla* | Lauraceae | Magnoliid comparison |
| *Magnolia kwangsiensis* | Magnoliaceae | Magnoliid comparison |
| *Persea americana* | Lauraceae | Magnoliid comparison |

These taxa allow gene family evolution in *Asimina triloba* to be evaluated in the context of both closely related Annonaceae and more distantly related magnoliid lineages.

---

# Repository Structure

The analysis scripts included in this repository are:

```text
01.download_data.sh
02.genome_assembly_stats.sh
03.filter_anchored_proteins.sh

1a.busco.sh
1b.busco_summaries.sh

2a.orthofinder.sh
2b.prepare_interproscan_input.sh

3a.interproscan.sh
3b.filter_interpro_output.sh

4a.UpsetPrepCafe.R

5a.cafe.sh
5b.cafe_zeroRoot.sh
5c.countfamilychanges.sh
5d.visualizecafechanges.R

6a.topgoprep.sh

7a.coldtolerance.Rmd
```

Large genomic datasets and intermediate computational outputs are not stored directly in the repository.

---

# Pipeline Steps

## 01. Data Download and Preparation

```text
01.download_data.sh
```

Downloads and organizes the genome and protein sequence data required for the comparative analysis.

The resulting data are organized into directories used by subsequent stages of the pipeline.

Example directories include:

```text
downloads/
genome_sequences/
protein_sequences/
unfiltered_protein_sequences/
```

---

## 02. Genome Assembly Statistics

```text
02.genome_assembly_stats.sh
```

Calculates summary statistics for the genome assemblies used in the comparative analysis.

This step provides an initial assessment of assembly characteristics before downstream comparative genomic analyses.

---

## 03. Chromosome-Anchored Protein Filtering

```text
03.filter_anchored_proteins.sh
```

Filters protein datasets to retain proteins associated with chromosome-level or pseudomolecule-level genome sequences.

This step reduces potential artifacts caused by proteins associated with unplaced or unanchored genomic scaffolds and provides a more standardized set of proteins for downstream gene family analyses.

The filtered protein datasets are subsequently used for BUSCO and OrthoFinder analyses.

---

# 1. Proteome Completeness Assessment

## 1a. BUSCO Analysis

```text
1a.busco.sh
```

BUSCO evaluates the completeness of each protein dataset using conserved single-copy orthologs.

This analysis is used to assess proteome quality and evaluate the effect of chromosome-anchored protein filtering.

Output:

```text
busco_results/
```

---

## 1b. BUSCO Summaries

```text
1b.busco_summaries.sh
```

Collects BUSCO results across species into a summary table for easier comparison.

Output:

```text
busco_summary_table.tsv
```

---

# 2. Orthogroup Inference

## 2a. OrthoFinder

```text
2a.orthofinder.sh
```

OrthoFinder clusters proteins into **orthogroups**, representing sets of genes descended from a common ancestral gene.

The analysis also estimates phylogenetic relationships among the included species.

Major outputs are stored within:

```text
orthofinder/
```

Important OrthoFinder outputs used downstream include files such as:

```text
Orthogroups.GeneCount.tsv
SpeciesTree_rooted.txt
```

The orthogroup gene-count matrix forms the basis for subsequent gene family evolution analyses.

---

## 2b. Prepare InterProScan Input

```text
2b.prepare_interproscan_input.sh
```

Extracts and organizes proteins from inferred orthogroups for functional annotation with InterProScan.

This step connects OrthoFinder orthogroups with the individual protein sequences that will be functionally annotated.

Output:

```text
interproscan_input/
```

---

# 3. Functional Annotation

## 3a. InterProScan

```text
3a.interproscan.sh
```

Runs InterProScan on proteins associated with the inferred orthogroups.

InterProScan is used to identify protein domains and functional annotations, including:

- Pfam domains
- PANTHER classifications
- InterPro accessions
- Gene Ontology (GO) terms

Output:

```text
interproscan_output/
```

---

## 3b. Filter InterProScan Output

```text
3b.filter_interpro_output.sh
```

Processes and filters the InterProScan results for downstream analyses.

This stage identifies usable functional annotations and creates cleaned GO annotation datasets.

Outputs include:

```text
include_orthogroups.txt
all_go_unique.tsv
all_go_unique_clean.tsv
```

Transposable-element-associated orthogroups are identified separately so that they can be excluded from downstream biological interpretation where appropriate.

Associated outputs include:

```text
transposon_associated_orthogroups.tsv
excluded_orthogroups_with_reason.tsv
```

---

# 4. Orthogroup Overlap and CAFE Preparation

```text
4a.UpsetPrepCafe.R
```

This R script summarizes patterns of orthogroup sharing among the six species and prepares gene-family counts for CAFE5.

The script generates **UpSet plots** to visualize shared and lineage-specific orthogroups and prepares the gene-family count matrix required for gene family evolution modeling.

Outputs include:

```text
cafe_gene_families.tsv
cafe_gene_families_filtered.tsv

upset_plot.pdf
upset_plot_all_orthogroups.pdf
upset_plot_all_orthogroups.png
```

An ultrametric species tree is used for the subsequent CAFE analysis:

```text
SpeciesTree_ultrametric.txt
```

---

# 5. Gene Family Evolution Modeling

## 5a. CAFE5 Analysis

```text
5a.cafe.sh
```

CAFE5 models changes in gene family size across the species phylogeny using a stochastic birth-death model.

The analysis asks whether observed gene family expansions or contractions along particular lineages are greater than expected under the estimated background rate of gene gain and loss.

Results are stored in:

```text
cafe_default/
```

---

## 5b. CAFE5 Zero-Root Analysis

```text
5b.cafe_zeroRoot.sh
```

Runs an additional CAFE5 analysis using the zero-root configuration.

This analysis allows gene families with zero inferred copies at the root to be retained for downstream examination.

Results are stored in:

```text
cafe_zeroRoot/
```

---

## 5c. Count Gene Family Changes

```text
5c.countfamilychanges.sh
```

Processes the CAFE results and summarizes gene family changes across individual lineages.

Gene families are classified according to whether they are:

- expanded
- contracted
- unchanged

Branch-specific CAFE probabilities are also used to distinguish rapidly evolving gene families from background gene family changes.

The resulting summaries allow expansion and contraction patterns to be compared across all six species.

---

## 5d. Visualize Gene Family Changes

```text
5d.visualizecafechanges.R
```

Generates figures summarizing gene family expansions and contractions across the phylogeny.

Outputs include:

```text
gene_family_changes.pdf
gene_family_changes.png
```

These figures provide a lineage-level overview of gene family evolution and highlight differences in expansion and contraction patterns among species.

---

# 6. GO Annotation Preparation

```text
6a.topgoprep.sh
```

Prepares Gene Ontology annotations for downstream functional analysis of rapidly evolving gene families.

This step connects CAFE-identified orthogroups with GO terms derived from the InterProScan annotations.

The resulting tables provide the link between:

```text
CAFE orthogroups
        ↓
Asimina triloba genes
        ↓
InterProScan annotations
        ↓
GO terms
```

These annotations are then used in the final *Asimina triloba*-focused analysis.

---

# 7. Cold-Tolerance Candidate Analysis

```text
7a.coldtolerance.Rmd
```

The final stage of the workflow is implemented as an **R Markdown analysis**.

R Markdown is used so that the analysis code, filtering decisions, figures, results, and biological interpretation can be documented together in a reproducible format.

This analysis focuses specifically on **rapidly evolving gene families in *Asimina triloba*** and investigates whether their functional annotations are associated with biological processes potentially relevant to adaptation to temperate environments.

Candidate functional categories investigated include processes related to:

- cold acclimation
- carbohydrate metabolism
- lipid metabolism
- oxidative stress
- membrane-associated processes
- transmembrane transport
- receptor and signaling pathways
- stress and defense responses

Rapidly evolving orthogroups are linked back to their associated *Asimina triloba* genes and InterProScan annotations, including:

- Gene Ontology terms
- Pfam domains
- PANTHER classifications
- InterPro annotations

Both **rapidly expanded and rapidly contracted gene families** are considered in the final analysis.

The `.Rmd` file serves as the reproducible source document for this analysis and can be rendered to HTML for interactive viewing of the results.

---

# Major Pipeline Outputs

| Output | Description |
|---|---|
| `busco_summary_table.tsv` | BUSCO completeness summary across species |
| `Orthogroups.GeneCount.tsv` | Gene counts for OrthoFinder orthogroups |
| `include_orthogroups.txt` | Orthogroups retained after annotation filtering |
| `all_go_unique.tsv` | GO annotations recovered from InterProScan |
| `all_go_unique_clean.tsv` | Cleaned GO annotation dataset |
| `cafe_gene_families.tsv` | Gene-family count matrix prepared for CAFE |
| `cafe_gene_families_filtered.tsv` | Filtered CAFE input matrix |
| `SpeciesTree_ultrametric.txt` | Ultrametric species tree used by CAFE |
| `gene_family_changes.pdf` | Visualization of lineage-specific gene family changes |
| `gene_family_changes.png` | PNG version of gene family change visualization |
| `upset_plot.pdf` | Orthogroup overlap among species |
| `upset_plot_all_orthogroups.pdf` | Orthogroup overlap across the complete dataset |
| `7a.coldtolerance.Rmd` | Final pawpaw cold-tolerance candidate analysis |


# Important Notes

- The analysis focuses on gene family evolution in *Asimina triloba* within a comparative magnoliid framework.
- Both **gene family expansions and contractions** are examined.
- Branch-specific CAFE probabilities are used to identify rapidly evolving gene families.
- Chromosome-anchored protein filtering is performed before orthogroup inference to reduce potential artifacts from unplaced genomic scaffolds.
- InterProScan annotations provide GO, Pfam, PANTHER, and InterPro functional information for downstream interpretation.
- Transposable-element-associated orthogroups are identified and excluded from biological interpretation where appropriate.
- The final cold-tolerance analysis considers both rapidly expanded and rapidly contracted *Asimina triloba* gene families.
- The final analysis is retained as an `.Rmd` file because it combines reproducible R code, figures, filtering decisions, and biological interpretation in a single document.
- Large genome files, protein FASTA files, and computationally intensive intermediate outputs are not intended to be stored directly in the GitHub repository.

---

