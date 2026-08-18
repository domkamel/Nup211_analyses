# Functional Analysis of the Fission Yeast Nuclear Basket Protein Nup211

## Rationale 

**Nup211** is a component of the nuclear basket in the fission yeast *Schizosaccharomyces pombe*.  
It belongs to a conserved family of nucleoporins that help form the "struts" of the nuclear basket structure. Other members of this family include the nucleoporin **TPR** in mammals, **Mlp1** and **Mlp2** in budding yeast, and **Alm1** and **Nup211** in *S. pombe*. While many of these proteins contribute to nucleocytoplasmic transport (e.g., mRNA quality control and export), numerous studies have also implicated them in regulating chromatin organization, transcription, DNA repair, and cell cycle progression. 

Our lab became particularly interested in fission yeast Nup211 because it is essential for viability, yet its molecular mechanisms and essential functions aren't well understood -- Nup211 regulates mRNA export, but how it does so isn't clear. Several studies have also shown that Nup211 associates with regions of chromatin, though whether it has roles outside of regulating nucleocytoplasmic transport remains an open question. 

We hypothesized that **Nup211 performs additional functions beyond mRNA export that are necessary for the viability of *S. pombe***.

---

## Our approach

Since Nup211 is essential, our lab generated a conditional depletion strain (*nup211-so*) in which *nup211* expression is controlled by the thiamine-repressible promoter P81nmt1. 

To investigate how Nup211 affects the transcriptome, We performed RNA-seq after depleting Nup211 protein for 16h and 48h. 

I also performed immunoprecipitation and proximity-labeling experiments to better understand its molecular mechanisms by identifying its interacting partners. 

---

## Contents

This repository contains the code used to analyze the RNA-seq, RT-qPCR, and mass spectrometry data collected during the course of my work. 

- [Repository structure](#repository-structure)
- [Data availability](#data-availability)
  - [Sequencing data](#sequencing-data-ncbi-geo)
  - [Processed RNA-seq files](#processed-rna-seq-files)
  - [Data found in this repository](#data-included-in-this-repository)
- [Thesis tables](#thesis-tables)
- [Requirements](#requirements)
- [Analyses](#analyses)
  - [RNA-seq - 16 h depletion (Salmon)](#rna-seq-16-h-depletion-salmon)
  - [RNA-seq - 48 h depletion (STAR)](#rna-seq-48-h-depletion-star)
  - [Pathway analyses with genekitr](#gene-set-testing-with-genekitr)
  - [CESR gene list comparison](#cesr-gene-list-comparison)
  - [RT-qPCR](#rt-qpcr)
  - [Interactome (Co-IP and TurboID)](#interactome-co-ip-and-turboid)
- [Reproducibility notes](#reproducibility-notes)
- [Publications](#publications)
- [License](#license)
- [Contact info](#contact)

---

## Repository structure

```
Nup211_analyses/
├── salmon_scripts/
│   ├── build_salmon_index.sh          # Decoy-aware index, ASM294v3
│   └── salmon_quant.sh                # Quantification, 16 h dataset
├── STAR_scripts/
│   └── 48h_star_alignment.sh          # Alignment, 48 h dataset
├── R_scripts/
│   ├── RNA-seq/                       # DESeq2, filtering, plots, GSEA, GO
│   ├── RT-qPCR/                       # Per-experiment relative quantification
│   └── Interactome/                   # Co-IP and TurboID mass spec analysis
├── data/
│   ├── rnaseq/                        # CESR reference gene lists
│   ├── rtqpcr/                        # Raw Cq exports and derived summary tables
│   └── mass_spec/                     # Protein identification tables
├── thesis_tables/                     # Dissertation appendix tables (.xlsx)
├── sessionInfo.txt                    # Package versions for all R scripts
├── .gitignore
├── LICENSE
└── README.md
```

---

## Data availability

### Sequencing data (NCBI GEO)

| Dataset                       | Timepoint | GEO Accession                                                             | Platform                |
| ----------------------------- | --------- | ------------------------------------------------------------------------- | ----------------------- |
| _nup211-so_ depletion RNA-seq | 48 h      | [GSE216106](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216106) | Illumina NovaSeq 6000   |
| _nup211-so_ depletion RNA-seq | 16 h      | [GSE342656](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE342656) | Illumina NovaSeq X Plus |

> Reviewer token for GSE342656 available upon request.

Raw FASTQ files are available through the linked Gene Expression Omnibus (GEO) records and their associated Sequence Read Archive (SRA) links. They are **not** stored in this repository.

### Processed RNA-seq files

The R scripts can be run without re-processing raw reads. The count matrices are hosted on GEO rather than in this repository. Download the following supplementary files before running the RNA-seq scripts:

| File                                | Source    | Description                           |
| ----------------------------------- | --------- | ------------------------------------- |
| nup211so_16h_counts.csv             | [GSE342656](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE342656) | Gene-level counts (16 h)              |
| DESeq2_results_protein_coding.csv   | [GSE342656](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE342656) | DESeq2 results, protein-coding subset |
| GSE216106_counts.csv.gz             | [GSE216106](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216106) | Gene-level counts (48 h)              |
| GSE216106_DEseq2Results.csv.gz      | [GSE216106](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216106) | Published 48 h DESeq2 results         |


> **Note on the 48 h dataset:** the 48 h script provided loads the counts table directly from GSE216106 rather than re-running `summarizeOverlaps` from `GenomicAlignments`. 


### Data included in this repository

| Example file path                                                    | Contents                                                                                    |
| -------------------------------------------------------------------- | ------------------------------------------------------------------------------------------- |
| data/rnaseq/cesr_induced.txt, cesr_repressed.txt                     | Core environmental stress response (CESR) gene lists                                        |
| data/rnaseq/cesr_induced_pombase.tsv, cesr_repressed_pombase.tsv     | The same lists mapped to PomBase systematic IDs                                             |
| data/rnaseq/cesr_induced.txt, cesr_repressed.txt                     | Core environmental stress response (CESR) gene lists from Chen et al. (2003), as gene names |
| data/rnaseq/cesr_induced_pombase.tsv, cesr_repressed_pombase.tsv     | The same lists mapped to PomBase systematic IDs, with gene products                         |
| data/rnaseq/Transcriptome_CESR_analysis_2x3.csv                      | Statistical analysis of CESR  genes at both time points                                     |
| data/rtqpcr/CT_values_*.xlsx                                         | Raw Ct values                                                                               |
| data/rtqpcr/shutoff_FC_(variable).csv, shutoff_mean_values_*.csv              | Fold-change and mean values produced by the RT-qPCR scripts                                 |
| data/rtqpcr/supplemental_*.csv                                       | Summary tables produced by the RT-qPCR scripts                                              |
| data/mass_spec/Nup211_IP_MS_Results.xlsx                             | Co-immunoprecipitation MS protein identifications                                           |
| data/mass_spec/TurboID_summary.xlsx                                  | TurboID proximity labeling MS protein identifications                                       |


---
## Thesis tables

Full-length versions of the dissertation appendix tables (gene set enrichment results, leading-edge gene lists, RT-qPCR summary statistics, mass spec hit lists, etc.) are found in [thesis_tables/](thesis_tables/). and are numbered to match the dissertation (`Table A10.` etc.).

---

## Requirements

**Command line**

- Salmon v1.10.3
- STAR v2.5.4b

**R (≥ 4.3)** -- Used R version 4.6.0 

See [sessionInfo.txt](sessionInfo.txt) for the complete list of package versions.

```r
# Prototype code to load/install all required packages in one shot 

cran <- c(
	"tidyverse", "magrittr", "scales", "generics", "matrixStats",
	"fstcore","cowplot", "patchwork", "ggrepel", "ggtext", "ggbreak",
	"ggExtra","ggsignif", "ggpubr", "ggvenn", "eulerr", "plotly", "gplots",
	"pheatmap","viridis", "RColorBrewer", "imager", "broom", "rstatix",
	"PoiClaClu","genekitr", "geneset", "flextable", "officer", "flexlsx",
	"openxlsx2","writexl", "readxl", "kableExtra", "knitr"
	)

bioc <- c(
	"DESeq2", "apeglm", "tximport", "EnhancedVolcano", "Glimma",
	"GenomicRanges", "GenomicAlignments", "GenomicFeatures", "txdbmaker",
	"Rsamtools", "Biostrings", "SummarizedExperiment", "BiocParallel",
	"AnnotationDbi", "AnnotationForge", "GEOquery", "DOSE"
	)

github <- c(pcr = "MahShaaban/pcr") # archived on CRAN

## Install anything missing

for (p in c("BiocManager", "remotes")) {
	if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
	}

new <- setdiff(cran, rownames(installed.packages()))
	if (length(new)) install.packages(new)

new <- setdiff(bioc, rownames(installed.packages()))
if (length(new)) BiocManager::install(new, update = FALSE, ask = FALSE)

new <- setdiff(names(github), rownames(installed.packages()))
for (p in new) remotes::install_github(github[[p]])

## Load

invisible(lapply(c(cran, bioc, names(github)), function(p) {
	suppressPackageStartupMessages(library(p, character.only = TRUE))
}))
```

---

## Analyses

### RNA-seq: 16 h depletion (Salmon)

[salmon_scripts/](salmon_scripts/) builds a decoy-aware index (k = 31) from `GCF_000002945.2_ASM294v3_rna.fna` (genome build ASM294v3) and quantifies transcript abundance using  `--validateMappings` and `-l A` (library type was inferred as ISR).

[R_scripts/RNA-seq/](R_scripts/RNA-seq/) reads the resulting count matrix into DESeq2 via `DESeqDataSetFromMatrix` using rounded integer counts. Protein-coding transcripts are selected using the PomBase `gene_type` annotation.

_Produces:_ Table A2 /  `nup211so_16h_counts.csv` file [GSE342656](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE342656) 

### RNA-seq: 48 h depletion (STAR)

[STAR_scripts/](STAR_scripts/) aligns with STAR v2.5.4b. Counting is done in R with `summarizeOverlaps` (GenomicAlignments) using `invertStrand`, against features built with `makeTxDbFromBiomart` using information from Ensembl Fungi.

*Produces:* Table A8 /  `GSE216106_DEseq2Results.csv.gz` file (GSE216106)

### Gene set testing with genekitr

Enrichment analysis uses [**genekitr**](https://www.genekitr.fun/), which performs over-representation analysis (`genORA`) and GSEA (`genGSEA`) on gene sets retrieved by its companion package [**geneset**](https://github.com/GangLiLab/geneset).

**Gene set source.** `geneset` does not bundle gene sets with the package. Each `getGO()` call downloads the gene set information from the `geneset` server, which is obtained from the official Gene Ontology database. The `geneset` backend server is updated on a **monthly basis**.

```r
# Get gene sets for S. pombe
sp_gs <- geneset::getGO(org = "Schizosaccharomyces pombe", ont = 'bp')
```

Gene sets for the dissertation analyses were initially downloaded on **2025-01-08**; subsequent reanalysis used freshly retrieved gene sets.

*Produces:* Figure 4.4 and Table A10.  


### CESR gene list comparison

Differentially expressed genes at both timepoints were compared against the core environmental stress response (CESR) gene lists defined by Chen et al. (2003). These genes are are commonly induced and repressed in response to multiple stress conditions.

Gene names in the published lists were mapped to PomBase systematic IDs before
comparison. Both the original and mapped versions are in [data/rnaseq/](data/rnaseq/). The script counts CESR-induced and CESR-repressed genes among the genes called up and down at 16 h and 48 h, and uses Fisher's exact test to assess whether the overlap exceeds what would be expected by chance.

*Produces:* Table A6

> Chen, D., Toone, W.M., Mata, J., Lyne, R., Burns, G., Kivinen, K., Brazma, A., Jones, N., and Bahler, J. (2003). Global transcriptional responses of fission yeast to environmental stress. _Mol Biol Cell_ **14**, 214-229. 10.1091/mbc.e02-08-0499. PMID: [**12529438**](https://pubmed.ncbi.nlm.nih.gov/12529438/)

### RT-qPCR

[R_scripts/RT-qPCR/](R_scripts/RT-qPCR/) contains several files:

|Script|Experiment|
|---|---|
|RTqPCR-16h.R|Validation of 16 h RNA-seq hits|
|RTqPCR-rescue.R|Rescue of the depletion phenotype|
|RTqPCR-655.R|CRISPR truncation strain 655 <!-- confirm wording -->|
|RTqPCR-863.R|CRISPR truncation strain 863 <!-- confirm wording -->|
|RTqPCR-leading-edge.R|Leading-edge gene panel|
|RTqPCR-leading-edge-pipeline_.R|Consolidation of target gene and sample condition combinations for the leading-edge panel|

Relative quantification uses the `pcr` package. The final dissertation gene panel is `chp1`, `mit1`, `swi6`, `for3`, `tea2`, `mus7`, `not1`.

_Produces:_ Tables A11 and A12

### Interactome (Co-IP and TurboID)

[R_scripts/interactome/](R_scripts/interactome/) processes mass spectrometry hit lists, cross-references them against PomBase, and calculates statistics like protein coverage and mean relative abundance.

*Produces:* Tables A13 and A14  

---

## Reproducibility notes

- **The two RNA-seq timepoints were processed on different pipelines.** 
	- The 16 h analysis used Salmon &rarr; `DESeqDataSetFromMatrix` with rounded counts and no length offsets. 
	- The 48 h analysis used STAR &rarr; `summarizeOverlaps`. 

- File paths in the scripts need to be adjusted to point to where data files live on your system.

---
## Publications

Publications associated with this work include: 
- Kamel, D., Sookdeo, A., Ikenouchi, A., and Zhong, H. (2024). Fission yeast essential nuclear pore protein Nup211 regulates the expression of genes involved in cytokinesis. _PLoS One_ **19**, e0312095. 10.1371/journal.pone.0312095. PMID: [**39666777**](https://pubmed.ncbi.nlm.nih.gov/39666777/)
- Kamel D. (2025). *Functional Analysis of the Fission Yeast Nuclear Basket Protein Nup211.* PhD dissertation, The Graduate Center, City University of New York. <!-- CUNY Academic Works / ProQuest link -->

---

## License

MIT - see [LICENSE](LICENSE).

---

## Contact

Domenick Kamel - domenick.kamel53@gc.cuny.edu

