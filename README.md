# Transcriptomic Profiling of Inflammatory vs Noninflammatory Aortic Aneurysms

## Overview
This repository contains the complete analysis workflow, results, and figures for my MSc thesis titled:

**“Transcriptomic Profiling Identifies Differential Gene Expression Signatures Between Inflammatory and Noninflammatory Aortic Aneurysms.”**

The study uses bulk RNA-sequencing data and bioinformatics analysis to identify molecular differences between inflammatory aortic aneurysms (IAA) and noninflammatory aortic aneurysms (NIAA), with a focus on immune activation, extracellular matrix remodeling, and vascular smooth muscle cell dysfunction.

---

## Dataset
- **Source:** NCBI Sequence Read Archive (SRA)
- **Accession:** SRP432690
- **Samples:**  
  - 10 Inflammatory aortic aneurysm (IAA) samples  
  - 10 Noninflammatory aortic aneurysm (NIAA) samples
- **Data type:** Human bulk RNA-seq (paired-end)

Only processed results and summary files are included in this repository. Raw FASTQ/BAM files are not uploaded due to size limitations.

---

## Bioinformatics Workflow
The transcriptomic analysis was performed using the following pipeline:

1. **Quality Control:** FastQC  
2. **Read Trimming:** Trimmomatic  
3. **Alignment:** HISAT2 (GRCh38 reference genome)  
4. **File Processing:** SAM → BAM using SAMtools  
5. **Transcript Assembly & Quantification:** StringTie  
6. **Differential Expression Analysis:** Ballgown (FPKM-based)  
7. **Filtering:** Variance-based filtering (row variance > 1)  
8. **Functional Enrichment:** DAVID (GO and KEGG)  
9. **Visualization:** Boxplots, Volcano plots, Heatmaps

---

## Key Results
- **Genes analyzed after filtering:** 1,165  
- **Differentially expressed genes (DEGs):** 147  
  - Upregulated in IAA: 126  
  - Downregulated in IAA: 21  

### Key Upregulated Genes
- **IL6, CXCL8, TNFAIP3, MMP9, SPP1**  
Indicating immune activation, cytokine signaling, and ECM degradation in inflammatory aneurysms.

### Key Downregulated Genes
- **ACTG2, CNN1, MYH11**  
Associated with vascular smooth muscle cell structure and contractility.

---

## Functional Enrichment Summary
DAVID enrichment analysis revealed significant involvement of immune-related pathways, including:
- TNF signaling pathway  
- IL-17 signaling pathway  
- NF-κB signaling pathway  
- Phagosome pathway  
- Cytokine–cytokine receptor interaction  

These pathways support an inflammation-driven mechanism in IAA compared to NIAA.

---

## Repository Structure

- `scripts/` – R scripts used for transcriptomic and differential expression analysis
- `results/` – Processed results, DEG tables, pathway enrichment outputs, and plots
- `docs/` – Final MSc thesis document


## Repository Contents
- `Ballgown_Analysis.R` – Differential expression analysis script  
- `DEG_summary_table.csv` – Summary of DEGs with statistics  
- `top_30_kegg_pathway_david.ods` – DAVID KEGG enrichment results  
- `*.png` – Visualization figures (boxplots, heatmaps, volcano plots)  
- `MSc thesis final doc.pdf` – Full thesis document  
- `README.md` – Project overview and documentation  

---

## How to Run the Analysis
The analysis was performed in **R** using the following key packages:
- `ballgown`
- `dplyr`
- `matrixStats`
- `ggplot2`
- `ggrepel`
- `pheatmap`
- `readr`

The main analysis script is provided in `Ballgown_Analysis.R`.

---


## Author
**Nighitha T N**  
MSc Bioinformatics (NGS Data Analytics)  
University of Kerala  

---

## License
This project is shared for academic and educational purposes. Please cite appropriately if used.
