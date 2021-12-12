# A comprehensive atlas of fetal splicing patterns in the brain of adult myotonic dystrophy type 1 patients

This repository stores R scripts and annotation files to reproduce the analysis of the publication "A comprehensive atlas of fetal splicing patterns in the brain of adult myotonic dystrophy type 1 patients". 

A preprint is available on bioRxiv: https://www.biorxiv.org/content/10.1101/2021.10.01.462715v1

## Prerequisites

Software:
* Git: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git
* Conda: https://docs.conda.io/en/latest/miniconda.html 

Data:
* Gene read counts for GTEx V8 can be downloaded from the public portal: https://gtexportal.org/home/datasets.
* Acquire raw RNA-Seq data for the following datasets:
  * [Otero et al. (2021)](https://www.cell.com/cell-reports/fulltext/S2211-1247(20)31623-5): [GSE157428](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157428) (Navigate to the SRA Run Selector)
  * BrainSpan Atlas of the Developing and Adult Human Brain: [phs000731.v2.p1](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000731.v2.p1)*
  * Genotype-Tissue Expression (GTEx) project: [phs000424.v8.p2](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000424.v8.p2)*
* We recommend to use this pipeline for data preprocessing: https://github.com/MDegener/RNAseq-pipeline.

*protected data that requires authorized access by dbGaP

## Getting started

1. Checkout this git repository: `git clone https://github.com/pacthoen/BrainDM1.git`

2. Navigate to the repository: `cd BrainDM1`

3. Create and activate conda environment with all R package dependencies
```
conda create -p r-env -c r -c bioconda -c conda-forge git rstudio r-base r-tidyr r-dplyr r-stringr r-purrr r-ggplot2 r-corrplot  bioconductor-rtracklayer bioconductor-biomaRt bioconductor-edgeR r-matrixTests r-psych r-data.table r-dunn.test r-cairo r-statmod r-gtools r-ppcor r-argparse r-r.utils r-venndiagram r-cowplot r-gghighlight r-ggrepel  r-ggpubr

conda activate r-env
```

4. Unzip annotation files: `gzip -d lib/SE.hg38.annotated.gff3.gz lib/gencode.v26.annotation.collapsed.gtf.gz`

5. Run the following R scripts in this order:
   1. [`scripts/prepareSampleMetadata.R`](https://github.com/pacthoen/BrainDM1/blob/main/scripts/prepareSampleMetadata.R): Combines metadata for all selected samples into one table
   2. [`scripts/prepareSummarizedData.R`](https://github.com/pacthoen/BrainDM1/blob/main/scripts/prepareSummarizedData.R): Creates a matrix of exon inclusion data for all selected samples 
   3. [`scripts/comparePSI.R`](https://github.com/pacthoen/BrainDM1/blob/main/scripts/comparePSI.R): Performs a group comparison of exon inclusion for all exon-skipping events

6. Now you can run any script in the `scripts/figure` or `scripts/tables` directory to reproduce the unedited content of the publication

## Overview of directories and content

| Directory     | Content                                                                                       |
|---------------|-----------------------------------------------------------------------------------------------|
| lib           |  Misc files for annotation and selection of samples                                           |
| results       |  Unedited output of the analysis scripts that are contained in the `/scripts` directory       |
| scripts       |  Scripts to create summarized data and to reproduce all figures and tables of the publication |              
