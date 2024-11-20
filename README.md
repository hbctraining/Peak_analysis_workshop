# Peak Analysis Worksop

| Audience | Computational Skills | Prerequisites | Duration |
:----------|:----------|:----------|:----------|
| Biologists | Intermediate | None | Introduction to R |


### Learning Objectives 
 * Understand peak data and different file formats generated from peak calling algorithms 
 * Assess various metrics used to assess the quality of peak calls
 * Compare peak calls across samples within a dataset
 * Create visualizations to evaluat peak annotations
 * Evaluate differentially enriched regions between two sample groups 

> These materials were developed for a trainer-led workshop, but are also amenable to self-guided learning.


### Lessons
* [Workshop schedule (trainer-led learning)](schedule/README.md)
* [Self-learning]()

### Description
This repository has teaching materials for a hands-on **Peak Analysis** workshop. This workshop will use the R statistical programming environment to evaluate files generated from peak calling of ChIP-seq (and related approaches i.e. CUT&RUN and ATAC-seq) data. We will provide participants with a suite of tools and a basic workflow beginning with quality metrics through to annotation and visualization. This workshop will introduce participants to:

* File formats for peak data
* Approaches to check peak quality and reproducibility across replicates
* Peak annotation methods and tools for visualization
* Differential peak enrichment analysis and functional analysis

> Working knowledge of R is required or completion of the [Introduction to R workshop](https://hbctraining.github.io/Intro-to-R/). 

**Note for Trainers:** Please note that the schedule linked below assumes that learners will spend between 3-4 hours on reading through, and completing exercises from selected lessons between classes. The online component of the workshop focuses on more exercises and discussion/Q & A.


### Dataset 
The R project for this workshop can be [downloaded with this link](https://www.dropbox.com/scl/fi/s9mxwd7ttqgjt040m6bm2/Peak_analysis.zip?rlkey=ceqbv4pyx59jxsoa0xoh9l6kb&st=q7rlclil&dl=1).

### Installation Requirements 

Download the most recent versions of R and RStudio for your laptop:

 - [R](http://lib.stat.cmu.edu/R/CRAN/) 
 - [RStudio](https://www.rstudio.com/products/rstudio/download/#download)
 
> **NOTE**: When installing the following packages, if you are asked to select (a/s/n) or (y/n), please select “a” or "y" as applicable.

(1) Install the below packages on your laptop from CRAN. You DO NOT have to go to the CRAN webpage; you can use the following function to install them:


```r
install.packages("BiocManager")
install.packages("tidyverse")
install.packages("pheatmap")
install.packages("UpSetR")
```

**Note that these package names are case sensitive!**


(2) Install the below packages from Bioconductor. Load BiocManager, then run BiocManager's `install()` function 7 times for the 7 packages:

```r
library(BiocManager)
install("insert_first_package_name_in_quotations")
install("insert_second_package_name_in_quotations")
& so on ...
```

Note that these package names are case sensitive!

```r
BiocManager::install("ChIPseeker")
BiocManager::install("DiffBind")
BiocManager::install("clusterProfiler")
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
BiocManager::install("IRanges")
BiocManager::install("GenomicRanges")

```

> **NOTE:** The library used for the annotations associated with genes (here we are using `TxDb.Mmusculus.UCSC.mm10.knownGene`) will change based on organism. The list of different organism packages are given [here](https://github.com/hbctraining/Training-modules/raw/master/DGE-functional-analysis/img/available_annotations.png).

(3) Finally, please check that all the packages were installed successfully by **loading them one at a time** using the `library()` function.  

```r
library(tidyverse)
library(pheatmap)
library(UpSetR)
library(ChIPseeker)
library(DiffBind)
library(clusterProfiler)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(IRanges)
library(GenomicRanges)
```

(4) Once all packages have been loaded, run sessionInfo().  

```r
sessionInfo()
```
