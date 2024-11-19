# Peak Analysis Worksop

| Audience | Computational Skills | Prerequisites | Duration |
:----------|:----------|:----------|:----------|
| Biologists | Intermediate | None | Introduction to R |


This repository has teaching materials for a hands-on **Peak Analysis** workshop. This workshop will use the R statistical programming environment to evaluate files generated from peak calling of ChIP-seq (and related approaches i.e. CUT&RUN and ATAC-seq) data. We will provide participants with a suite of tools and a basic workflow beginning with quality metrics through to annotation and visualization. This workshop will introduce participants to:

* File formats for peak data
* Approaches to check peak quality and reproducibility across replicates
* Peak annotation methods and tools for visualization
* Differential peak enrichment analysis and functional analysis

Working knowledge of R is required or completion of the [Introduction to R workshop](https://hbctraining.github.io/Intro-to-R/). 

**Note for Trainers:** Please note that the schedule linked below assumes that learners will spend between 3-4 hours on reading through, and completing exercises from selected lessons between classes. The online component of the workshop focuses on more exercises and discussion/Q & A.

> These materials were developed for a trainer-led workshop, but are also amenable to self-guided learning.

### Learning Objectives - UPDATE
 * ...
 * 


### Lessons
* [Workshop schedule (trainer-led learning)](schedule/README.md)
* [Self-learning]()

  
### Dataset - ADD LINKS

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
ChIPseeker
DiffBind
clusterProfiler
TxDb.Mmusculus.UCSC.mm10.knownGene
IRanges
GenomicRanges

```

> **NOTE:** The library used for the annotations associated with genes (here we are using `TxDb.Mmusculus.UCSC.mm10.knownGene` and `EnsDb.Hsapiens.v75`) will change based on organism (e.g. if studying mouse, would need to install and load `TxDb.Mmusculus.UCSC.mm10.knownGene`). The list of different organism packages are given [here](https://github.com/hbctraining/Training-modules/raw/master/DGE-functional-analysis/img/available_annotations.png).

(3) Finally, please check that all the packages were installed successfully by **loading them one at a time** using the `library()` function.  

```r
library(tidyverse)
library(ChIPseeker)
library(DiffBind)
library(clusterProfiler)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(EnsDb.Hsapiens.v75)
```

(4) Once all packages have been loaded, run sessionInfo().  

```r
sessionInfo()
```
