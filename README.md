# Analyses of peak call data generated by high-throughput sequencing

| Audience | Computational Skills | Prerequisites | Duration |
:----------|:----------|:----------|:----------|
| Biologists | Intermediate | None | Introduction to R |

> NOTE: **This workshop is currently under development.**
> Materials in thsi repo are not currently maintained


This workshop will focus on using the R statistical programming environment to evaluate files generated from peak calling of ChIP-seq (and related approaches i.e. CUT&RUN and ATAC-seq) data. We describe the different file formats encountered when working with peaks and use various Bioconductor packages to look at concordance between replicates, peak quality and annotate regions using nearest gene approaches. We demonstrate the use of DiffBind to evaluate changes in binding patterns between groups of samples, and how to explore genomic regions of interest as tracks in the Integrated Genome Viewer (IGV). We also briefly touch on various tools for motif-based sequence analyses.


* File formats for ChIP-seq
* Peak concordance between replicates
* ChIPQC?
* Peak annotation and visualization
* Differential enrichment analysis
* Peak visualization using a genome viewer (IGV)


### Lessons
**[Click here]() for links to lessons and the suggested schedule**

### Dataset

### Installation Requirements

Download the most recent versions of R and RStudio for your laptop:

 - [R](http://lib.stat.cmu.edu/R/CRAN/) 
 - [RStudio](https://www.rstudio.com/products/rstudio/download/#download)
 
> **NOTE**: When installing the following packages, if you are asked to select (a/s/n) or (y/n), please select “a” or "y" as applicable.

(1) Install the below packages on your laptop from CRAN. You DO NOT have to go to the CRAN webpage; you can use the following function to install them:


```r
install.packages("BiocManager")
install.packages("tidyverse")
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
AnnotationDbi
TxDb.Hsapiens.UCSC.hg19.knownGene
EnsDb.Hsapiens.v75
org.Hs.eg.db
```

> **NOTE:** The library used for the annotations associated with genes (here we are using `TxDb.Hsapiens.UCSC.hg19.knownGene` and `EnsDb.Hsapiens.v75`) will change based on organism (e.g. if studying mouse, would need to install and load `TxDb.Mmusculus.UCSC.mm10.knownGene`). The list of different organism packages are given [here](https://github.com/hbctraining/Training-modules/raw/master/DGE-functional-analysis/img/available_annotations.png).

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
