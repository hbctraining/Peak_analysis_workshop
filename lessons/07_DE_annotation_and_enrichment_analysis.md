---
title: "Annotation and differential enrichment analysis"
author: "Heather Wick, Upendra Bhattarai, Meeta Mistry"
date: "Aug 16th, 2024"
---

Contributors: Heather Wick, Upendra Bhattarai, Meeta Mistry

Approximate time: 

## Learning Objectives

* Annotate and visualize the differentially bound regions.
* Perform over-representation analysis using clusterProfiler on the significant genes from DiffBind.
* Discuss the Functional analysis approaches and the biological insights from the analysis.



## Overview

Functional enrichment analysis will determine whether some functions are enriched in the differentially bound sites. Here, we will map the differentially bound sites to a functional annotation database and visualize the enrichment. Functional enrichment analysis gives insights on the collective function of the group of genes rather than individuals.

This is the next step after we get a list of differentially bound sites in the comparison.
Functional enrihment analysis generally follows following three decision steps.

1. Use of different gene identifiers and gene descriptions for functional annotation (Entrezid, Uniprot, KEGG, etc.)
2. Over-Representation Analysis (ORA) vs Gene Set Enrichment Analysis.
3. Using R or Web-based tools.
4. Interpreting the results.

In this session, we will peform Over-representation analysis (ORA) using an R Bioconductor package called `clusterProfiler`.

# Over-representation analysis

Over-representation analysis is used to determine whether the biological functions or pathways represented in the list of interesting genes occur more than expected by chance (over-represented) compared to the biological functions or pathway in the complete list of genes. Most genes in the genome have some pre-existing annotation associated with it which has been compiled through a combination of manual curation and computational algorithms. There are a number of existing databases which define genes using a controlled vocabulary and then categorize genes into groups (gene sets) based on shared function, or involvement in a pathway, or presence in a specific cellular location etc. A very commonly used gene annotataion resource is the Gene Ontology (GO) database, and is what we will use in our workflow.

## Hypergeometric test
The statistical test that will determine whether something is actually over-represented is the Hypergeometric test.

Hypergeometric distribution is a probability distribution that describes the probability of some number of genes (k) being associated with "Functional category 1", for all genes in our gene list (n), compared to drawing some number of genes (M) associated with "Functional category 1" from a population of all of the genes in entire genome (N).

The p-value can be calculated as:

$$ P(X = k) = \frac{\binom{K}{k} \binom{N - K}{n - k}}{\binom{N}{n}} $$

This test will result in an adjusted p-value (after multiple test correction) for each category tested.

# Running ORA with clusterProfiler

Now lets take our significantly differentially bound sites and their nearest gene annotation to determine if there are any GO terms over-represented in the list of our gene of interest. 

First open a R script in our RStudio and load the GRange object saved with the result from DiffBind analysis.


```{r}
#Libraries to load if not already loaded
library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(clusterProfiler)

res_all <- readRDS("res_all.rds")

```

Annotate the with ChipSeeker
This GRange object contains the result output from the Diffbind analysis. It has the information about the sites those were analysed with genomic coordinates and the other statistics like Fold change, pvalue and FDR for each of those sites.
Before we proceed with the functional analysis, we need to annotate the site loci with the nearest gene name. To do that we will use chipSeeker.

```{r}
annot_res_all <- annotatePeak(res_all, tssRegion = c(-3000, 3000), TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, annoDb = "org.Hs.eg.db")
```

Note: You can plot the annotation using chipSeeker to see the distribution of genomic features in our result file.

```{r}
plotAnnoPie(annot_res_all)
plotAnnoBar(annot_res_all)
upsetplot(annot_res_all)
plotDistToTSS(annot_res_all)
```

To run ORA analysis we need a background dataset and a query dataset. In our case, we will use all the sites in our analysis as a background and significantly upregulated sites in cKO vs WT as a query set for the hypergeometric test.


Lets convert the annotated GRange object to a dataframe.
```{r}
annot_res_all_df <- as.data.frame(annot_res_all)
```

Create a background dataset for the hypergeometric testing.

```{r}
background_set <- as.character(annot_res_all_df$geneID)
```

Extract gene list for the significantly upregulated genes in cKO vs WT to prepare a query set for ORA.

```{r}
sigUP <- dplyr::filter(annot_res_all_df, FDR < 0.1, Fold > 0)
sigUp_genes <- as.character(sigUp$geneId)
```
Now we can perform ORA with Gene Ontology (GO) dataset as follows.
```{r}
go_ORA_Up <- enrichGO(gene = sigUp_genes,
                      universe = background_set,
                      keyType = "ENTREZID",
                      OrgDb = org.Hs.eg.db,
                      ont = "ALL",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      readable = TRUE)
```

Note: Note 1: The different organisms with annotation databases available to use with for the OrgDb argument can be found here.

Note 2: The keyType argument may be coded as keytype in different versions of clusterProfiler.

Note 3: The ont argument can accept either "BP" (Biological Process), "MF" (Molecular Function), and "CC" (Cellular Component) subontologies, or "ALL" for all three.

Lets save the ORA results

```{r}
go_ORA_Up_df <- data.frame(go_ORA_Up)
write.csv(go_ORA_Up_df, "results/GO_ORA_clusterProfiler_cKO_vs_WT_Upregulated.csv
```

# Exploring results from ORA analysis.

Let's take a look at what terms are identified as over-represented in the genes up-regulated in cold conditions.
```{r}
View(go_ORA_Up_df)
```
In the first few columns we see the GO identifier and the descriptive term name. In the next two columns that follow, we observe GeneRatio and BgRatio. These values allows us to compare the overlaps to the background.

BgRatio: M/N
The total number of genes in the GO term gene set (M), divided by the total number of genes in universe (N)

GeneRatio: k/n
The total number of genes in our sig DE gene set which overlap with the GO term gene set (k), divided by the total number of genes in our sig DE gene set that overlap with the universe gene set (n).

Other columns of interest are the p.adjust column (by which results are ordered by default), and the geneID column which lists the gene symbols of the overlapping genes.

<p align="center">
<img src="../img/Go_ORA_up.png"  width="600">
</p>

Exercise:
1. Carryout the ORA for the down regulated sites in cKO vs WT result.
2. Do you find any significantly over represented terms?

# Plotting the ORA results

There are multiple options to plot the ORA results through clusterProfiler. Lets explore few of them:

## Dotplot
It shows statistics associated with a user-selected top number (default=10) of significant terms. The color of the dots represent the p-adjusted values for these terms, and size of the dots corresponds to the total count of sig DE genes annotated with the GO term (count). This plot displays the enriched GO terms ordered by gene ratio, not p-adjusted value.

<p align="center">
<img src="../img/ORA_dotplot_go.png"  width="600">
</p>


~~~~~~~~~~~~~~~~~~~`
- Annotate DE regions with target genes.
- IGV validation of differentially enriched regions
- Functional analysis
    [old material](https://github.com/hbctraining/Peak_analysis_workshop/blob/main/lessons/OLD_ChipSeeker_analysis.md#functional-enrichment-r-based-tools)
    [webased functional analysis](https://github.com/hbctraining/Peak_analysis_workshop/blob/main/lessons/OLD_web_based_functional_analysis.md)
***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
