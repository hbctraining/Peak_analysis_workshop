---
title: "Peak annotation and visualization"
author: "Upendra Bhattarai, Meeta Mistry"
date: "Aug 16th, 2024"
---

Contributors: Heather Wick, Upendra Bhattarai, Meeta Mistry

Approximate time: 45 minutes

## Learning Objectives

* Annotate peaks with genomic features using Chipseeker
* Visualize annotations and compare peak coverage between experimental groups


## Peak annotation 
Understanding the biological questions addressed by ChIP-seq experiments begins with annotating the genomic regions we have identifed as peaks with genomic context. In order to interpret these binding regions, a number of different peak annotation tools exist. Some examples include [Homer](http://homer.ucsd.edu/homer/motif/), [GREAT](https://pmc.ncbi.nlm.nih.gov/articles/PMC4840234/) (a web-based tool), [ChIPpeakAnno](https://www.bioconductor.org/packages/devel/bioc/vignettes/ChIPpeakAnno/inst/doc/ChIPpeakAnno.html) and [ChIPseeker](https://www.bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html). Because many cis-regulatory elements are close to transcription start sites of their targets, it is common to associate each peak to its nearest gene, either upstream or downstream. Annotation tools  **apply methods to calculate the nearest TSS to the given genomic coordinates and annotates the peak with that gene**. However, problems exist in regions where multiple genes are located in close proximity. Different tools address this complex issue using different approaches and this can result in varying results.

<p align="center">
<img src="../img/nearest_gene_image.png"  width="800">
</p>

_Image source: Welch R.P. et al, Nucleic Acids Research, 2014 [doi: 10.1093/nar/gku463ChIP](https://www.researchgate.net/publication/262812725_ChIP-Enrich_Gene_set_enrichment_testing_for_ChIP-seq_data)_


## Annotating peaks 
In this workshop we will use an R Bioconductor package called **[ChIPseeker](https://bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html) to annotate peaks, visualize features, and compare profiles**. Some features of ChIPseeker include:

* Comparing results in batch; across replicates or between experimental groups
* Perform functional annotations, and infer cooperative regulation
* Support querying the GEO database to compare experimental datasets with publicly available ChIP-seq data and offers statistical testing for significant overlaps among datasets.

### Setting up
Let's open up a new script file and call it `peak_annotation.R`. Add a title to the script and as usual we will begin with loading required libraries:

```
## Peak Annotation using ChIPseeker

# Load libraries
library(ChIPseeker)
library(clusterProfiler)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(tidyverse)
```

> **NOTE:** The `readPeakFile()` function allows us to load peaks from a BED file. It will convert the BED file into GRanges object. **However, we already have all of our peak data loaded in our environment as GRanges object from the previous lesson.**

First, we will need to assign the **annotation database** to a variable. Bioconductor provides the annotation database `TxDb` for commonly used genome versions to use for annotation in ChIPseeker, e.g. `TxDb.Mmusculus.UCSC.mm10.knownGene`, `TxDb.Hsapiens.UCSC.hg38.knownGene`, but users can also prepare prepare their own database object using UCSC Genome Bioinformatics and BioMart database in R with `makeTxDbFromBiomart()` and `makeTxDbFromUCSC.TxDb()` functions.

Annotation of the peaks to the nearest gene and for various genomic characteristics is performed by `annotatePeak()` function. By default TSS region is defined as -3kb to +3kb, however users can define this region. The result of the annotation comes in csAnno (a special format for ChIP-seq annotation). This can be converted to GRanges with `as.GRanges()` format and to data frame with `as.data.frame()` function.

Lets annotate our peakset:

```{r}
# Set the annotation database
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Annotate Peaks 
annot_WT1 <- annotatePeak(WT_H3K27ac_ChIPseq_REP1, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Mm.eg.db")

# View the result
annot_WT1@anno %>% 
  data.frame() %>% View()
```

> Note: The parameter `annoDb` is optional, if provided extra information about the annotation including gene symbol, genename, ensembl/entrezid will be added in extra columns.


<p align="center">
<img src="../img/anno_resullt_table.png"  width="800">
</p>


The annotation output file retains all the information from the peak file and also reports the annotation information. The nearest genes, their position and strand information along with the distance from peak to the TSS of its nearest gene is also reported.  For annotating genomic regions, `annotatePeak()` will not only give the gene information but also reports detailed information when genomic region is Exon or Intron. For example, ‘Exon (ENSMUST00000144339.1/ENSMUST00000144339.1, exon 1 of 3’, means that the peak overlaps with the 1st of 3 exons that transcript possesses.

To view a genomic annotation summary, we can just type out the variable name. We see that a fairly large proportion of our peaks for WT replicate 1 are in promotor regions, but there are also equallly high percentages in intronic and distal intergenic regions. 

```{r}
annot_WT1
```

```{r, output}
Annotated peaks generated by ChIPseeker
100559/100570  peaks were annotated
Genomic Annotation Summary:
              Feature   Frequency
9    Promoter (<=1kb) 20.81862389
10   Promoter (1-2kb)  3.84848696
11   Promoter (2-3kb)  3.17226703
4              5' UTR  0.20982707
3              3' UTR  1.83374934
1            1st Exon  1.59806681
7          Other Exon  2.97437325
2          1st Intron 11.96809833
8        Other Intron 25.65061307
6  Downstream (<=300)  0.07856084
5   Distal Intergenic 27.84733341
```

## Annotation visualization
There are several functions for visualization provided by ChIPseeker package to effectively visualize annotation of various genomic features. We will go through some of these options below.

### Piechart
This function provides a visualization of the suammry we generated earlier, in the form of a pie chart.

```{r}
# Piechart
plotAnnoPie(annot_WT1)
```

<p align="center">
<img src="../img/annot.pie1.png"  width="500">
</p>

### Barplot
We can plot the exact same data, but using a stacked barplot instead.

```{r}
# Barplot
plotAnnoBar(annot_WT1)
```

<p align="center">
<img src="../img/barplot_annot_wt1.png"  width="550">
</p>


### UpsetR

Annotation overlaps can be visualized by upsetR plot. Here, we use a function from ChIPseeker that grabs the required data and formats it to be compatible with UpSetR and draws the plot. With this plot we can see that there are many peaks that contain more than one annotation. We can observe the counts for various combinations of annotations.

```{r}
upsetplot(annot_WT1)
```

<p align="center">
<img src="../img/upset_anno_wt1.png"  width="600">
</p>


### Distribution of TF-binding loci with respect to TSS

The distance between the peak and the TSS of the nearest gene is also reported in the annotation output and can be visualzed with a barplot.

```{r}
plotDistToTSS(annot_WT1)
```

<p align="center">
<img src="../img/dist_tss_annot_wt1.png"  width="550">
</p>


## Visualizing multiple samples
These are all great ways to visualize the information from our annotation table, however so far we have only done this for a single sample. It would be very helpful to create similar plots after collating annotations across all samples in our dataset. In this way, we can **assess consistencies across replicates within a group and compare samples between groups.**

In order to do this, we will first combine the GRanges objects for each samples into a list and then annotate each sample using a for loop:

```{r}

# Create a list of GRanges objects
samples_list <- list(
  WT1 = WT_H3K27ac_ChIPseq_REP1,
  WT2 = WT_H3K27ac_ChIPseq_REP2,
  WT3 = WT_H3K27ac_ChIPseq_REP3,
  cKO1 = cKO_H3K27ac_ChIPseq_REP1,
  cKO2 = cKO_H3K27ac_ChIPseq_REP2,
  cKO3 = cKO_H3K27ac_ChIPseq_REP3)

# Annotate each sample
for(s in names(samples_list)){
  annot <- annotatePeak(samples_list[[s]], tssRegion=c(-3000, 3000),
                        TxDb=txdb, annoDb="org.Mm.eg.db")
  assign(paste0("annot_", s), annot)
}

```


### Chip peak annotation Feature distribution all samples
Annotation feature distribution of all the samples can also be plotted together, making a list of annotation.


plotAnnoBar(annot_list)
```
<p align="center">
<img src="../img/Chipseq_annot_all.png"  width="600">
</p>


### Heatmap profile of ChIP binding to TSS regions
To plot the profile of peaks binding to TSS region, we need to prepare the TSS regions. These are the flanking sequence of the TSS sites, here we are settting up the flanking regions of 2000bp upstream and 2000bp downstream. Then align the peaks that are mapping to these regions, and generate the tagMatrix. This tagMatrix can be visualized as a heatmap with the function `tagHeatmap()`.

```{r}
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
promoter <- getPromoters(TxDb = txdb, upstream = 2000, downstream = 2000)
tagMatrix_ckoR1 <- getTagMatrix(ckoR1, windows = promoter)
tagHeatmap(tagMatrix_ckoR1)
```

<p align="center">
<img src="../img/Chipseq_tagheatmap1.png"  width="600">
</p>


> Note: Tag matrices can also be created for other genomic regions and visualized as heatmaps.

### Average profile of ChIP peaks binding to TSS region

```{r}
plotAvgProf(tagMatrix_ckoR1, xlim=c(-2000, 2000), 
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
```

<p align="center">
<img src="../img/Chipseq_avg.profile1.png"  width="600">
</p>

> Note: Confidence intervals, estimated using bootstrap methods, are supported for characterizing binding profiles.



### Functional enrichment analysis

As we have the nearest gene annotation for the peaks, we can perform functional enrichment analysis to investigate the predominant biological theme among the genes. Depending on the biological questions we can use, clusterProfiler for GO, KEGG enrichment analysis, DOSE to investigate Disease ontology (DO), or ReactomePA for reactome pathway.

In this session, we are using Reactome database. 

```{r}
library(ReactomePA) #load the library if you have not already.
reac_ckoR1 <- enrichPathway(as.data.frame(annot.ckoR1)$geneId, organism = "mouse")
```

**Dotplot**

Dotplot is an easy way to visualize the enriched categories:

```{r}
dotplot(reac_ckoR1)
```

<p align="center">
<img src="../img/reac_cko1.png"  width="600">
</p>


## Exercise:

Carryout similar visualization and annotation for a replicate of WT sample `WT_Rep1`.
Work along the codes below:

Use `readPeakFile()` function to import the peakset as GRange object.

```{r}
wtR1 <- readPeakFile("data/macs2/narrowPeak/WT_H3K27ac_ChIPseq_REP1_peaks.narrowPeak")

```

### Peak coverage plot
Use `covplot()` to visualize the coverage of the peaks in the genome.
You can look at a different chromosome or the whole genome by removing `chrs=` options. 

```{r}
covplot(wtR1, weightCol = "V5", chrs = c("chr1", "chr2"))
```

Image below is for chromosome 1 and 2.

<p align="center">
<img src="../img/chipseq_covplot4.png"  width="600">
</p>

### Heatmap profile of ChIP binding to TSS regions

Generate the tagMatrix and plot the heatmap.


```{r}
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
promoter <- getPromoters(TxDb = txdb, upstream = 2000, downstream = 2000)
tagMatrix_wtR1 <- getTagMatrix(wtR1, windows = promoter)
tagHeatmap(tagMatrix_wtR1)
```

<p align="center">
<img src="../img/Chipseq_tagheatmap4.png"  width="600">
</p>


### Average profile of ChIP peaks binding to TSS region

```{r}
plotAvgProf(tagMatrix_wtR1, xlim=c(-2000, 2000), 
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
```

<p align="center">
<img src="../img/Chipseq_avg.profile4.png"  width="600">
</p>


### Peak annotation

Annotate the peakset using `annotatePeak()` function

```{r}
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
annot.wtR1 = annotatePeak(wtR1, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Mm.eg.db")
```

Genomic annotation summary:

```{r}
annot.wtR1
```

```{r, output}
Annotated peaks generated by ChIPseeker
100559/100570  peaks were annotated
Genomic Annotation Summary:
 
Feature                     Frequency
<fctr>                          <dbl>
9	Promoter (<=1kb)	       20.81762945		
10	Promoter (1-2kb)	       3.84749252		
11	Promoter (2-3kb)	       3.17226703		
4	5' UTR	                  0.20883263		
3	3' UTR	                  1.83275490		
1	1st Exon	                1.59806681		
7	Other Exon	              2.97337881		
2	1st Intron	             11.96909277		
8	Other Intron	           25.65459084		
6	Downstream (<=300)	      0.07856084
5	Distal Intergenic	      27.84733341					
```

### Annotation visualization
 

Use Pie, Bar, and upsetR plots to visualize the genomic annotation.

#### Piechart


```{r}
plotAnnoPie(annot.wtR1)
```

<p align="center">
<img src="../img/annot.pie4.png"  width="600">
</p>

#### Barplot

```{r}
plotAnnoBar(annot.wtR1)
```

<p align="center">
<img src="../img/Chipseq_annotbar4.png"  width="600">
</p>


#### UpsetR


```{r}
upsetplot(annot.wtR1)
```

<p align="center">
<img src="../img/Chipseq_upset4.png"  width="600">
</p>

### Distribution of TF-binding loci with respect to TSS

The distance between peak and the TSS of the nearest gene is reported in the annotation output and can be visualzed with barplot.

```{r}
plotDistToTSS(annot.wtR1)
```

<p align="center">
<img src="../img/Chipseq_tss4.png"  width="600">
</p>


### Functional enrichment analysis


Let's perform the enrichment analysis as before 

```{r}
reac_wtR1 <- enrichPathway(as.data.frame(annot.wtR1)$geneId, organism = "mouse")
```

**Dotplot**

Dotplot is an easy way to visualize the enriched categories:

```{r}
dotplot(reac_ckoR1)
```

<p align="center">
<img src="../img/reac_wt1.png"  width="600">
</p>

---


## Dataset comparison

We need the tagMatrix and annotations from other samples to carryout the overall comparison and visualization.

#### Sample CKO_Rep2

Importing peak sets
```{r}
ckoR2 <- readPeakFile("data/macs2/narrowPeak/cKO_H3K27ac_ChIPseq_REP2_peaks.narrowPeak"
```
Generate TagMatrix:

```{r}
promoter <- getPromoters(TxDb = txdb, upstream = 2000, downstream = 2000)
tagMatrix_ckoR2 <- getTagMatrix(ckoR2, windows = promoter)
```

Annotate:

```{r}
annot.ckoR2 = annotatePeak(ckoR2, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Mm.eg.db")
```

<details>
    <summary>Annotation results</summary>


Annotated peaks generated by ChIPseeker

83085/83089  peaks were annotated

Genomic Annotation Summary:
 
| Feature               | Frequency     |
|-----------------------|---------------|
| Promoter (<=1kb)      | 23.41577902   |
| Promoter (1-2kb)      | 3.77083709    |
| Promoter (2-3kb)      | 3.06433171    |
| 5' UTR               | 0.22025636    |
| 3' UTR               | 1.82343383    |
| 1st Exon             | 1.66335680    |
| Other Exon           | 3.00174520    |
| 1st Intron           | 11.68321598   |
| Other Intron         | 24.45327075   |
| Downstream (<=300)   | 0.08064031    |
| Distal Intergenic    | 26.82313294   |			


</details>


#### Sample CKO_Rep3

Importing peak sets

```{r}
ckoR3 <- readPeakFile("data/macs2/narrowPeak/cKO_H3K27ac_ChIPseq_REP3_peaks.narrowPeak"

```
Generate TagMatrix:

```{r}
#promoter <- getPromoters(TxDb = txdb, upstream = 2000, downstream = 2000)
tagMatrix_ckoR3 <- getTagMatrix(ckoR3, windows = promoter)
```
Annotate:

```{r}
annot.ckoR3 = annotatePeak(ckoR3, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Mm.eg.db")
```

<details>
    <summary>Annotation results</summary>

 
Annotated peaks generated by ChIPseeker

91046/91053  peaks were annotated

Genomic Annotation Summary:
 
| Feature               | Frequency     |
|-----------------------|---------------|
| Promoter (<=1kb)      | 22.96092085   |
| Promoter (1-2kb)      | 3.91560310    |
| Promoter (2-3kb)      | 3.11271226    |
| 5' UTR               | 0.19660391    |
| 3' UTR               | 1.72879643    |
| 1st Exon             | 1.57612635    |
| Other Exon           | 2.81725721    |
| 1st Intron           | 11.52494344   |
| Other Intron         | 24.53594886   |
| Downstream (<=300)   | 0.08237594    |
| Distal Intergenic    | 27.54871164   |			


</details>


#### Sample WT_Rep2

Importing peak sets
```{r}
wtR2 <- readPeakFile("data/macs2/narrowPeak/WT_H3K27ac_ChIPseq_REP2_peaks.narrowPeak"

```

Generate TagMatrix:

```{r}
#promoter <- getPromoters(TxDb = txdb, upstream = 2000, downstream = 2000)
tagMatrix_wtR2 <- getTagMatrix(wtR2, windows = promoter)
tagHeatmap(tagMatrix_wtR2)
```

Annotate:

```{r}
annot.wtR2 = annotatePeak(wtR2, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Mm.eg.db")
```

<details>
    <summary>Annotation results</summary>

 
Annotated peaks generated by ChIPseeker

111961/111978  peaks were annotated

Genomic Annotation Summary:
 
| Feature               | Frequency     |
|-----------------------|---------------|
| Promoter (<=1kb)      | 19.52644224   |
| Promoter (1-2kb)      | 3.83079822    |
| Promoter (2-3kb)      | 3.11894320    |
| 5' UTR               | 0.22418521    |
| 3' UTR               | 1.85957610    |
| 1st Exon             | 1.55500576    |
| Other Exon           | 3.08768232    |
| 1st Intron           | 12.18727950   |
| Other Intron         | 26.07425800   |
| Downstream (<=300)   | 0.09378266    |
| Distal Intergenic    | 28.44204678   |


</details>



#### Sample WT_Rep3

Importing peak sets
```{r}
wtR3 <- readPeakFile("data/macs2/narrowPeak/WT_H3K27ac_ChIPseq_REP3_peaks.narrowPeak"
```

Generate TagMatrix:
```{r}
#promoter <- getPromoters(TxDb = txdb, upstream = 2000, downstream = 2000)
tagMatrix_wtR3 <- getTagMatrix(wtR3, windows = promoter)
```

Annotate:

```{r}
annot.wtR3 = annotatePeak(wtR3, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Mm.eg.db")
```

<details>
    <summary>Annotation results</summary>

 
Annotated peaks generated by ChIPseeker

81024/81030  peaks were annotated

Genomic Annotation Summary:
 
| Feature               | Frequency     |
|-----------------------|---------------|
| Promoter (<=1kb)      | 25.30731635   |
| Promoter (1-2kb)      | 3.72111967    |
| Promoter (2-3kb)      | 2.85470972    |
| 5' UTR               | 0.18883294    |
| 3' UTR               | 1.64025474    |
| 1st Exon             | 1.58965245    |
| Other Exon           | 2.62144550    |
| 1st Intron           | 11.02142575   |
| Other Intron         | 23.74234795   |
| Downstream (<=300)   | 0.09626777    |
| Distal Intergenic    | 27.21662717   |


</details>



### Average profiles of all the samples

We can plot average profile of all the samples together. For that we need to make a list of list of tagMatrix.

```{r, message=FALSE, warning=FALSE, results='hide'}
tagMatrix_list <- list(
  cKO_Rep1 = tagMatrix_ckoR1,
  cKO_Rep2 = tagMatrix_ckoR2,
  cKO_Rep3 = tagMatrix_ckoR3,
  WT_Rep1 = tagMatrix_wtR1,
  WT_Rep2 = tagMatrix_wtR2,
  WT_Rep3 = tagMatrix_wtR3
)

plotAvgProf(tagMatrix_list, xlim = c(-2000, 2000), facet = "row")
```

<p align="center">
<img src="../img/Chipseq_av.profiles_all.png"  width="600">
</p>





### Chip peak annotation distrubution of loci with respect to TSS all samples
Lets plot the peak annotation distribution of loci with respect to TSS for all the samples together.

We can supply the same annotation list that we made above.

```{r}
plotDistToTSS(annot_list)
```
<p align="center">
<img src="../img/Chipseq_tss_all.png"  width="600">
</p>

### Functional profiles comparison

The `clusterProfiler` package provides `compareCluster()` function for comparing functional annotations gene clusters and can be adopted to compare different ChIP peak experiments.

```{r}
genes <- lapply(annot_list, function(i) as.data.frame(i)$geneId)
names(genes) <- sub("_", "\n", names(genes))

compPathway <- compareCluster(geneCluster = gene, 
                         fun = "enrichPathway",
                         organism = "mouse",
                         pvalueCutoff = 0.05, 
                         pAdjustMethod = "BH")
```

**Dotplot**

```{r}
dotplot(compPathway, showCategory = 10, title = "Reactome Pathway Enrichment Analysis")
```
<p align="center">
<img src="../img/reac_all.png"  width="600">
</p>

### Exercise:
1. Carryout functional comparison using KEGG database.
   

***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
