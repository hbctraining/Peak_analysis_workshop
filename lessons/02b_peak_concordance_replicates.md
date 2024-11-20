---
title: "Concordance of peaks across replicates"
author: "Heather Wick, Upendra Bhattarai, Meeta Mistry"
date: "Aug 13th, 2024"
---

Contributors: Heather Wick, Upendra Bhattarai, Meeta Mistry

Approximate time: 

## Learning Objectives

* Evaluate signal enrichment vs rank plots for samples in dataset
* Plot PCA and hierarchical clustering to assess inter-sample variability
  
## Concordant peak calls between samples

In the previous lesson, we evaluated quality metrics concerning peaks and reads in peaks in individual samples. But these aren't the only ways to measure quality in our data set. We did look for consistency of these metrics across all of our samples, but now it's **time for a closer look at our samples to see how they comapre for samples within treatment groups and between groups**. To do this we will be visualizing our data using two different types of data: 1) read count distribution across the genome and 2) siganl enrichment within regions called as peaks.

## Read count density 
Two commonly used methods for evaluating sample similarity is the use of Principal Components Analysis (PCA) and hierarchical clustering. In order to implement both of these, we need the **appropriate inputs**. The requirement is a data matrix. In the case of ChIP-seq or related data, this **matrix would have samples in the columns and genomic regions in the rows**. This matrix can be created using a command-line tool called [`multiBamSummary`](https://deeptools.readthedocs.io/en/develop/content/tools/multiBamSummary.html) from the [deepTools suite](https://deeptools.readthedocs.io/en/develop/index.html) for exploring deep sequencing data. Since this is an R-based workshop, **we have created the matrix for you**, but have provided the code in the drop-down below if you wanted to use it on your own data.


<details>
<summary><b>Click here for the code to create a read count matrix using deepTools</b></summary>
<br>The command <code>multiBamSummary</code> computes the read coverages for genomic regions for two or more BAM files. The analysis can be performed for the entire genome by running the program in ‘bins’ mode. If you want to count the read coverage for specific regions only (for example only consensus regions), you can use the BED-file mode instead. <br><br>
The standard output of multiBamSummary is a compressed numpy array (.npz), which can be used with other functions in deepTools to create PCA plots and correlation heatmaps. <b>Instead we opt to get a tab-delimited file to allow us flexibility to create our own plots by using the <code>--outRawCounts</code> parameter.</b><br>

The command to run this is:</br></br>
<pre>
multiBamSummary bins \
	  --bamfiles cKO_H3K27ac_ChIPseq_REP1.mLb.clN.sorted.bam cKO_H3K27ac_ChIPseq_REP2.mLb.clN.sorted.bam cKO_H3K27ac_ChIPseq_REP3.mLb.clN.sorted.bam \
               cKO_H3K27ac_input_REP1.mLb.clN.sorted.bam cKO_H3K27ac_input_REP2.mLb.clN.sorted.bam cKO_H3K27ac_input_REP3.mLb.clN.sorted.bam \
               WT_H3K27ac_ChIPseq_REP1.mLb.clN.sorted.bam WT_H3K27ac_ChIPseq_REP2.mLb.clN.sorted.bam WT_H3K27ac_ChIPseq_REP3.mLb.clN.sorted.bam \
               WT_H3K27ac_input_REP1.mLb.clN.sorted.bam WT_H3K27ac_input_REP2.mLb.clN.sorted.bam WT_H3K27ac_input_REP3.mLb.clN.sorted.bam \
	  --outFileName multiBamsummary.npz \
	  --labels cKO_1 cKO_2 cKO_3 cKO_i1 cKO_i2 cKO_i3 WT_1 WT_2 WT_3 WT_i1 WT_i2 WT_i3 \
	  -p 6 \
	  --outRawCounts multiBAMsummary.tab
</pre></br>

</details>

We have created two versions of the read count matrix: one that contains ChIP and input samples and another that contains ChIP samples only. Both of these files can be found in the `data/multBamSummary` folder. Let's read in the data:

```
# Read in data

```

> #### Should we normalize the data matrix?


## Principal Component Analysis (PCA)

Principal Component Analysis (PCA) is a technique used to emphasize variation and bring out strong patterns in a dataset (dimensionality reduction). This is a very important technique used in the QC and analysis of ChIPseq.

If you've done any RNAseq or single cell analysis, you're likely familiar with the concept of PCA, as this technique is used to evaluate. If you'd like to refamiliarize yourself on the details of how PCA is calculated, we recommend you read our materials [[here]](https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/principal_component_analysis.html) or [[watch this useful video from StatQuest]](https://www.youtube.com/watch?v=_UVHneBUBW0&ab_channel=StatQuestwithJoshStarmer).

### Interpreting PCA plots
Essentially, if two samples have similar levels of expression peak enrichment that contribute significantly to the variation represented by a given PC (Principal Component), they will be plotted close together on the axis that represents that PC. Therefore, we would expect that biological replicates to have similar scores (because our expectation is that the same genes are changing) and cluster together. This is easiest to understand by visualizing some example PCA plots.

We can run PCA to evaluate the variation amongst our samples and whether or not the greatest sources of variation in the data (PC1 and PC2) can be attributed to the factors of interest in this experiment.

We will be plotting the Variance Stabalized Transformation (VST) counts from our dds object and usin a package called degPCA which automatically performs the PCA dimensionality reduction and plots the results. To start, we will look at both PC1 vs PC2 and PC3 vs PC4 colored by the condition "genotype", since that is the primaro factor of interest for us:

```
library(degPCA)
pca1v2 <- degPCA(assays(dds)$vst, coldata_for_dds,
              condition = "genotype", data = T)[["plot"]]
pca3v4 <- degPCA(assays(dds)$vst, coldata_for_dds,
              condition = "genotype", data = T, pc1="PC3", pc2="PC4")[["plot"]]

pca1v2
pca3v4

```

####INSERT pc1v2 AND pca3v4 SIDE BY SIDE

We are fortunate that our data separates on PC1, which accounts for the most amount of variance in the data, by genotype. If this were not the case, we would want to consider coloring our data points by other aspects of the metadata (`coldata_for_dds` object). In our case, that would be mostly technical factors, but if our samples were processed in different batches, or performed on different dates; or if the samples were from tissues with different sexes or other features, these would be important features to label our plot by when looking at PCA. If any factor contributed to a large amount of variance in the data, we might need to take further steps to take these factors into account

## Hierarchical clustering

Inter-correlation analysis (ICA) is another way to look at how well samples cluster by plotting the correlation between the peak regions of the samples. We can also add metadata to this plot to see if any factors explain any of the clustering on any level.

First, we will compute the correlation of our VST-normalized peak counts using the Pearson correlation ####INCLUDE EXPLANATION? RATIONAL? INFO ABOUT COR() AND OTHER OPTIONS?:

```
vst_cor <- cor(assays(dds)$vst)
```

Then, we will subset our metadata. Here, we are just selecting genotype, but it is a good idea to plot other things, like batch, or other technical factors, especially if we saw our data segregate by them on any of the PCs

```
colma=coldata_for_dds %>% as.data.frame()
rownames(colma) <- colma$sample
colma <- colma[rownames(vst_cor), ]
colma <- colma %>% dplyr::select(genotype)
```

Finally, we can make the heatmap itself:

```
p <- pheatmap(vst_cor, 
         annotation = colma,
         #annotation_colors = anno_colors,
         show_rownames = T, 
         show_colnames = T, 
         )
p
```

####INSERT HEATMAP HERE

As expected given what we saw in the PCA plots, our samples cluster nicely by genotype. If we had more samples and plotted more metadata, we might also be able to see whether batch or other biological or technical factor affected the clustering, and on what level.



####INSERT IMAGE HERE


## Signal concordance across peaks



## Signal enrichment vs peak rank

One way to evaluate concordance of peaks between samples is a signal enrichment vs peak rank plot. This shows the rank of each peak vs the strength of each peak for each repilcate. It will help us evaluate the number of peaks we would retain if thresholding by peak enrichment:

```
ggplot(peaks, aes(x = peak_rank, y = peak_enrichment, color = sample)) + 
  geom_line() +
  xlab("Peak rank") + ylab("Peak enrichment")
```

It is also valuable to see how this differs between replicates within a sample group. By adding faceting to our plot, we can more easily look at samples within each genotaype to evaluate sample similarity:

```
ggplot(peaks, aes(x = peak_rank, y = peak_enrichment, color = sample)) + 
  geom_line() +
  facet_grid(. ~ genotype) +
  xlab("Peak rank") + ylab("Peak enrichment")
```

####INSERT PR v PE PLOT HERE

Our data looks pretty consistent between samples until enrichment is low ####COMMENT ON MEANING BHIND THIS

## Histogram of quality scores

Here, we plot a histogram of peak signal values for each sample. This plot can be used to help determine a minimum value for peak enrichment that can be used for filtering. 

```
ggplot(peaks, aes(x = peak_enrichment, fill = .data[[params$factor_of_interest]])) + 
  geom_histogram(aes(peak_enrichment)) +
  scale_fill_cb_friendly() +
  xlab("Peak enrichment")
```

####INSERT IMAGE HERE

####WRITE NOTE ABOUT HOW OUR DATA LOOKS AND WHETHER WE SHOULD DO ANY FILTERING



***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
