---
title: "From sequence reads to peak calls"
author: "Meeta Mistry, Upendra Bhattarai, Heather Wick"
date: "Aug 13th, 2024"
---

Contributors: Meeta Mistry, Upendra Bhattarai, Will Gammerdinger

Approximate time: 20 minutes

## Learning Objectives

* Describe the workflow from sequence reads to peak calls


## From sequence reads to peak calls

In this lesson, we highlight for you the important steps involved in a typical workflow for the analysis of ChIP-seq or related data. If you are looking for more in-depth information on the background and theory for each step, we suggest looking at the [Understanding chromatin biology using high throughput sequencing](https://hbctraining.github.io/Intro-to-ChIPseq-flipped/schedule/) workshop materials. 

> **NOTE:** When starting out with your experiment, there are numerous quality control considerations at the bench when preparing your samples. We have highlighted some of the important points [within this lesson](). Investing time early in the experiement to ensure good quality samples will pay off with meaningful and reproducible results down the line.


<p align="center">
<img src="../img/chipseq_peakworkflow_sept2021.png" width="600">
</p>


## Sequence data QC
The FASTQ files you obtain from the sequencing facility will first need to be assessed for quality. Here, we look at metrics like ... 

At this stage we are flagging samples with values deviating from expected ranges. If quality metrics are very poor, this can results in removal of the sample from downstream analysis.
**Specific details/images on what we asssess can be taken from here: https://hbctraining.github.io/Intro-to-ChIPseq-flipped/lessons/troubleshooting_chipseq_partI.html**

## Alignment to genome 
Next, we take our high quality reads an map them to genome. We do this using Bowtie2...


## Filtering BAM files
* mention duplicates, multi-mappers, blacklisted regions
* 

## Peak calling
Peak calling, the next step in our workflow, is a computational method used to identify areas in the genome that have been enriched with aligned reads as a consequence of performing a ChIP-sequencing experiment.

For ChIP-seq experiments, what we observe from the alignment files is a **strand asymmetry with read densities on the +/- strand, centered around the binding site**. The 5' ends of the selected fragments will form groups on the positive- and negative-strand. The distributions of these groups are then assessed using statistical measures and compared against background (input or IgG samples) to determine if the site of enrichment is likely to be a real binding site.

<p align="center">
<img src="../img/plos_chipseq_arrow.png" width = "600">
</p>

*Image source: [Wilbanks and Faccioti, PLoS One 2010](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0011471)*





***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
