---
title: "Peak annotation and visualization"
author: "Heather Wick, Upendra Bhattarai, Meeta Mistry"
date: "Aug 16th, 2024"
---

Contributors: Heather Wick, Upendra Bhattarai, Meeta Mistry

Approximate time: 

## Learning Objectives

* Annotate peaks with genomic features using Chipseeker
* Visualize the annotation and compare peak coverage between experimental groups.
* Explore the biological content using functional enrichment analysis



## Overview

In this tutorial we use ChIPseeker, a Bioconductor package in R to annotate peaks with genomic feaures, explore various ways to visualize them and compare the peak profile and coverage between the experimental groups. We will also perform functional annotation and analysis to explore the biology in our dataset.


Notes:

- Annotation for individual samples and for consensus peaks for each sample group.
- Peak annotation of genomic features (not on target genes) [see](https://hbctraining.github.io/Intro-to-ChIPseq/lessons/12_functional_analysis.html)
    - Peak distance from TSS and other genomic feaures.
    - ChiPSeeker plots [see](https://github.com/hbctraining/Peak_analysis_workshop/blob/main/lessons/OLD_ChipSeeker_analysis.md#visualization-with-chipseeker)
    - Talk about Deeptools as an alternative using bam files. point to [chromatin biology workshop material](https://github.com/hbctraining/Intro-to-ChIPseq-flipped/blob/main/lessons/09_data_visualization.md)
***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
