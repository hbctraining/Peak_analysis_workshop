---
title: "Peak visualization using a genome browser"
author: "Heather Wick, Upendra Bhattarai, Meeta Mistry"
date: "Aug 16th, 2024"
---

Contributors: Heather Wick, Upendra Bhattarai, Meeta Mistry, Will Gammerdinger

Approximate time: 

## Learning Objectives

* Visualize peaks in IGV (Integrative Genomics Viewer).
* Describe the file formats used for visualiztion.
* Use R package (graphViz/rtracklayer) to plot and study differential peak enrichment between experimental groups.


## Overview

<p align="center">
<img src="../img/IGV_logo.jpeg"  width="150">
</p>

When working with next-generation sequecing data, it can be helpful to visualize your analysis. [Integrate Genomics Viewer (IGV)](https://igv.org/) is an open-source desktop genome visualization tool, which supports visualization for a wide variety of file formats including:

- BED
- bigWig
- VCF
- BAM/SAM
- wig
- bedGraph
- narrowPeak/broadPeak
- GFF3
- GTF
- More and with a complete list of suppported file formats can be found [here](https://igvteam.github.io/igv-webapp/fileFormats.html)

Below we will learn about the file formats used to visualize peaks in IGV and use it to explore peaks in our experiment samples. We will also use graphViz/rtracklayer to plot the peak signals in R. 

## IGV

### Setting Up

Open IGV on your computer. The first thing we will need to use the appropriate reference genome (mm10). There should be a dropdown menu in the top-left of IGV, which let's you select your preferred reference genome. If your selected reference genome is already **mm10**, then you can skip the next few steps. Otherwise, left-click on the dropdown menu:

<p align="center">
<img src="../img/IGV_default_with_caption.png"  width="600">
</p>

You can see a few options for refernce genomes to select from. If you see **mm10**, then you can select it. Otherwise, left-click on "Click for more ..."

<p align="center">
<img src="../img/IGV_reference_dropdown_with_caption.png"  width="600">
</p>

Find "Mouse mm10" from the menu and left-click "OK"

<p align="center">
<img src="../img/IGV_reference_menu_with_caption.png"  width="600">
</p>

***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
