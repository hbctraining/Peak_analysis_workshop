---
title: "Peak visualization using a genome browser"
author: "Heather Wick, Upendra Bhattarai, Meeta Mistry"
date: "Aug 16th, 2024"
---

Contributors: Heather Wick, Upendra Bhattarai, Meeta Mistry, Will Gammerdinger

Approximate time: 

## Learning Objectives

* Visualize peaks in IGV (Integrative Genomics Viewer)
* Qualitatively assess whether DiffBind correctly assess differential binding
* Use R package (graphViz/rtracklayer) to plot and study differential peak enrichment between experimental groups


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

### Selecting the mm10 reference genome

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

### Download files for visualization

Next, we need to download some files that we will be visualizing in IGV. Right-click [this link](https://www.dropbox.com/s/7pzb1yvpgopzar9/bigWig.zip?st=jgzjx01b&dl=0) and select "Save Link As..." to download the ZIP compressed directory of the BigWig files from Dropbox. Place this file within the `data` folder of your `Peak_analysis`. Repeat this process with [this link](https://www.dropbox.com/scl/fi/lwwm49gfgh433pzhxzdnu/DiffBind.zip?rlkey=b8xwzcftsn0dyjha1uzpdnn3q&st=jfyqpmco&dl=0) for the BED file of differentially called peaks from DiffBind and also place this file within the `data` folder of your `Peak_analysis`.

<details><summary><b>Click here to see how to create bigWig files</b></summary>

We will need to use the output from the <code>picard</code>`'s <code>CollectAlignmentSummaryMetrics</code> tool once again. As a reminder the code for running this would be:<br>

<pre>
# Run picard CollectAlignmentSummaryMetrics for a sample
java -jar picard.jar CollectAlignmentSummaryMetrics \
  --INPUT $INPUT_BAM \
  --REFERENCE_SEQUENCE $REFERENCE \
  --OUTPUT $OUTPUT_METRICS_FILE
</pre><br>

We will be interested in the value associated with <code>PF_READS_ALIGNED</code>. This is the number of your mapped reads. We will use this number to create a scaling factor in the next step to create a bedGraph file.<br>

The full documentation for <code>picard CollectAlignmentSummaryMetrics</code> can be found <a href="https://gatk.broadinstitute.org/hc/en-us/articles/360040507751-CollectAlignmentSummaryMetrics-Picard">here</a>.<br>

In order to create a bedGraph file we will use <code>bedtools</code>'s <code>genomecov</code> tool.<br>

<pre>
$SCALE_FACTOR=`awk 'BEGIN { print 1000000 / $MAPPED_READ_COUNT }'`
  
bedtools genomecov \
  -ibam $INPUT_SORTED_BAM \
  -bg \
  -scale $SCALE_FACTOR | \
  sort -k1,1 -k2,2n \
  > $OUTPUT_FILE
</pre><br>

We will need to create a bash variabel to hold a scale factor to a miilion which we will call in the <code>bedtools</code> command. The command to do this is:<br>
<code>$SCALE_FACTOR=`awk 'BEGIN { print 1000000 / $MAPPED_READ_COUNT }'`</code><br>

<ul><li><code>bedtools genomecov</code> - Calls the <code>bedtools</code>'s <code>genomecov</code> tool
  <li><code>-ibam</code> - Input sorted BAM file</li>
  <li><code>-bg</code> - Output depth in bedGraph format</li>
  <li><code>-scale</code> - Scale factor used for scaling the data</li>
  <li><code>sort -k1,1 -k2,2n</code> - The output file is unsorted and needs to be sorted by chromosome and chromosome start location
  <li><code>&gt; $OUTPUT_FILE</code> - Write to an ouptut file</li>
</ul><br>

The full documentation for <code>bedtools genomecov</code> can be found <a href="https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html">here</a>.

At this point it is possible to load bedGraph formatted files into IGV, but they are much larger than BigWig files, so we will convert our bedGraph files into BigWig files using a tool called <code>bedGraphToBigWig</code>:<br>

<pre>
bedGraphToBigWig \
    $BEDGRAPH_INPUT \
    $CHROMOSOME_SIZES_FILE \
    $BIGWIG_OUTPUT 
</pre><br>

We can break this command down as:<br>

<ul><li><code>bedGraphToBigWig</code> - Calls the <code>bedGraphToBigWig</code> tool</li>
  <li><code>$BEDGRAPH_INPUT</code> - The bedGraph input file</li>
  <li><code>$CHROMOSOME_SIZES_FILE</code> - A tab-delimited file with chromosomes in the first column and their associated sizes in the second column. For more information on how to create this file, click on the dropdown below called "Click here to see how to create a chromosome sizes file"</li>
  <li><code>$BIGWIG_OUTPUT</code> - The BigWig output file</li>
</ul><br>

The download for <code>bedGraphToBigWig</code> can be found <a href="https://github.com/ENCODE-DCC/kentUtils">here</a>.<br>

<details><summary><b>Click here to see how to create a chromosome sizes file</b></summary>
There are several ways to make a tab-delimited file with the chromosomes in the first column and their asscoiated sizes in the second column. One way is to use the <code>samtools</code> package <code>faidx</code> to create a FASTA index file. The documentation to run this cool can be found <a href="https://www.htslib.org/doc/samtools-faidx.html">here</a>. The command to run the <code>samtools faidx</code> tool is:

<pre>
  samtools faidx \
    $REFERENCE_GENOME_FASTA
</pre><br>

We can break this command into two parts:<br>

<ul><li><code>samtools faidx</code> - Call the <code>samtools faidx</code> tool</li>
  <li><code>$REFERENCE_GENOME_FASTA</code> - The reference genome FASTA file</li>
</ul><br>

The full documentation for using <code>samtools faidx</code> can be found <a href="https://www.htslib.org/doc/samtools-faidx.html">here</a>.<br>

However, this output will have a few more column than you need. You only need the first two columns, so we can use <code>awk</code> to parse out the first two columns:

<pre>
  awk \
  -v OFS='\t' \
  '{print $1, $2}' \
  $FASTA_INDEX_FILE \
  > $CHROMOSOME_SIZES_FILE
</pre>

This command is composed of a few parts:

<ul><li><code>awk</code> - Calling <code>awk</code></li>
  <li><code>-v OFS='\t'</code> - Output as tab-delimited</li>
  <li><code>'{print $1, $2}'</code> - Print the first and second columns</li>
  <li><code>$FASTA_INDEX_FILE</code> - Input FASTA index file that was made with the above <code>samtools faidx</code> command</li>
  <li><code>&gt; $CHROMOSOME_SIZES_FILE</code> - Output Chromosome sizes file that we can use in out <code>bedGraphToBigWig</code> command</li>
</ul>
<hr />
</details>
We discuss an alternative way of generating BigWig files in our <a href="https://hbctraining.github.io/Intro-to-ChIPseq-flipped/lessons/08_creating_bigwig_files.html">Chromatin Biology materials</a>.
<hr />
</details>

### Load a track

Now that we have the data that we would like to visualize, let's go ahead and load it into IGV. Many file formats you load into IGV will load as a single "track", or horizonatla row of genomic data. However some, like BAM/SAM, will load as multiple tracks. Let's go ahead and load the BigWig track for our cKO IP replicate 3 sample.

In order to load a track, left-click on "File" in the top-left and select "Load from File...":

<p align="center">
<img src="../img/IGV_File_load_track_with_caption.png"  width="600">
</p>

Next, navigate in your file browser to the file you'd like to load, then left-click to select it and left-click "Open". In this case we are trying to open the file called "cKO_H3K27ac_ChIPseq_REP3.bigWig", which should be inside your "bigWig" directory within your "data" directory:

<p align="center">
<img src="../img/IGV_select_track_with_caption.png"  width="600">
</p>

After loading the cKO IP replicate 3 BigWig track, your IGV session should look like:

<p align="center">
<img src="../img/IGV_loaded_track_with_caption.png"  width="600">
</p>

### Exercise 

1. Load more Bigwig tracks - WT IP and Input

2. Add DiffBind BED

### Moving tracks


### Exercise

move to the order we show

### Zooming in and out

### Jumping to regions

### Modifying Tracks

#### Adjust track height

height



#### Adjust data range

#### Adjust color

#### Rename Tracks

#### Load from URL
Remove

### Save a Session

### Load a Session



***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
