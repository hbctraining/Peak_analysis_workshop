---
title: "Peak overlaps"
author: "Heather Wick, Upendra Bhattarai, Meeta Mistry"
date: "Aug 13th, 2024"
---

Contributors: Heather Wick, Upendra Bhattarai, Meeta Mistry

Approximate time: 

## Learning Objectives

* Describe IRanges and GRanges in R along with some basic functions
* Identify overlaps between replicates and visualize with VennDiagrams and upSetR
* Extract consensus regions across replicates for a sample group
  
## Peak Overlaps

So far we have looked at quality control metrics for individual samples, as well as statistical concordance between samples. An additional way to look at sample similarity is to look at peak overlap between samples; that is to ask, what peaks are in common between replicate samples within our treatment groups? Looking at peak overlaps serves two purposes: 

1) It is another way of measuring consistency between our samples
2) We can create a set of **consensus peaks** (peaks in common between samples within a treatment group), in which we are more confident, as these are less likely to be miscalls due to background noise or other technical variation. These consensus peaks can be used in downstream visualization and analysis.

### Setup 
Let's begin by creating a new script for this lesson, and let's call it `peak_overlaps_analyis.R`. We can add a header to our script and start with a section to load the required libraries:

```
## Finding overlapping peaks

# Load libraries
library(ChIPpeakAnno)
library(UpSetR)
```

## Essential tools: IRanges and GenomicRanges

You may be familiar with **bedtools** as a useful command line too for manipulating bed files, including finding overlap of genomic regions. Whenever we are doing anything involving overlap of genomic ranges in R, two additional essential tools are of great help: **IRanges** and **GenomicRanges**. These packages allow us to convert bed files and other, more complex and/or binary genomic coordinate files, such as bam files, narrowPeak files and bigWigs into objects in R, and come with a number of different functions that allow us to find overlaps, exclusions, or nearest genomic features, among other things.

IRanges and GRanges are data structures that can be used to solve a variety of problems, typically related to annotating and visualizing the genome. These data structures are very fast and efficient. Both packages contain very extensive and useful vignettes which we have linked below. If you find yourself stuck with these data structures at any point, we encourage you to review them!

### IRanges: the minimal representation of a range in a single space

An IRanges object in R is a very simple representation of a coordinate in a single space (chromosome, in our case), with a `start`, and `end` and a `width`. To construct an IRanges object, we call the IRanges constructor. Ranges are normally specified by passing two out of the three parameters: start, end and width:

```
# Example IRanges
ir <- IRanges(start=1:5, width=5:1)
ir

IRanges object with 5 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         1         5         5
  [2]         2         5         4
  [3]         3         5         3
  [4]         4         5         2
  [5]         5         5         1

```
> **NOTE:** The [IRanges vignette](https://bioconductor.org/packages/release/bioc/vignettes/IRanges/inst/doc/IRangesOverview.pdf), is a great place to learn more about IRanges objects and how to manipulate them.

### GenomicRanges, or GRanges: ranges in multiple spaces

A GRanges object is a little more complex. It lets us store IRanges in multiple spaces (ie multiple chromosomes). In addition to chromosome, a GRanges object also indicates the strand for each region. These objects can also hold additional meta data. GRanges provide a way to store and manipulate sets of genomic regions. In our example, we will be using it to store the peak calls from each of our samples. Let's start with simple example to create a GRanges object using the `GRanges()` constructor:

```
gr <- GRanges(ranges=IRanges(start=c(100, 200), end=c(199, 299)), 
              seqnames=c("chr2L", "chr3R"),
              strand=c("+", "-"))
gr

GRanges object with 2 ranges and 0 metadata columns:
      seqnames    ranges strand
         <Rle> <IRanges>  <Rle>
  [1]    chr2L   100-199      +
  [2]    chr3R   200-299      -
  -------
  seqinfo: 2 sequences from an unspecified genome; no seqlengths
```

You can see that one of the required inputs is an IRanges object, and so functionality in this package is very dependent on the basics of IRanges. For more detailed information on GenomicRanges, we enocurage you to browse through the [GRanges vignette](https://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html).


Once you have your genomic coordinate data stored in one of these data structures, there are many functions that allow you to easily manipulate the data. There are **functions for basic interval operation** like `shift()`, `reduce()`, `flank()`, `intersect()` and so much more. We have linked for you a [helpful cheatsheet](https://rpubs.com/Pazz/bioc_cheat_sheet) which describes commonly used functions and in some use cases. In this lesson we will first convert our peak files into GRanges and then we will use a package called [ChipPeakAnno](https://bioconductor.org/packages/release/bioc/html/ChIPpeakAnno.html) which will provide wrapper functions that allow us to easily operate on our peak data and pull out the infformation we need.

## Reading in narrowPeak files as Granges objects
To create GRanges objects for our peak files, we are going to use the function `toRanges()` from the ChIPpeakAnno package. Note that there are other alternative packages that have functionality to do this. This function requires at minimum:

* Genomic coordinate data (as a path to file, or a data frame)
* Format (narrowPeak, broadPeak, BED)

We will provide file paths that we had stored in the `sample_files` variable in the [previous lesson](02b_sample_similarity.md#signal-concordance-across-peaks), and we will **use a for loop to create a GRanges object for each sample**. Note that we are assigning the object to overwrite the variables which had the peak data stored as data frames.

```
# If you cannot find sample_files in your environment, uncomment the code below
# Get all narrowpeak file names and path
# sample_files <- list.files(path = "./data/macs2/narrowPeak/", full.names = T)


# Reassign vars so that they are now GRanges instead of dataframes
for(r in 1:length(sample_files)){
  obj <- ChIPpeakAnno::toGRanges(sample_files[r], format="narrowPeak", header=FALSE)  
  assign(vars[r], obj)
}
```

Now let's take a quick look at one of the samples:

```
WT_H3K27ac_ChIPseq_REP1
GRanges object with 100570 ranges and 5 metadata columns:
                                                  seqnames          ranges
                                                     <Rle>       <IRanges>
       WT_H3K27ac_ChIPseq_REP1_peak_1                 chr1 3094273-3094445
       WT_H3K27ac_ChIPseq_REP1_peak_2                 chr1 3095208-3095464
       WT_H3K27ac_ChIPseq_REP1_peak_3                 chr1 3113322-3113761
       WT_H3K27ac_ChIPseq_REP1_peak_4                 chr1 3119377-3120158
       WT_H3K27ac_ChIPseq_REP1_peak_5                 chr1 3120430-3120942

                                    strand |     score signalValue    pValue
                                       <Rle> | <integer>   <numeric> <numeric>
       WT_H3K27ac_ChIPseq_REP1_peak_1      * |        28     3.71719   4.55519
       WT_H3K27ac_ChIPseq_REP1_peak_2      * |        89     6.27302  10.90470
       WT_H3K27ac_ChIPseq_REP1_peak_3      * |        97     7.01319  11.73230
       WT_H3K27ac_ChIPseq_REP1_peak_4      * |       277    13.74740  30.10350
       WT_H3K27ac_ChIPseq_REP1_peak_5      * |       193    10.71350  21.52090

                                     qValue      peak
                                      <numeric> <integer>
       WT_H3K27ac_ChIPseq_REP1_peak_1   2.87376       100
       WT_H3K27ac_ChIPseq_REP1_peak_2   8.96056       111
       WT_H3K27ac_ChIPseq_REP1_peak_3   9.76906        75
       WT_H3K27ac_ChIPseq_REP1_peak_4  27.77500       291
       WT_H3K27ac_ChIPseq_REP1_peak_5  19.34280       160
```

We see that each peak is stored with it's ranges and all associated infomration generated by MACS2.

## Find overlapping peaks

Once we have read in our peak files, we can look for overlapping genomic ranges using the `findOverlappingPeaks()` function:

```
Will need to take code from section below and make it work for this but not using the peaks object
```

### Venn Diagram

### UpsetR plot

####THIS CODE USES PEAKS OBJECT BUT WE WILL NEED TO GET IT TO WORK ON OUR GRANGES OBJECTS

```
for (current_sample_group in unique(peaks$sample_group)){
  cat("## ", current_sample_group, "\n")

  peaks_sample_group <- peaks %>% filter(sample_group == current_sample_group)
  
  peaks_sample_group_granges <- sapply(
    unique(peaks_sample_group$sample), 
    function(current_sample) {
      ChIPpeakAnno::toGRanges(
        peaks_sample_group %>% filter(sample == current_sample),
        format = ifelse(grepl('broadPeak', peaks_dir), 'broadPeak', 'narroPeak')
      )
    }
  )
  
  # maxgap defaults to -1 which means that two peaks will be merged if they overlap by    at least 1 bp
  # connectedpeaks examples (https://support.bioconductor.org/p/133486/#133603), if 5     peaks in group1 overlap with 2 peaks in group 2, setting connectedPeaks to "merge"      will add 1 to the overlapping counts
  overlaps <- findOverlapsOfPeaks(peaks_sample_group_granges, connectedPeaks = 'merge')
  
  set_counts <- overlaps$venn_cnt[, colnames(overlaps$venn_cnt)] %>% 
    as.data.frame() %>% 
    mutate(group_number = row_number()) %>%
    pivot_longer(!Counts & !group_number, names_to = 'sample', values_to = 'member') %>%
    filter(member > 0) %>%
    group_by(Counts, group_number) %>% 
    summarize(group = paste(sample, collapse = '&'))
  
  set_counts_upset <- set_counts$Counts
  names(set_counts_upset) <- set_counts$group

  p <- upset(fromExpression(set_counts_upset), order.by = "freq", text.scale = 1.5)
  print(p)
  
  cat('\n\n')

}
```

## Finding Consensus peaks




***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
