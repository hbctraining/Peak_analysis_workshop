---
title: "Differential binding analysis"
author: "Heather Wick, Upendra Bhattarai, Meeta Mistry"
date: "Aug 16th, 2024"
---

Contributors: Heather Wick, Upendra Bhattarai, Meeta Mistry

Approximate time: 

## Learning Objectives

* Carryout data exploratory analysis and visualization for sample concordance.
* Learn to use DiffBind for differential binding analysis
* Evaluate the results and export the output.




## Overview

### SHORT INTRO ON DIFFERENT WAYS OF CARRYING OUT DIFFERENTIAL BINDING ANALYSIS FROM OLD MATERIAL ###

In this session we will go over through the DiffBind pipeline to perform some data exploratory analysis and carryout the differential binding analysis. We will then evaluate the results and save those in files.

DiffBind is an R Biocondutor package that is used for identifying sites that are differentially enriched between sample groups. It works primarily with sets of peak calls ('peaksets'), which are sets of genomic intervals representing candidate protein binding sites for each sample. It includes functions that support the processing of peaksets, including overlapping and merging peak sets across an entire dataset, counting sequencing reads in overlapping intervals in peak sets, and identifying statistically significantly differentially bound sites based on evidence of binding affinity (measured by differences in read densities). We will discuss the importance of each step but for more detailed information please take a look at the DiffBind vignette ###LINK###.

## Setting up
You will need to have all of the required input files to follow this lesson. The instructions for this can be found in the "Setting up" section.

## Reading in Peaksets
DiffBind pipeline starts with reading in a peaksets derived either from peak callers such as MACS or using other criterion like genomic windows. We need to generate a metadata file with one line for each peakset and in the format that is compatible with DiffBind. 

Note: If multiple peak callers are used for comparison purposes each sample would have more than one line in the sample sheet. A merging function will generate a consensus peakset for the experiment. 

Once the peaksets are read in, a merging function finds all overlapping peaks and derives a single set of unique genomic intervals covering all the supplied peaks (a consensus peakset for the experiment). A region is considered for the consensus set if it appears in more than two of the samples. This consensus set represents the overall set of candidate binding sites to be used in further analysis.


## Metadata format
We have metadata prepared lets read it see the format DiffBind requires. 


```{r}
samples <- read.csv("data/metadata.csv")
view(samples)
names(samples)
```

An experiment will have multiple samples. Each sample needs a unique SampleID. A comparative analysis requires at leaset two samples in a class. Classes are indicated in the metadata as **Factor**, **Tissue**, **Condition** and/or **Treatment**

- Tissue - This is a designation for the cell type, tissue type, or some other indication of the biological source of the material.
  
- Factor - This is usually what protein the antibody was targeting, such as a transcription factor or a histone mark.

- Condition - Indicates an experimental condition, such as WT or Mutant.
  
- Treatment - Can denote how the cells were treated.

Not all of these classes are required, but there should be at leaset one class to compare. eg in our datasheet we have Tissue, Factor, and Condition.

Metadata should also comprise (point to) aligned sequencing reads (generally in bam format). Each sample needs an aligned sequencing library, but may include the following:

- bamReads: This points to the primary aligned file for the sample from Chip experiment (or other assays like ATAC)
- bamControl: This is an optional set of control reads associated with the sample or sample class. For ChIP experiments, this is most often an Input control (ChIP run without an antibody), or a ChIP run with a non-specific antibody. ATAC experiment usually do not have a control.
- SpikeIn: This is also an option set of spike-in reads for normalization.

Id for the input control is given in `ControlID` column.

The analysis presented here requires a called peaks for each sample. The aligned reads are used to call the peaks with the peak calling software such as MACS. It is also possible to perform a quantitative analysis without calling the peaks. For example, using windows that cover the entire genome like with `csaw` package or if the regions of interest are known in advance. The called peaks here are denoted in the `Peaks` column.

The first few steps in the pipeline can be computationally intensive, so this has been done for you and the output saved in the data folder as an `rds` object. Lets go over through these steps (Do NOT RUN)

```{r}
sample=read.csv("metadata.csv") # reading the metadata as before

dbObj=dba(sampleSheet=sample, scoreCol=5) # creating the diffbind object with dba function, scoreCol refers to the peak signal score in the peak files. If used MACS software scores are generally in 5th column

dbObj=dba.count(dbObj, bParallel=FALSE) # This is the most computationally intensive part


```

First step in the pipeline is to read in the metadata that we prepared. The second step involves creating the diffbind object. We do that using `dba` function. This function takes the metadata as `sampleSheet` and `scoreCol` referes to the peak singal score in the peak files. If used MACS software, the scores are generally in the 5th column.

After creating the dba object, next is to filter out the blacklists and greylists regions. Which are the problematic regions ###DEFINE PROBLEMATIC### in the genome. `Blacklist` regions for many reference genomes have been identified as part of the ENCODE project and can be accessed through the `dba.blacklist()` function in DiffBind. If the control samples are available, one can prepare regions to be excluded specific to the experiment called `Greylists` using the `GreyListChIP` package.

We have ran our data through nf-core pipeline, which employs blacklist removal in its pipeline so we are going to 


[old material](https://github.com/hbctraining/Peak_analysis_workshop/blob/main/lessons/OLD_diffbind_differential_peaks.md)
- Needs significant updating due to updates to diffbind
- Go through entire workflow; share dbObj with them because we will compute the affinity binding matrix on cluster (provide code and script)
- Some Viz such as MA plots, volcano and heatmaps


***

*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
