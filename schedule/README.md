# Workshop Schedule

> **NOTE:** The *Basic Data Skills* [Introduction to R](https://hbctraining.github.io/Intro-to-R-flipped/) workshop is a prerequisite. If you would like some practice wih R prior to taking this workshop, please work through this [R refresher lesson](https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/R_refresher.html).


## Pre-reading:

* Please **study the contents** and **work through all the exercises** within the following lessons:
  * [Workflow overview: From sequenced reads to peaks](../lessons/00a_peak_calling_workflow_review.md)
  * [Existing workflows for ChIP-seq analysis](../lessons/00b_peak_calling_with_nfcore.md)
  * Download the compressed R Project that we will be using by right-clicking on [this link](https://www.dropbox.com/scl/fi/s9mxwd7ttqgjt040m6bm2/Peak_analysis.zip?rlkey=ceqbv4pyx59jxsoa0xoh9l6kb&st=q7rlclil&dl=1) and selecting **"Save Link As..."** and download the compressed R Project to your desired location. Double-click on the compressed ZIP file in order to uncompress it.
 
  
## Day 1

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 09:30 - 09:45 | [Workshop Introduction](../lectures/workshop_intro_slides.pdf) | Will |
| 09:45 - 10:15 | Pre-reading discussion | Meeta |
| 10:15 - 11:00 | [Understanding peaks and peak file formats](../lessons/01_Introduction_to_peak_files.md)  | Will |
| 11:00- 11:05 | Break|  |
| 11:05 - 12:00 | [Assessing peak quality metrics](../lessons/02a_peak_quality_metrics_assesment.md) | Meeta |


### Before the next class:

I. Please **study the contents** and **work through all the code** within the following lessons:
   1. [Assessing sample similarity and identifying potential outliers](../lessons/02b_sample_similarity.md)
      <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br>One step in the QC of samples is to see how samples compare to one another. Generally, we expect replicates from each sample group to be more similar to each other and dissimilar to replicates from a different sample group. Here, we use read density (counts across the genome) and peak signal data to check if it meets our expectations.<br><br>In this lesson you will:<br>
             - Create PCA plots and inter-sample correlation heatmaps<br>
             - Evaluate plots to identify potential outliers and other effects<br>
             - Create visualiations using signal data from peaks to identify proposed thresholds for downstream analysis<br><br>
        </details>
   

   2. [Concordance across replicates using peak overlaps](../lessons/03_peak_overlaps.md)
      <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br>A quantitative way of evaluating how similar replicates are is to identify how many of the same peaks were called in each replicate. Biological replicates will inevitably exhibit some amount of variability, but the hope is that the majority of our peaks are identified in each sample. By looking at peak overlaps we can identify and remove a weaker replicate and/or use the overlap to create a consensus set of peaks. <br><br>In this lesson, we will:<br>
             - Discuss IRange and GRanges data structures in R<br>
             - Compute peak overlaps and create visualizations for the results<br><br>
        </details>
        

3. **Complete the exercises**:
   * Each lesson above contains exercises; please go through each of them.
   * Copy over your solutions into the [Google Forms](https://forms.gle/PMaZvtMWy92AhBEd7) the **day before the next class**.


### Questions?
* ***If you get stuck due to an error*** while runnning code in the lesson, [email us](mailto:hbctraining@hsph.harvard.edu) 

## Day 2

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 09:30 - 10:00 | Self-learning review | All |
| 10:00 - 10:45 | [Peak annotation and visualization using ChIPseeker](../lessons/04_peak_annotation_and_visualization.md)| Will |
| 10:45- 10:55 | Break|  |
| 10:55 - 12:00 | [Differential enrichment analysis using DiffBind](../lessons/05_diffbind_analysis.md) | Meeta |


### Before the next class:

I. Please **study the contents** and **work through all the code** within the following lessons:
   1. [Peak visualization using IGV](../lessons/06_peak_visualization_with_igv.md)
      <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br>Now that we have identified regions that are differentially enriched, it would be good to perform a qualitative assessment. To do this we will take a look at the data in IGV, a genome browser and see what read density looks like in significant regions.<br><br>In this lesson, we will:<br>
             - Learn how to navigate IGV and introduce various features<br>
             - Evaluate significant regions from DiffBind<br><br>
        </details>
   

   2. [Annotation and functional analysis of DE regions](../lessons/07_DE_annotation_and_enrichment_analysis.md)
      <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br>To gain biological insight from the genomic coordinates identified as differentially bound, we need to map them back to genomic features and see if there is some over-representation of target genes in specific pathways.  <br><br>In this lesson, we will:<br>
             - Use ChIPseeker to annotate the DE regions<br>
             - Perform functional analysis on the DE target genes<br><br>
        </details>
        
3. **Complete the exercises**:
   * The Functional Analysis lesson above contains exercises; please go through each of them.
   * Copy over your solutions into the [Google Forms](https://forms.gle/vivYK1LraGt2BUxu7) the **day before the next class**.

### Questions?
* ***If you get stuck due to an error*** while runnning code in the lesson, [email us](mailto:hbctraining@hsph.harvard.edu) 

***

## Day 3

| Time |  Topic  | Instructor |
|:-----------:|:----------:|:--------:|
| 09:30 - 10:30 | Self-learning review | All |
| 10:30 - 11:15 | [Motif analysis/discovery](../lessons/08_motif_analysis.md) | Meeta |
| 11:15- 11:25 | Break|  |
| 11:05 - 11:45 | Discussion Q&A | All |
| 11:45 - 12:00 | [Wrap-up](../lectures/Workshop_wrapup.pdf) | Will |


## Answer keys

* [Day 1 exercises](../homework/Day_1_answer_key.md)

* [Day 2 exercises](../homework/Day_2_answer_key.md)

## Resources

- [Comprehensive assessment of differential ChIP-seq tools guides optimal algorithm selection](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02686-y)
- [GRanges Tutorial](https://research.stowers.org/cws/CompGenomics/Tutorial/GRanges/guide.html)
- [ChIPpeakAnno Vignette](https://www.bioconductor.org/packages/devel/bioc/vignettes/ChIPpeakAnno/inst/doc/ChIPpeakAnno.html)
- [Introduction to Chromatin Biology Workshop](https://hbctraining.github.io/Intro-to-ChIPseq-flipped/)

****

*These materials have been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

***

*These materials have been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
