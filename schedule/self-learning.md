# Introduction to Peak Analysis

## Learning Objectives

- Describe peak data and different file formats generated from peak calling algorithms
- Assess various metrics used to assess the quality of peak calls
- Compare peak calls across samples within a dataset
- Create visualizations to evaluate peak annotations
- Evaluate differentially enriched regions between two sample groups

## Installations

### On your desktop

1. [R](https://www.r-project.org/)
2. [RStudio](https://posit.co/download/rstudio-desktop/)
3. [Integrative Genomics Viewer (IGV)](https://software.broadinstitute.org/software/igv/)
4. [The listed R packages](../../Peak_analysis_workshop#installation-requirements)

### On your HPCC (if not using Harvard's O2 cluster)

#### Required
1. [`Nextflow`](https://www.nextflow.io/) version 24.11.0-edge

#### Alternative to Nextflow
1. [`samtools`](http://www.htslib.org) version 1.15.1
2. [`bedtools`](https://bedtools.readthedocs.io/en/latest/index.html) version 2.30.0
3. [`Picard`](https://broadinstitute.github.io/picard/) version 2.27.5
5. [`phantompeakqualtools`](https://github.com/kundajelab/phantompeakqualtools) version 1.2.2
6. [`deepTools`](https://deeptools.readthedocs.io/en/develop/index.html) version 3.5.6
7. [`bedGraphToBigWig`](https://github.com/ENCODE-DCC/kentUtils) version 302.1

> ***NOTE:*** If you are not working on the O2 cluster and are using different versions of these software programs, these packages may still work with the provided commands. However, this workshop was designed on these versions specifically, so you may need to tweak some of the commands if you use different versions of this software.

## Lessons

1. [Workflow overview: From sequenced reads to peaks](../lessons/00a_peak_calling_workflow_review.md)
2. [Existing workflows for ChIP-seq analysis](../lessons/00b_peak_calling_with_nfcore.md)
3. [Understanding peaks and peak file formats](../lessons/01_Introduction_to_peak_files.md)
4. [Assessing peak quality metrics](../lessons/02a_peak_quality_metrics_assesment.md)
5. [Assessing sample similarity and identifying potential outliers](../lessons/02b_sample_similarity.md)
6. [Concordance across replicates using peak overlaps](../lessons/03_peak_overlaps.md)
7. [Peak annotation and visualization using ChIPseeker](../lessons/04_peak_annotation_and_visualization.md)
8. [Differential enrichment analysis using DiffBind](../lessons/05_diffbind_analysis.md)
9. [Peak visualization using IGV](../lessons/06_peak_visualization_with_igv.md)
10. [Annotation and functional analysis of DE regions](../lessons/07_DE_annotation_and_enrichment_analysis.md)
11. [Motif analysis/discovery](../lessons/08_motif_analysis.md)


> ***NOTE:*** If you aren't working on Harvard's O2 cluster the directory structure for the HPCC that you are using is likely different and you will need to modify paths to work within your HPCC's directory structure.

## Answer key

- [Day 1 exercises](../homework/Day_1_answer_key.md)
- [Day 2 exercises](../homework/Day_2_answer_key.md))
***

*These materials have been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
