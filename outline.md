Outline for peak analysis

**• Day 1**
**• Intro to File formats for ChIP-seq**
	• https://www.dropbox.com/scl/fi/04she5ab5j8l9zpa1ajr7/Workflows_and_fileformats.pdf?rlkey=hl5joq03cjg1xxpvv09a0gbc5&st=1wor80gl&dl=0
**• Peak concordance between replicates**
	• Should link to additional ref file or talk about multiqc output
	• PCA -- see QC_cutnrun_additional.Rmd
		○ (want to reassess/discuss best methods for this -- use peaks? VST/rlog? Top variable or all peaks?)
	• Correlation heatmap -- see QC_cutnrun_additional.Rmd
	• Rank vs enrichment  -- see QC_cutnrun_additional.Rmd
	• Histogram of quality scores -- see QC_cutnrun_additional.Rmd
**• Additional QC**
	• Peak overlap between replicates? -- see QC_cutnrun_additional.Rmd
		○ upsettR vs Venn
  
**• Day 2**
**• Peak annotation and visualization**
	• Intro to genomicRanges
	• Peak annotation -- see QC_cutnrun_additional.Rmd
	• Peak distance from TSS (or other genomic features) -- see QC_cutnrun_additional.Rmd
	• Consider including genometricCorr or other statistical evaluation of proximity of peaks to genomic regions of interest
		○ Need to check if this is commonly used or if there is another tool which might be best practice
	• ChIPSeeker plots -- see QC_cutnrun_additional.Rmd
		○ Mention Deeptools as an alternative using bam files, and can link to chromatin biology workshop materials
**• Differential enrichment analysis**
	• https://github.com/hbctraining/Peak_analysis_workshop/blob/main/lessons/08_diffbind_differential_peaks.md
		○ Needs significant updating due to updates to diffbind
**• Peak visualization using a genome viewer (IGV)**
	• https://github.com/hbctraining/Peak_analysis_workshop/blob/main/lessons/11_qualitative_assessment_IGV.md
	
**Additional reference file:**
-multiqc and other basic metrics to look out for![image](https://github.com/hbctraining/Peak_analysis_workshop/assets/33556230/b3e747d9-8a8f-41f7-9c73-05531362b49c)
