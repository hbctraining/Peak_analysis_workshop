#!/bin/bash

#SBATCH -J deeptools
#SBATCH -t 0-01:00
#SBATCH -p priority
#SBATCH -c 6
#SBATCH	--mem 5G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

## Step-1:

#multiBamSummary bins \
#	  --bamfiles bam/cKO_H3K27ac_ChIPseq_REP1.mLb.clN.sorted.bam bam/cKO_H3K27ac_ChIPseq_REP2.mLb.clN.sorted.bam bam/cKO_H3K27ac_ChIPseq_REP3.mLb.clN.sorted.bam \
#               bam/cKO_H3K27ac_input_REP1.mLb.clN.sorted.bam bam/cKO_H3K27ac_input_REP2.mLb.clN.sorted.bam bam/cKO_H3K27ac_input_REP3.mLb.clN.sorted.bam \
#               bam/WT_H3K27ac_ChIPseq_REP1.mLb.clN.sorted.bam bam/WT_H3K27ac_ChIPseq_REP2.mLb.clN.sorted.bam bam/WT_H3K27ac_ChIPseq_REP3.mLb.clN.sorted.bam \
#               bam/WT_H3K27ac_input_REP1.mLb.clN.sorted.bam bam/WT_H3K27ac_input_REP2.mLb.clN.sorted.bam bam/WT_H3K27ac_input_REP3.mLb.clN.sorted.bam \
#	  --outFileName multiBamsummary.npz \
#	  --labels cKO_1 cKO_2 cKO_3 cKO_i1 cKO_i2 cKO_i3 WT_1 WT_2 WT_3 WT_i1 WT_i2 WT_i3 \
#	  -p 6 \
#	  --outRawCounts multiBAMsummary.tab


## Step-2

#plotCorrelation --corData multiBamsummary.npz \
#               --plotFile H3k27ac_pearson_corr_bin.png \
#               --outFileCorMatrix H3k27ac_pearson_corr_matrix_bin.txt \
#               --whatToPlot heatmap \
#               --corMethod pearson

## Step-3

#plotPCA -in multiBamsummary.npz -o H3K27Ac_deeptools_PCA.png
