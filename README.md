# SpatialResolution
##### Author: Stephanie Zimmerman
This repository contains code used for analyses and visualizations for the paper Spatially resolved whole transcriptome profiling in human and mouse tissue using Digital Spatial Profiling published in Genome Research

Each folder contains code underlying the analysis of a figure in the manuscript. Data is found in the data/ folder and the data for each analysis is as follows:

Fig 1 - correlation to RNAseq/correlation_to_RNAseq.R and Fig 2 - sensitivity and specificity/sens+spec_RNAseq.R
human:
geomx_data_file <- "data/Human_11CPA_CollapsedTargetCounts.xlsx"
rnaseq_data_file <- "data/Human_updated_RNAseq_ccle_TPM_GM.zip"

mouse:
geomx_data_file <- "data/mouse_CPA_CollapsedTargetCounts.xlsx"
rnaseq_data_file <- "data/Mouse_Int_RNASeq_NSTG.csv"

Fig 2 - sensitivity and specificity/sens+spec_RNAscope.R
geomx_data_file <- "data/Human_23CPA_TargetCounts.xlsx"
rnascope_data_file <- "data/RNAscope_spot_counting.txt"
