## Use the heatmap, but taking into account only the top estimable parameters

# run:
# knitr::opts_chunk$set(dpi=300,fig.width=7)
# knitr::spin("report_1_2_heatmaps_estimable_pars_clinic_data_annot.R")
# Then:
# knitr::spin("report_1_2_heatmaps_estimable_pars_clinic_data_annot.R",format = "Rhtml")



library(reshape2)
library(pheatmap)
library(plyr)
library(dplyr)

par_global_estimability = readRDS("./data/models/pkn_v4_midas_v4/parameter_variance/par_global_estimability.RDS")
par_global_variation = readRDS("./data/models/pkn_v4_midas_v4/parameter_variance/par_global_variation.RDS")
par_stats = readRDS("./data/models/pkn_v4_midas_v4/par_stats.RDS")


estimable_parameters = as.character(par_global_estimability$parameter[69:139])

mean_mat_2 = dcast(filter(par_stats,parameter %in% estimable_parameters),cell_line~parameter,value.var = "CL_median")

hmdata_mean_mat2 = mean_mat_2
rownames(hmdata_mean_mat2) = hmdata_mean_mat2$cell_line
hmdata_mean_mat2$cell_line = NULL

pheatmap(hmdata_mean_mat2,colorRampPalette(brewer.pal(n = 9, name ="Reds"))(100),cluster_cols = F)



lumi_basal_annot = read.csv("./data/Cell_Lines_LuminalBasal_calling.csv",sep = ";")
CNV_EGFR = readRDS("./data/CNV_EGFRandERRB2.rds")
CNV_EGFR = as.data.frame(CNV_EGFR)
CNV_EGFR$Cell.Line = rownames(CNV_EGFR)

clinical_annotation = merge(lumi_basal_annot,CNV_EGFR,all = T)
rownames(clinical_annotation) = clinical_annotation$Cell.Line

phm.annot.row = data.frame(luminal.basal=factor(clinical_annotation$Basal.Luminal),
						   egfr=factor(clinical_annotation$EGFR),
						   HER2=factor(clinical_annotation$ERBB2))
rownames(phm.annot.row) = clinical_annotation$Cell.Line

#pheatmap::pheatmap(cellline_parameters,annotation_row = phm.annot.row)

out1 = pheatmap::pheatmap(hmdata_mean_mat2,
						  annotation_row = phm.annot.row,
						  colorRampPalette(brewer.pal(n = 7, name ="Reds"))(100),
						  fontsize_row =  6,
						  fontsize_col = 6,
						  cluster_cols = F,
						  cutree_rows = 8,
						  clustering_distance_rows = "manhattan",
						  main="Clustering model parameters ",
						  annotation_colors = list(luminal.basal=c("Basal"="#E41A1C", "Luminal"="#377EB8", "Normal"="#4DAF4A"),
						  						 egfr=c(`-1`="#984EA3",`0`="#FF7F00",`1`="#FFFF33"),
						  						 HER2=c(`-1`="#A65628",`0`="#F781BF",`1`="#999999")))


# clustering feature with clinical data

edge_str = readRDS("./data/models/pkn_v4_midas_v4/features/feature_table_2_3_full_timecourse_edge_strength.RDS")
load("./data/models/pkn_v4_midas_v4/features/feature_table_2_edge_mean_activity.RData")

iegfr_13 = filter(edge_str,treatment=="iEGFR",time=="13")
rownames(iegfr_13) = iegfr_13$cell_line
iegfr_13$treatment = NULL
iegfr_13$time = NULL
iegfr_13$cell_line = NULL

out2 = pheatmap::pheatmap(iegfr_13,
						  annotation_row = phm.annot.row,
						  colorRampPalette(brewer.pal(n = 7, name ="Reds"))(100),
						  fontsize_row =  6,
						  fontsize_col = 6,
						  cluster_cols = F,
						  cutree_rows = 8,
						  clustering_distance_rows = "manhattan",
						  main="Clustering edge strength (iEGFR @13)",
						  annotation_colors = list(luminal.basal=c("Basal"="#E41A1C", "Luminal"="#377EB8", "Normal"="#4DAF4A"),
						  						 egfr=c(`-1`="#984EA3",`0`="#FF7F00",`1`="#FFFF33"),
						  						 HER2=c(`-1`="#A65628",`0`="#F781BF",`1`="#999999")))


rownames(feature_table_2) = feature_table_2$cell_line
feature_table_2$cell_line = NULL

out2 = pheatmap::pheatmap(feature_table_2,
						  annotation_row = phm.annot.row,
						  colorRampPalette(brewer.pal(n = 7, name ="Reds"))(100),
						  fontsize_row =  6,
						  fontsize_col = 6,
						  cluster_cols = F,
						  cutree_rows = 8,
						  clustering_distance_rows = "euclidean",
						  main="Clustering mean edge strength",
						  annotation_colors = list(luminal.basal=c("Basal"="#E41A1C", "Luminal"="#377EB8", "Normal"="#4DAF4A"),
						  						 egfr=c(`-1`="#984EA3",`0`="#FF7F00",`1`="#FFFF33"),
						  						 HER2=c(`-1`="#A65628",`0`="#F781BF",`1`="#999999")))
