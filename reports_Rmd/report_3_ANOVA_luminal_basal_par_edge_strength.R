# anova analysis of parameters and edge strength between luminal and basal cell lines
library(plyr)
library(dplyr)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)


# load data
sampling_df_qc = readRDS("./data/models/pkn_v4_midas_v4/parameter_samplings_QC.RDS")
lumi_basal_annot = read.csv("./data/Cell_Lines_LuminalBasal_calling.csv",sep = ";")


mean_par_value = ddply(sampling_df_qc,.(cell_line),function(df){
	colMeans(df[,2:140])
})

best_par_value = ddply(sampling_df_qc,.(cell_line),function(df){
	df[which.min(df$fobj),2:140]
})

### 1. ANOVA analysis of parameter values between luminal and basal ---------------

# with mean data
aov_data = mean_par_value
aov_data = best_par_value

# with the parameter sample
#aov_data = sampling_df_qc[,2:141]

aov_data = merge(aov_data, lumi_basal_annot[,c("Cell.Line","Basal.Luminal")],by.x = "cell_line",by.y = "Cell.Line")
aov_data = melt(aov_data,id.vars = c("Basal.Luminal","cell_line"),variable.name = "parameter",value.name = "value",factorsAsStrings = T)


res = aov_data %>% group_by(parameter) %>%
	do(Model = aov(value ~ Basal.Luminal, data=.))


res_clear = broom::tidy(res,Model)
res_clear = as.data.frame(res_clear)
res_clear = filter(res_clear,term!="Residuals")

res_clear$p.value.adj = p.adjust(res_clear$p.value,method="BY")

signif_parameters = as.character(res_clear[res_clear$p.value<0.05,"parameter"])

ggplot(filter(aov_data,parameter %in% signif_parameters),aes(Basal.Luminal,value)) + geom_boxplot() +   facet_wrap(~parameter) + theme_bw()
ggplot(filter(aov_data,parameter %in% signif_parameters),aes(Basal.Luminal,value)) +
	geom_violin(draw_quantiles = c(.25,.5,.75)) +
	stat_summary(aes(group=Basal.Luminal), fun.y=mean, geom="point",
				 fill="black", shape=21, size=3, position = position_dodge(width = .9)) +
	facet_wrap(~parameter) + theme_bw()



## 1.2 post-hoc analysis of anova -----
poshoc_tHSD = lapply(res$Model,TukeyHSD)
pval_tHSD = ldply(poshoc_tHSD,function(test_res)test_res$Basal.Luminal[,4])
diff_tHSD = ldply(poshoc_tHSD,function(test_res)test_res$Basal.Luminal[,1])
rownames(pval_tHSD) = res$parameter
rownames(diff_tHSD) = res$parameter

pval_tHSD = pval_tHSD[rownames(pval_tHSD) %in% signif_parameters,]
diff_tHSD = diff_tHSD[rownames(diff_tHSD) %in% signif_parameters,]


gg = pheatmap(t(diff_tHSD),
		 display_numbers =  ifelse(t(pval_tHSD)>0.05,"",ifelse(t(pval_tHSD)>0.01,"*",ifelse(t(pval_tHSD)>0.001,"**","***"))),
		 cluster_rows = F,cluster_cols = F, main="ANOVA: parameter differences across luminal/basal/normal cell lines")
if(FALSE) ggsave(plot = gg$gtable,filename = "./supp_info/results_for_march_2019/figures/ANOVA_best_parameters_luminal_basal_normal.pdf",width = 8,height = 3.5)




### 2. ANOVA analysis of edge strength  between luminal and basal ---------------

rm(list=ls())
# load data
edge_strength = readRDS("./data/models/pkn_v4_midas_v4/features/feature_table_2_3_full_timecourse_edge_strength.RDS")
lumi_basal_annot = read.csv("./data/Cell_Lines_LuminalBasal_calling.csv",sep = ";")


aov_data = edge_strength
aov_data = merge(aov_data, lumi_basal_annot[,c("Cell.Line","Basal.Luminal")],by.x = "cell_line",by.y = "Cell.Line")
aov_data = melt(aov_data,id.vars = c("treatment","time","Basal.Luminal","cell_line"),measure.vars = 4:87 ,variable.name = "interaction",value.name = "value",factorsAsStrings = T)


res = aov_data %>% group_by(interaction) %>%
	do(Model = aov(value ~ Basal.Luminal, data=.))





res_clear = broom::tidy(res,Model)
res_clear = as.data.frame(res_clear)
res_clear = filter(res_clear,term!="Residuals")

res_clear$p.value.adj = p.adjust(res_clear$p.value,method="BY")

signif_parameters = as.character(res_clear[res_clear$p.value<0.05,"interaction"])

ggplot(filter(aov_data,parameter %in% signif_parameters),aes(Basal.Luminal,value)) + geom_boxplot() +   facet_wrap(~parameter) + theme_bw()
ggplot(filter(aov_data,parameter %in% signif_parameters),aes(Basal.Luminal,value)) +
	geom_violin(draw_quantiles = c(.25,.5,.75)) +
	stat_summary(aes(group=Basal.Luminal), fun.y=mean, geom="point",
				 fill="black", shape=21, size=3, position = position_dodge(width = .9)) +
	facet_wrap(~parameter) + theme_bw()



## 2.2 post-hoc analysis of anova -----
poshoc_tHSD = lapply(res$Model[res_clear$p.value.adj<0.05],TukeyHSD)
pval_tHSD = ldply(poshoc_tHSD,function(test_res)test_res$Basal.Luminal[,4])
diff_tHSD = ldply(poshoc_tHSD,function(test_res)test_res$Basal.Luminal[,1])
rownames(pval_tHSD) = res$interaction
rownames(diff_tHSD) = res$interaction

pval_tHSD = pval_tHSD[rownames(pval_tHSD) %in% signif_parameters,]
diff_tHSD = diff_tHSD[rownames(diff_tHSD) %in% signif_parameters,]


# filter for at least one difference is largerthan 0.2 :
max_dist = apply(diff_tHSD,1,function(x)abs(max(x)))

high_diff_tHSD = diff_tHSD[max_dist>0.1,]
high_pval_tHSD = pval_tHSD[max_dist>0.1,]

gg = pheatmap(t(diff_tHSD),
			  display_numbers =  ifelse(t(pval_tHSD)>0.05,"",ifelse(t(pval_tHSD)>0.01,"*",ifelse(t(pval_tHSD)>0.001,"**","***"))),
			  cluster_rows = F,cluster_cols = T, main="ANOVA: edge strength across luminal/basal/normal cell lines")
if(FALSE) ggsave(plot = gg$gtable,filename = "./supp_info/results_for_march_2019/figures/ANOVA_best_parameters_luminal_basal_normal.pdf",width = 8,height = 3.5)



gg = pheatmap(t(high_diff_tHSD),
			  display_numbers =  ifelse(t(high_pval_tHSD)>0.05,"",ifelse(t(high_pval_tHSD)>0.01,"*",ifelse(t(high_pval_tHSD)>0.001,"**","***"))),
			  cluster_rows = F,cluster_cols = T, main="ANOVA: edge strength across luminal/basal/normal cell lines")
if(FALSE) ggsave(plot = gg$gtable,filename = "./supp_info/results_for_march_2019/figures/ANOVA_edge_str_luminal_basal_normal.pdf",width = 8,height = 3.5)


## vulcano:

vulcano_data = merge(melt(cbind(interaction=rownames(pval_tHSD),pval_tHSD),variable.name = "Luminal.Basal",value.name = "p_value"),
melt(cbind(interaction=rownames(diff_tHSD),diff_tHSD),variable.name = "Luminal.Basal",value.name = "difference"))
ggplot(vulcano_data,aes(difference,-log10(p_value))) + geom_point(aes(col=interaction))


### 2.1 ANOVA analysis of parameter values for EGFR and HER2 ----------------------
lumi_basal_annot = read.csv("./data/Cell_Lines_LuminalBasal_calling.csv",sep = ";")
CNV_EGFR = readRDS("./data/CNV_EGFRandERRB2.rds")
CNV_EGFR = as.data.frame(CNV_EGFR)
CNV_EGFR$Cell.Line = rownames(CNV_EGFR)
clinical_annotation = merge(lumi_basal_annot,CNV_EGFR,all = T)

# with mean data
aov_data = best_par_value
aov_data$cell_line = rownames(best_par_value)

# with the parameter sample
aov_data = sampling_df_qc[,2:141]

aov_data = merge(aov_data, clinical_annotation[,c("Cell.Line","Basal.Luminal","EGFR","ERBB2")],by.x = "cell_line",by.y = "Cell.Line")
aov_data = melt(aov_data,id.vars = c("Basal.Luminal","EGFR","ERBB2","cell_line"),variable.name = "parameter",value.name = "value",factorsAsStrings = T)

aov_data$EGFR = factor(aov_data$EGFR)
aov_data$ERBB2 = factor(aov_data$ERBB2)


#res = aov_data %>% group_by(parameter) %>%
#	do(Model = aov(value ~ Basal.Luminal+ERBB2+EGFR, data=.))
res = aov_data %>% group_by(parameter) %>%
	do(Model = aov(value ~ EGFR, data=.))

res = aov_data %>% group_by(parameter) %>%
	do(Model = aov(value ~ ERBB2, data=.))

res_clear = broom::tidy(res,Model)
res_clear = as.data.frame(res_clear)
res_clear = filter(res_clear,term!="Residuals")

res_clear$p.value.adj = p.adjust(res_clear$p.value,method="BY")

res_clear = res_clear[order(res_clear$statistic,decreasing = T),]

ggplot(filter(aov_data,parameter %in% res_clear$parameter[1:12]),aes(EGFR,value)) +
	geom_boxplot() +
	stat_summary(aes(group=EGFR), fun.y=mean, geom="point",
				 fill="black", shape=21, size=3, position = position_dodge(width = .9)) +
	facet_wrap(~parameter) + theme_bw() + ggtitle("ANOVA by EGFR",subtitle = "top 12 parameters, highest F-statistics")# + scale_y_log10()


ggplot(filter(aov_data,parameter %in% res_clear$parameter[1:12]),aes(ERBB2,value)) +
	geom_boxplot() +
	stat_summary(aes(group=ERBB2), fun.y=mean, geom="point",
				 fill="black", shape=21, size=3, position = position_dodge(width = .9)) +
	facet_wrap(~parameter) + theme_bw() + ggtitle("ANOVA by ERBB2",subtitle = "top 12 parameters, highest F-statistics")# + scale_y_log10()


#res_clear = filter(res_clear,term=="Basal.Luminal")
# with p-value adjustment:
signif_parameters = as.character(res_clear[res_clear$p.value.adj<0.05,"parameter"])

# without p-value adjustment
signif_parameters = as.character(res_clear[res_clear$p.value<0.05,"parameter"])


ggplot(filter(aov_data,parameter %in% signif_parameters),aes(Basal.Luminal,value)) + geom_boxplot() +   facet_wrap(~parameter) + theme_bw()
ggplot(filter(aov_data,parameter %in% signif_parameters),aes(EGFR,value)) + geom_boxplot() +   facet_wrap(~parameter) + theme_bw()
ggplot(filter(aov_data,parameter %in% signif_parameters),aes(ERBB2,value)) + geom_boxplot() +   facet_wrap(~parameter) + theme_bw()


## pos hoc analysis of anova
poshoc_tHSD = lapply(res$Model,TukeyHSD)
pval_tHSD = ldply(poshoc_tHSD,function(test_res)test_res$Basal.Luminal[,4])
diff_tHSD = ldply(poshoc_tHSD,function(test_res)test_res$Basal.Luminal[,1])
rownames(pval_tHSD) = res$parameter
rownames(diff_tHSD) = res$parameter

pval_tHSD = pval_tHSD[rownames(pval_tHSD) %in% signif_parameters,]
diff_tHSD = diff_tHSD[rownames(diff_tHSD) %in% signif_parameters,]

pheatmap(t(diff_tHSD),
		 display_numbers =  ifelse(t(pval_tHSD)>0.05,"",ifelse(t(pval_tHSD)>0.01,"*",ifelse(t(pval_tHSD)>0.001,"**","***"))),cluster_rows = F,cluster_cols = F)


