# cluster the cell lines based on parameters:
# - clustering parametrers
# - clustering based on parameter samples (mfu-sampler)
# - clustering with clinical data (luminal, basal, EGFR, HER2)


library(multiCellNOpt)
library(plyr)
library(dplyr)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)


## parameters of fitted models -------------------------------------------------
# this was first approach while waiting for MFusampling results, later we consider
# ensemble of model parameters
# load the parameters from the fitted models
# cluster the cells based on these parameters
calibrated_model_files = list.files("./data/models/pkn_v4_midas_v4/outputs/","*.RDS",full.names = T)

calibrated_models = lapply(calibrated_model_files,function(f){
	m = readRDS(f)
	# update the class definition to the newest version
	#m2 = logicODEModel$new(oldModel = m)
}
)

names(calibrated_models) = gsub(".RDS","",basename(calibrated_model_files))
# get parameters

estim_pars = ldply(calibrated_models,function(M){
	M$ode_parameters$parValues[M$ode_parameters$index_opt_pars]
},.id = "cell_line")

# plot raw parameters on heatmap
heatmap_data = estim_pars
rownames(heatmap_data) = heatmap_data$cell_line
heatmap_data$cell_line = NULL

if(FALSE){
	pheatmap::pheatmap(t(heatmap_data), cluster_rows =  F,colorRampPalette(brewer.pal(n = 7, name ="Reds"))(100), )
}
# ==> we dont see strong clustering, let see the distribution of the parameters over cell lines


estim_pars_long = melt(estim_pars,id.vars = "cell_line",value.name = "parameter_value",variable.name = "parameter")


if(FALSE){
	ggplot(filter(estim_pars_long,parameter=="PIP3_k_AKT_S473")) + geom_density(aes(x=parameter_value))

	ggplot(filter(estim_pars_long,parameter=="PIP3_k_AKT_S473"),aes(x=parameter,y=parameter_value)) + geom_violin() + geom_jitter()
	ggplot(estim_pars_long,aes(x=parameter,y=parameter_value)) + geom_violin() + geom_jitter() + facet_wrap(~parameter,scales = "free_x")

	ddply(estim_pars_long,.(parameter),summarise, sd_par=sd(parameter_value))
}

# checking the sampling results ------------------------------------------------

# 1. load results
if(FALSE){
	sampling_files = list.files("./data/models/pkn_v4_midas_v4/mfuSampling/","mfusampler_.*.RDS",full.names = T)

	samplings = lapply(sampling_files,function(f){
		m = readRDS(f)
		# update the class definition to the newest version
		#m2 = logicODEModel$new(oldModel = m)
	}
	)
	names(samplings) = gsub(".RDS","",basename(sampling_files))

	# concatenate
	sampling_df = ldply(samplings,function(df){
		return(df[!is.na(df[,1]),])
	},.id="sampling_id",.inform = T)


	M = calibrated_models[[1]]
	colnames(sampling_df) = c("sampling_id",  c(M$ode_parameters$parNames[M$ode_parameters$index_opt_pars],gsub("^1","x0",M$ode_parameters$x0Names)))
	sampling_df$cell_line = gsub("mfusampler_rundid_[0-9]+_","",sampling_df$sampling_id)
	sampling_df$run_id = gsub("_","",stringr::str_extract(sampling_df$sampling_id, "_[0-9]_"))


	# 2. add objective function
	sampling_df = adply(sampling_df,.margins = 1,function(x){

		fobj = calibrated_models[[x$cell_line]]$objectiveFunction(as.numeric(x[1,2:140]))
		return(fobj)
	},.progress = progress_text())


	colnames(sampling_df)[[143]] = "fobj"

	# add relative objective value
	sampling_df = ddply(sampling_df, .(cell_line),function(x) {
		x$rel_fobj = x$fobj/min(x$fobj)
		return(x)}
	)

	if(FALSE) saveRDS(sampling_df,"./data/models/pkn_v4_midas_v4/mfuSampling/00.estim_params_collected_v2.RDS")
}else{
	sampling_df = readRDS("./data/models/pkn_v4_midas_v4/mfuSampling/00.estim_params_collected_v2.RDS")
}

## plots
## 3. check the distribution of samples by objective function

# fobj distribution
ggplot(sampling_df,aes(cell_line,fobj)) + geom_jitter() +
	coord_flip() +
	ggtitle("Distribution of MfU sampler results for each cell lines") +
	xlab("cell line") +
	ylab("objective function value")+ theme_bw()

# relative fobj distribution
ggplot(sampling_df,aes(cell_line,rel_fobj)) + geom_jitter(size=.3) +
	coord_flip() +
	ggtitle("Distribution of MfU sampler results for each cell lines") +
	xlab("cell line") +
	ylab("relative objective function value")+ theme_bw()

if(FALSE) ggsave("./figures/cell_line_parameter_variance/samplings_rel_obj_function_distr.pdf",width = 8.5,height = 8.5)

# TODO: generate some sort of goodness of fit.


### Quality control

sampling_df_qc = sampling_df[sampling_df$rel_fobj<1.1,]

### cluster cell lines


### t-SNE plots ------------------------------------------------------------------

library(Rtsne)
library(ggrepel)

rtsne.out = Rtsne(sampling_df_qc[,2:140],verbose = TRUE,perplexity=100)

rtsne.gg = cbind(as.data.frame(rtsne.out$Y),"celline"=sampling_df_qc$cell_line)
rtsne.gg = filter(rtsne.gg,celline!="HCC70_2")
colnames(rtsne.gg)[1:2] = c("tSNE.X","tSNE.Y")
rtsne.gg$label = as.character(rtsne.gg$celline)
rtsne.gg$label[duplicated(rtsne.gg$label)] = ""

# t_SNE colored by cell line
ggplot(rtsne.gg,aes(tSNE.X,tSNE.Y,col=celline)) + geom_point()+
	geom_text(aes(label=label),color="black") + ggtitle("t-SNE based on parameter sampling") + theme_bw() +  guides(label=FALSE,color=FALSE)

if(FALSE) ggsave("./figures/cell_line_parameter_variance/tSNE_par_sampling_perplex_100.pdf",width = 10.5,height = 8.5)

# t-SNE with clinical annotation (luminal/basal)
lumi_basal_annot = read.csv("./data/Cell_Lines_LuminalBasal_calling.csv",sep = ";")

rtsne.gg = merge(rtsne.gg,lumi_basal_annot,by.x = "celline",by.y = "Cell.Line",all.x = T)


### Cluster feature and annotate with clinical data -- Done
ggplot(rtsne.gg,aes(tSNE.X,tSNE.Y,col=Basal.Luminal)) +
	geom_point()+geom_text(aes(label=label),color="black") +
	ggtitle("t-SNE based on parameter sampling") + theme_bw() +
	guides(label=FALSE)

if(FALSE)ggsave("./figures/cell_line_parameter_variance/tSNE_par_sampling_perplex_100_lumi_bas_normal.pdf",width = 10.5,height = 8.5)



### HEATMAPS --------------------------------------------------------------------
# cluster cellline parameters (mean/best fitting parameters) annotated by clinical data

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

mean_par_value = ddply(sampling_df_qc,.(cell_line),function(df){
	colMeans(df[,2:140])
})

best_par_value = ddply(sampling_df_qc,.(cell_line),function(df){
	df[which.min(df$fobj),2:140]
})




rownames(best_par_value) = best_par_value$cell_line
best_par_value$cell_line = NULL

#pheatmap::pheatmap(cellline_parameters,annotation_row = phm.annot.row)

out = pheatmap::pheatmap(best_par_value,
						 annotation_row = phm.annot.row,
						 colorRampPalette(brewer.pal(n = 7, name ="Reds"))(100),
						 fontsize_row =  6,
						 fontsize_col = 6,
						 cutree_rows = 9,
						 main="Clustering model parameters (based on best fitting parameters)",
						 annotation_colors = list(luminal.basal=c("Basal"="#E41A1C", "Luminal"="#377EB8", "Normal"="#4DAF4A"),
						 						 egfr=c(`-1`="#984EA3",`0`="#FF7F00",`1`="#FFFF33"),
						 						 HER2=c(`-1`="#A65628",`0`="#F781BF",`1`="#999999")))

if(FALSE) ggsave("./figures/cell_line_parameter_variance/best_parameter_heatmap_clinical.pdf",width = 12.5,height = 8.5,plot = out$gtable)


rownames(mean_par_value) = mean_par_value$cell_line
mean_par_value$cell_line = NULL


out = pheatmap::pheatmap(mean_par_value,
						 annotation_row = phm.annot.row,
						 colorRampPalette(brewer.pal(n = 7, name ="Reds"))(100),
						 fontsize_row =  6,
						 fontsize_col = 6,
						 cutree_rows = 9,
						 main="Clustering model parameters (based on mean parameters)",
						 annotation_colors = list(luminal.basal=c("Basal"="#E41A1C", "Luminal"="#377EB8", "Normal"="#4DAF4A"),
						 						 egfr=c(`-1`="#984EA3",`0`="#FF7F00",`1`="#FFFF33"),
						 						 HER2=c(`-1`="#A65628",`0`="#F781BF",`1`="#999999")))
if(FALSE) ggsave("./figures/cell_line_parameter_variance/mean_parameter_heatmap_clinical.pdf",width = 12.5,height = 8.5,plot = out$gtable)


# TODO: define distance based on t-score of samples?

# scaled_pars = log10(as.matrix(best_par_value+1e-8))
# scaled_pars = (scaled_pars - mean(scaled_pars))/sd(scaled_pars)
#
# out = pheatmap::pheatmap( scaled_pars,
# 						 annotation_row = phm.annot.row,
# 						 colorRampPalette(brewer.pal(n = 11, name ="RdYlBu"))(100),
# 						 fontsize_row =  6,
# 						 fontsize_col = 6,
# 						 main="Clustering model parameters (based on mean parameters)",
# 						 annotation_colors = list(luminal.basal=c("Basal"="#E41A1C", "Luminal"="#377EB8", "Normal"="#4DAF4A"),
# 						 						 egfr=c(`-1`="#984EA3",`0`="#FF7F00",`1`="#FFFF33"),
# 						 						 HER2=c(`-1`="#A65628",`0`="#F781BF",`1`="#999999")))

### PCA based on parameters and colored by clinical data -----------------------
pca_res = prcomp(sampling_df_qc[,2:140])

head(rtsne.gg)

rtsne.gg = cbind(rtsne.gg,pca_res$x[which(sampling_df_qc$cell_line!="HCC70_2"),1:2])

ggplot(rtsne.gg,aes(PC1,PC2,col=Basal.Luminal)) +
	geom_point(size=0.3)+ geom_text_repel(aes(label=label),color="black") +
	ggtitle("PCA based on parameter sampling") + theme_bw() +
	guides(label=FALSE)
if(FALSE)ggsave("./figures/cell_line_parameter_variance/PCA_par_sampling_lumi_bas_normal.pdf",width = 10.5,height = 8.5)



### ANOVA analysis of parameter values between luminal and basal ---------------

# with mean data
aov_data = mean_par_value
aov_data$cell_line = rownames(mean_par_value)

# with the parameter sample
aov_data = sampling_df_qc[,2:141]

aov_data = merge(aov_data, lumi_basal_annot[,c("Cell.Line","Basal.Luminal")],by.x = "cell_line",by.y = "Cell.Line")
aov_data = melt(aov_data,id.vars = c("Basal.Luminal","cell_line"),variable.name = "parameter",value.name = "value",factorsAsStrings = T)


res = aov_data %>% group_by(parameter) %>%
	do(Model = aov(value ~ Basal.Luminal, data=.))


res_clear = broom::tidy(res,Model)
res_clear = as.data.frame(res_clear)
res_clear = filter(res_clear,term!="Residuals")

res_clear$p.value.adj = p.adjust(res_clear$p.value,method="BY")

signif_parameters = as.character(res_clear[res_clear$p.value<0.05,"parameter"])

ggplot(filter(aov_data,parameter %in% signif_parameters[[1]]),aes(Basal.Luminal,value)) + geom_boxplot() +   facet_wrap(~parameter) + theme_bw()
ggplot(filter(aov_data,parameter %in% signif_parameters[[1]]),aes(Basal.Luminal,value)) +
	geom_violin(draw_quantiles = c(.25,.5,.75)) +
	stat_summary(aes(group=Basal.Luminal), fun.y=mean, geom="point",
				 fill="black", shape=21, size=3, position = position_dodge(width = .9)) +
	facet_wrap(~parameter) + theme_bw()


# TOP 12
res_clear = res_clear[order(res_clear$statistic,decreasing = T),]

ggplot(filter(aov_data,parameter %in% res_clear$parameter[1:12]),aes(Basal.Luminal,value)) +
	geom_boxplot() +
	stat_summary(aes(group=Basal.Luminal), fun.y=mean, geom="point",
				 fill="black", shape=21, size=3, position = position_dodge(width = .9)) +
	facet_wrap(~parameter) + theme_bw()+ ggtitle("ANOVA by Basal/Luminal",subtitle = "top 12 parameters, highest F-statistics")# + scale_y_log10()


gplots::plotmeans(filter(aov_data,parameter=="GSK3B_k_PIP3")$value~filter(aov_data,parameter=="GSK3B_k_PIP3")$Basal.Luminal)

# show significance on plot usubg ggpubr
# my_comparisons <- list( c("Basal", "Luminal"), c("Luminal", "Normal"), c("Basal", "Normal") )
# ggboxplot(filter(aov_data,parameter=="GSK3B_k_PIP3"), x = "Basal.Luminal", y = "value",
# 		  color = "Basal.Luminal", palette = "npg")+
# 	# Add pairwise comparisons p-value
# 	stat_compare_means(comparisons = my_comparisons, label.y = c(10, 12, 14))+
# 	stat_compare_means(label.y = 18)



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


### ANOVA analysis of parameter values for EGFR and HER2 ----------------------
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







## 4. check distribution of parametrs ---

#parameter_df_wide = reshape2::dcast(sampling_df,sampling_id+cell_line+run_id~parameter+fobj+rel_fobj,value.var = "parameter_value")
parameter_df_long = reshape2::melt(sampling_df,id.var = c("sampling_id","cell_line","run_id","fobj","rel_fobj"),value.name = "parameter_value",variable.name = "parameter")
head(parameter_df)




ggplot(parameter_df, aes(parameter,parameter_value)) + geom_jitter() + facet_wrap(~parameter,scales = "free_x")

ggplot(filter(parameter_df,parameter=="GSK3B_k_PIP3")) + geom_violin(aes(cell_line,parameter_value))+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("GSK3B_k_PIP3")


dir.create("./figures/cell_line_parameter_variance/")
pdf("./figures/cell_line_parameter_variance/par_distribution_violins.pdf",width = 15,height = 4)
parameters = as.character(unique(parameter_df$parameter))
# i =24
pb = progress::progress_bar$new(total = length(parameters))
for(i in 1: length(parameters)){
	pb$tick()
	gg = ggplot(filter(parameter_df,parameter==parameters[[i]])) +
		geom_violin(aes(cell_line,parameter_value))+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle(parameters[[i]])
	print(gg)
}

dev.off()

# parameter_stats: contains the statistics of parameters for each cell lines
parameter_stats = ddply(parameter_df, .(cell_line, parameter),function(df){
	data.frame(mean_value=mean(df$parameter_value),
			   median_value = median(df$parameter_value),
			   sd_value = sd(df$parameter_value),
			   interquartile.top = quantile(df$parameter_value,probs = 0.75),
			   interquartile.bottom = quantile(df$parameter_value,probs = 0.25)
	)
})
# parameter_stats_global: parameter statistics computed over all cell lines.
parameter_stats_global = ddply(parameter_df, .(parameter),function(df){
	data.frame(mean_value=mean(df$parameter_value),
			   median_value = median(df$parameter_value),
			   sd_value = sd(df$parameter_value),
			   interquartile.top = quantile(df$parameter_value,probs = 0.75),
			   interquartile.bottom = quantile(df$parameter_value,probs = 0.25)
	)
})





### check difference between clinical groups.



# apply clinical data colouring on t-SNE


