# Analysis of parameter variance for QTL analysis
# rank parameters by in-cell-line variance
# rank parameters by inter-cell-line variance of the median
# biaxial plot shows parameters to use.



library(multiCellNOpt)
library(plyr)
library(dplyr)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)


sampling_df = readRDS("./data/models/pkn_v4_midas_v4/mfuSampling/00.estim_params_collected_v2.RDS")

### Quality control
sampling_df_qc = sampling_df[sampling_df$rel_fobj<1.1,]

# save data for marco:
if(FALSE) saveRDS(sampling_df_qc,"./data/models/pkn_v4_midas_v4/parameter_samplings_QC.RDS")

sampling_df_qc_m = melt(sampling_df_qc,measure.vars = c(2:140),variable.name = "parameter" )


# par_stats: data.frame containing some statistics for each parameter estimated in each cell-line:
# fields:
#	- CL_mean: mean estimated value of parameter for each cell-line
#	- CL_median: median estimated value of the parameter for each cell-line
#	- CL_sd: standard deviation of the parameter value for each cell-line
#   - CL_cov: sd/mean aka coefficient of variation of the parameter in each cell line
par_stats = ddply(sampling_df_qc_m,.(cell_line,parameter),function(df){
	data.frame(CL_mean=mean(df$value),
			   CL_median=median(df$value),
			   CL_sd=sd(df$value))
})
par_stats$CL_cov =par_stats$CL_sd/par_stats$CL_mean

# save parameter stats for Marco:
if(FALSE) saveRDS(par_stats,"./data/models/pkn_v4_midas_v4/par_stats.RDS")


### Parameter estimability -----------------------------------------------------
# Rank the parameters based their estimability: how well their values were estimated
# in each cell-line the accuracy of the parameter estimate is measured by cov value. (par_stats$cov)
# we check for each parameter, how this par_stats$cov changes across cell-lines.
# an overall small value means that the parameter can be estimated in many cell-lines.

# par_global_estimability (of the coeff of var.): how the coef of variation changes across cell-lines, i.e. new  stats on is it's distribution.

par_global_estimability = ddply(par_stats,.(parameter),summarise,
								median_CL_cov=median(CL_cov),
								CL_cov_1q=quantile(CL_cov,0.05),
								CL_cov_3q=quantile(CL_cov,0.95),
								CL_cov_range =quantile(CL_cov,0.95)-quantile(CL_cov,0.05) )
# RANKING:
# order in INCREASING order: on the top of the list are the most estimables
# DONT CHANGE order or downstream breaks!!!!!
par_global_estimability = par_global_estimability[order(par_global_estimability$median_CL_cov,decreasing = F),]
par_global_estimability$parameter = factor(par_global_estimability$parameter,levels = par_global_estimability$parameter)


# plot the values in different forms:
ggplot(rbind(head(par_global_estimability,20),tail(par_global_estimability,20)),aes(x=parameter,y=median_CL_cov)) + geom_point()+
	geom_errorbar(aes(ymin=CL_cov_1q,ymax=CL_cov_3q)) + coord_flip() + theme_bw() +
	ggtitle("Parameter's estimation accuracy (top/bottom 20/20)",subtitle = "median and .9 quantile range \nestimated across all cell lines") +
	xlab("parameters") + ylab("coefficient of variation")

ggplot(filter(sampling_df_qc_m,parameter %in% c(as.character(par_global_estimability$parameter[c(1:2,138:139)]),"MEK12_k_ERK12","PKC_k_NFkB"))) +
	geom_boxplot(aes(cell_line,value)) + facet_wrap(~parameter,nrow=1) + coord_flip() + theme_bw() +
	xlab("cell line") + ylab('est. parameter values') + ggtitle("Top and least estimable parameters")

sampling_df_qc_m$parameter = factor(sampling_df_qc_m$parameter,levels = par_global_estimability$parameter)

ggplot(filter(sampling_df_qc_m,parameter %in% c(as.character(par_global_estimability$parameter[c(1:20,119:139)])))) +
	geom_violin(aes(parameter,value),scale = "width",trim = T) + coord_flip() + theme_bw() +
	xlab("cell line") + ylab('est. parameter values') + ggtitle("Top and least estimable parameters")


### Parameter variance across cell-lines ---------------------------------------
# rank the parameters how strongly they vary between cell-lines.
# large difference between cell-lines probably indicates different signaling routes.


par_global_variation = ddply(par_stats,.(parameter),summarise,
							 intCL_median=median(CL_median),  # median of the cellline medians: estimated values over all cell line
							 intCL_sd=sd(CL_median),   # standard deviation fo the cell-line medians
							 intCL_1q=quantile(CL_median,0.05),
							 intCL_3q=quantile(CL_median,0.95),
							 intCL_range=quantile(CL_median,0.95)-quantile(CL_median,0.05))  # range of the estimated values

# across cell variability RANKING
# On the top are the parameteres, that varies the most across cell-lines
# DONT CHANGE order or downstream breaks!!!!!
par_global_variation = par_global_variation[order(par_global_variation$intCL_range,decreasing = T),]
par_global_variation$parameter = factor(par_global_variation$parameter,levels = par_global_variation$parameter)


ggplot(par_global_variation,aes(x=parameter,y=intCL_median)) +
	geom_errorbar(aes(ymin=intCL_1q,ymax=intCL_3q)) + coord_flip() + theme_bw() +
	ggtitle("Parameter estimate changing across cell lines",subtitle = ".9 quantile range \nestimated across all cell lines") +
	xlab("parameters") + ylab("range of parameter value across cell lines")

ggplot(rbind(head(par_global_variation,30),tail(par_global_variation,30)),aes(x=parameter,y=intCL_median)) +
	geom_errorbar(aes(ymin=intCL_1q,ymax=intCL_3q)) + coord_flip() + theme_bw() +
	ggtitle("Parameter estimate changing across cell lines",subtitle = ".9 quantile range \nestimated across all cell lines") +
	xlab("parameters") + ylab("range of parameter value across cell lines")


parameters = as.character(par_global_estimability$parameter)

if(FALSE)  {
	saveRDS(par_global_estimability,"./data/models/pkn_v4_midas_v4/parameter_variance/par_global_estimability.RDS")
	saveRDS(par_global_variation,"./data/models/pkn_v4_midas_v4/parameter_variance/par_global_variation.RDS")
}

estimability_rank = match(parameters,par_global_estimability$parameter )
variability_rank = match(parameters,par_global_variation$parameter )

ranking = data.frame(parameter=parameters,
					 estimability_rank=estimability_rank,
					 variability_rank=variability_rank)

ggplot(ranking,aes(estimability_rank,variability_rank,fill=estimability_rank+variability_rank)) +
	geom_point() + geom_label_repel(aes(label=parameter),size = 2) +
	xlab("estimability rank (lower better)") + ylab("variability rank (lower the stronger)")+
	scale_fill_gradient(low='light blue', high='grey50') + guides(fill="none")



score_based_list = merge(par_global_estimability[,c("parameter","median_CL_cov")],par_global_variation[,c("parameter","intCL_range")])

gg_par = ggplot(score_based_list,aes(median_CL_cov,intCL_range,fill=median_CL_cov-intCL_range)) +  #
	geom_point() + geom_label_repel(aes(label=parameter),size = 2) +
	xlab("estimability score (lower better)") + ylab("inter-cellline variability score (larger, stronger)")+
	scale_fill_gradient(low='orange', high='grey90') + guides(fill="none") + theme_bw()

if(FALSE)  ggsave(plot =gg_par,filename = "./figures/cell_line_parameter_variance/par_estimability_vs_variance.pdf",width = 7,height = 7)

score_based_list$parameter_type = NA
score_based_list$parameter_type[grep("_k_",score_based_list$parameter)] = "k"
score_based_list$parameter_type[grep("x0",score_based_list$parameter)] = "x0"
score_based_list$parameter_type[grep("tau_",score_based_list$parameter)] = "tau"

gg_per_par = ggplot(score_based_list,aes(median_CL_cov,intCL_range,color=median_CL_cov-intCL_range)) +  #
	geom_point() + geom_text_repel(aes(label=parameter),size = 2) +
	xlab("estimability score (lower better)") + ylab("inter-cellline variability score (larger, stronger)")+
	scale_color_gradient(low='firebrick', high='grey20') + guides(color="none") + theme_bw() + facet_wrap(~parameter_type)

if(FALSE)  ggsave(plot =gg_per_par,filename = "./figures/cell_line_parameter_variance/par_estimability_vs_variance_perpartype.pdf",width = 13,height = 6)

