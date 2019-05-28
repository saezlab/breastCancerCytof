# Goodness of fit again: compare model error to biological reproducibility
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(purrr)
library(multiCellNOpt)
# we compare the RMSE of the proteins versus the biological differences between A/B time courses.



######## 1. EXTRACT DATA STATISTICS ############################################

if(FALSE){
	# loads the data from the data_preparation pipeline and interpolates the A/B time courses, adds standard deviation
	celline_data_melted_ext  = readRDS("./data/all_cellline_scaled_AB_timecourse_data_RDS")

	celline_data_melted_ext = filter(celline_data_melted_ext, treatment != "full")

	base_time = c(0,5.5,7,9,13,17,23,30,40,60)
	celline_data_melted_ext$time = as.numeric(as.character(celline_data_melted_ext$time))

	# we need to interpolate the time course A and B

	celline_data_melted_ext_interp = ddply(celline_data_melted_ext,.(cell_line,treatment,markers),function(df){
		# df = filter(celline_data_melted_ext,cell_line=="HCC1428",treatment=="iPI3K",markers=="p-S6")

		df = df[order(df$time),]
		timecourse_A = filter(df,time_course=="A")
		timecourse_B = filter(df,time_course=="B")

		y_interp_A = approx(timecourse_A$time,timecourse_A$scaled_signal,base_time, method = "linear",rule=2)
		y_interp_B = approx(timecourse_B$time,timecourse_B$scaled_signal,base_time, method = "linear",rule=2)


		interp_data = data.frame(time=base_time,signal_A=y_interp_A$y,signal_B=y_interp_B$y)

	},.progress = progress_text())


	head(celline_data_melted_ext_interp)
	celline_data_melted_ext_interp$AB_mean_signal = 0.5*(celline_data_melted_ext_interp$signal_A + celline_data_melted_ext_interp$signal_B)

	# we compute the sample standard deviation (n = 2)
	celline_data_melted_ext_interp$sd_signal = sqrt((celline_data_melted_ext_interp$signal_A - celline_data_melted_ext_interp$AB_mean_signal)^2 +
														(celline_data_melted_ext_interp$signal_B - celline_data_melted_ext_interp$AB_mean_signal)^2)
	# plots
	hist(celline_data_melted_ext_interp$sd_signal,100)
	plot(celline_data_melted_ext_interp$signal_A,celline_data_melted_ext_interp$signal_B)

	if(FALSE) saveRDS(celline_data_melted_ext_interp,"supp_info/results_for_march_2019/data/celline_data_melted_ext_interp.RDS")

}else{
	celline_data_melted_ext_interp = readRDS("supp_info/results_for_march_2019/data/celline_data_melted_ext_interp.RDS")
}


## Compare RMSE difference between time course A and B for each cell line and each protein

# A vs B
data_rmse_r2_by_prots = ddply(celline_data_melted_ext_interp,.(cell_line,markers),summarize,
	  						 rmse = sqrt(sum((signal_A-signal_B)^2)/length(signal_A)),
	  						 r2 = cor(signal_A,signal_B)^2,
							 r = cor(signal_A,signal_B))

#A/B vs mean
data_rmse_r2_by_prots = ddply(celline_data_melted_ext_interp,.(cell_line,markers),summarize,
							  rmse = sqrt(
							  	(sum((signal_A-AB_mean_signal)^2) + sum((signal_B-AB_mean_signal)^2))/(2*length(signal_A))
							  	),
							  r2 = cor(signal_A,signal_B)^2,
							  r = cor(signal_A,signal_B))





######## 2. EXTRACT MODEL STATISTICS ###########################################

# use the best parameters or use the best parameter vector from the initial fit?
use_bestPar = FALSE

calibrated_model_files = list.files("./data/models/pkn_v4_midas_v4/outputs/","*.RDS",full.names = T)
calibrated_models = lapply(calibrated_model_files,function(f){
	m = readRDS(f)
	# update the class definition to the newest version
	#m2 = logicODEModel$new(oldModel = m)
}
)
names(calibrated_models) = gsub(".RDS","",basename(calibrated_model_files))

sampling_df = readRDS("./data/models/pkn_v4_midas_v4/mfuSampling/00.estim_params_collected_v2.RDS")
sampling_df_qc = sampling_df[sampling_df$rel_fobj<1.1,]

# we report the statistics only on the best model from each cell-line



get_annot_df_from_cnolist = function(signal_list,data_list,timepoints,experiments){

	treatment = rep(experiments, ncol(signal_list[[1]])*length(signal_list))
	time = rep(timepoints, each=ncol(signal_list[[1]])*nrow(signal_list[[1]]))
	markers = rep(rep(colnames(signal_list[[1]]), each=nrow(signal_list[[1]])),length(signal_list))
	ggdata = data.frame(sim=unlist(signal_list),dat=unlist(data_list),treatment,time,markers)
}


if(use_bestPar){
	best_parameters = ddply(sampling_df_qc,.(cell_line),function(df) df[which.min(df$fobj),])


	sim_vs_data_all = adply(best_parameters,.expand = F,.margins = 1,function(x){
		#x = best_parameters[1,]

		M = calibrated_models[[x$cell_line]]
		#M$plotFit()

		M$ode_parameters$parValues[M$ode_parameters$index_opt_pars] = as.numeric(x[1,M$ode_parameters$parNames[M$ode_parameters$index_opt_pars]])
		M$ode_parameters$x0Values[M$ode_parameters$index_opt_x0] = as.numeric(x[1,gsub("^1","x0",M$ode_parameters$x0Names[M$ode_parameters$index_opt_x0])])
		M$updateInitialValueEstimates(ode_parameters = M$ode_parameters,initialConditions = M$initialConditions)

		experiments = c("control",   "EGF", "iEGFR",  "iMEK12", "imTOR",  "iPI3K",  "iPKC")
		ggdata = get_annot_df_from_cnolist(signal_list =  M$simulate()$signals, data_list = M$exps$signals, timepoints=M$exps$timepoints, experiments = experiments)
		ggdata$cell_line = x$cell_line
		return(ggdata)
	},.progress = progress_text())
}else{
	sim_vs_data_all <- tibble(models = calibrated_models,cell_line=names(calibrated_models)) %>%
		mutate( simdata = map(models,function(M){

			experiments = c("control",   "EGF", "iEGFR",  "iMEK12", "imTOR",  "iPI3K",  "iPKC")
			ggdata = get_annot_df_from_cnolist(signal_list =  M$simulate()$signals, data_list = M$exps$signals, timepoints=M$exps$timepoints, experiments = experiments)
		} )) %>% unnest(simdata)

}


randModel_vs_data_all <- tibble(models = calibrated_models,cell_line=names(calibrated_models)) %>%
	mutate( simdata = map(models,function(M){
		# initialise a model with random parameters, but with the real data and SIF
		 Mrand = multiCellNOpt::logicODEModel$new(SIFfile = M$pkn,exps = M$exps)
		experiments = c("control",   "EGF", "iEGFR",  "iMEK12", "imTOR",  "iPI3K",  "iPKC")
		iter = 0
		while(TRUE){
			Mrand$initODE_parameters(random = T)
			ggdata = get_annot_df_from_cnolist(signal_list =  Mrand$simulate()$signals, data_list = Mrand$exps$signals, timepoints=Mrand$exps$timepoints, experiments = experiments)
			iter = iter + 1
			if(!all(is.na(ggdata$sim)) | iter > 10) break()
		}

		ggdata
	} )) %>% unnest(simdata)




# sim vs mean data

model_rmse_r2_by_prots = ddply(sim_vs_data_all,.(cell_line,markers),summarize,
						 rmse = sqrt(sum((sim-dat)^2)/length(sim)),
						 r2 = cor(sim,dat)^2,
						 r = cor(sim,dat),
						 data_sd = sd(dat))

random_model_rmse_r2_by_prots = ddply(randModel_vs_data_all,.(cell_line,markers),summarize,
							   rmse = sqrt(sum((sim-dat)^2)/length(sim)),
							   r2 = cor(sim,dat)^2,
							   r = cor(sim,dat),
							   data_sd = sd(dat))



######## 3. COMPARE MODEL AND DATA STATISTICS ##################################


model_rmse_r2_by_prots$type = "model_vs_data"
random_model_rmse_r2_by_prots$type = "random_model_vs_data"
data_rmse_r2_by_prots$type = "A/B_vs_mean"
data_rmse_r2_by_prots$markers = as.character(data_rmse_r2_by_prots$markers)
data_rmse_r2_by_prots$markers = gsub("p-","",data_rmse_r2_by_prots$markers)

colnames(model_rmse_r2_by_prots)[[2]] = "markers"
model_rmse_r2_by_prots$markers = as.character(model_rmse_r2_by_prots$markers)


combined_stats = rbind(data_rmse_r2_by_prots,model_rmse_r2_by_prots[,-6],random_model_rmse_r2_by_prots[,-6])


## 3.1 Based on RMSE ####
# compare: model vs data RMSE error and time_course_A vs time_course_B RMSE

# order cell lines by A_VS_B RMSE
tmp = filter(combined_stats,type=="A/B_vs_mean")
tmp = ddply(tmp,.(cell_line),summarise,median_rmse =median(rmse) )
tmp = tmp[order(tmp$median_rmse),]
combined_stats$cell_line = factor(combined_stats$cell_line,levels = tmp$cell_line)

ggplot(combined_stats,aes(x=cell_line,y=rmse,col=type,fill=type)) +
	geom_boxplot()  + coord_flip() + theme_bw() +
	ggtitle("Compare model and experimental reproducibility",subtitle = "model fitting error versus biological reproduction error") +
	xlab("Cell lines") + ylab("Root mean square error (RMSE)") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

if(FALSE) ggsave("./supp_info/results_for_may_2019/figures/RMSE_modeling_vs_biological_mean_exp_nogrid.pdf",width = 8,height = 7.5)


ggplot(combined_stats,aes(x=rmse,col=type,fill=type)) +
	geom_density(alpha=0.5)  +
	 xlab("Root mean square error (RMSE)") + theme_bw()
if(FALSE) ggsave("./supp_info/results_for_may_2019/RMSE_modeling_vs_biological_density.pdf",width = 8,height = 3)


# compare: model vs data RMSE error and time_course_A vs time_course_B RMSE -- overall

ggplot(combined_stats,aes(x=type,y=rmse,fill=type)) +
	geom_boxplot()  + coord_flip() + theme_bw() +
	ggtitle("Compare model and experimental reproducibility",subtitle = "model fitting error versus biological reproduction error") +
	xlab("Cell lines") + ylab("Root mean square error (RMSE)") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

if(FALSE) ggsave("./supp_info/results_for_may_2019/RMSE_randModel_modeling_vs_biological_mean_exp_nogrid.pdf",width = 8,height = 3)




## 3.1.1 Based on RMSE ####
# compare: plot the ratio of model_RMSE and data_RMSE
library(tidyverse)
data_rmse_by_prots <- data_rmse_r2_by_prots %>% as_tibble() %>% rename(data_rmse = rmse) %>% select(cell_line,markers,data_rmse)
model_rmse_by_prots <- model_rmse_r2_by_prots %>% as_tibble() %>% rename(model_rmse = rmse) %>% select(cell_line,markers,model_rmse)

RMSE_ratio  <- data_rmse_by_prots %>% full_join(x = .,y=model_rmse_by_prots,by=c("cell_line","markers"))

RMSE_ratio <- RMSE_ratio %>% mutate(rmse_ratio=model_rmse/data_rmse)

RMSE_ratio %>%
	ggplot(aes(x=cell_line,y=rmse_ratio)) +
	geom_boxplot()  + coord_flip() + theme_bw() +
	ggtitle("Model and Data RMSE ratio",subtitle = "modelling error compared to data reproducability") +
	xlab("Cell lines") + ylab("RMSE_model / RMSE_data") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_hline(yintercept = 1,color="navyblue")
if(FALSE) ggsave("./supp_info/results_for_march_2019/figures/RMSE_ratio_modeling_vs_biological_boxplot.pdf",width = 8,height = 7.5)

##### 3.1.2 Pick some example and plot ####
library(tidyr)

# merge the model time courses and the data
celline_data_melted_ext_interp$markers = as.character(celline_data_melted_ext_interp$markers)
celline_data_melted_ext_interp$markers = gsub("p-","",celline_data_melted_ext_interp$markers)

m_timecourse_data = merge(celline_data_melted_ext_interp,
						  filter(sim_vs_data_all,treatment!="control"),all = T)

get_min = function(x,y){


	df=data.frame(x,y)
	apply(df, 1, FUN=min)
}

get_max = function(x,y){
	df=data.frame(x,y)
	apply(df, 1, FUN=max)
}

# selects cell-line, and plots the model sim vs data
gen_RMSE_fig <- function(RMSE_cell_line,m_timecourse_data,rmse_stats){

	cl_timecourse_data = filter(m_timecourse_data,cell_line==RMSE_cell_line, markers==rmse_stats$markers)

	cl_timecourse_data %>% as_tibble() %>%
		filter(complete.cases(.))%>%
		ggplot() +
		#geom_ribbon(aes(time,ymin=get_min(signal_A,signal_B),ymax=get_max(signal_A,signal_B),fill=treatment), alpha=0.2) +
		geom_line(aes(time,sim,col=treatment))+
		geom_point(aes(time,dat,col=treatment,shape="mean_signal")) +
		geom_point(aes(time,signal_A,col=treatment,shape="signal_A")) +
		geom_point(aes(time,signal_B,col=treatment,shape="signal_B")) +
		geom_pointrange(aes(time,dat,ymin=get_min(signal_A,signal_B),ymax=get_max(signal_A,signal_B),col=treatment)) +
		# # geom_point(aes(time,AB_mean_signal,col=treatment),shape=2) +
		scale_shape_manual(values=c("signal_A" = 2,"signal_B" = 4,"mean_signal"=1))+
		theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
		ggtitle(paste0(rmse_stats$cell_line,", marker ", rmse_stats$markers,", model RMSE: ",round(rmse_stats$rmse[[1]],digits = 3),", data RMSE: ",round(rmse_stats$rmse[[2]],digits = 3) ))+
		ylab("Signal") + xlab("Time [mins]") +coord_cartesian(ylim=c(0,1))
}

# select cell_line sample 1
RMSE_cell_line = "ZR751"
combined_stats = as_tibble(combined_stats)
model_stat = combined_stats %>% filter(cell_line==RMSE_cell_line,type=="model_vs_data") %>%  top_n(1,rmse)
data_stat = combined_stats %>% filter(cell_line==RMSE_cell_line,type=="A/B_vs_mean",markers==model_stat$markers)
gg = gen_RMSE_fig(RMSE_cell_line,m_timecourse_data,rbind(model_stat,data_stat))
plot(gg)
if(FALSE) ggsave(plot = gg,filename = "./figures/for_paper/fitted_trajectory_sample_1_errorbar_v2.pdf",width = 6,height = 5)
#if(FALSE) ggsave(plot = gg,filename = "./figures/for_paper/fitted_trajectory_sample_1_with_ribbon.pdf",width = 6,height = 5)
#if(FALSE) ggsave(plot = gg,filename = "./figures/for_paper/fitted_trajectory_sample_1_with_errorbar.pdf",width = 6,height = 5)


# select cell_line sample 2
RMSE_cell_line = "AU565"
combined_stats = as_tibble(combined_stats)
model_stat = combined_stats %>% filter(cell_line==RMSE_cell_line,type=="model_vs_data") %>%  top_n(1,rmse)
data_stat = combined_stats %>% filter(cell_line==RMSE_cell_line,type=="A/B_vs_mean",markers==model_stat$markers)
gg = gen_RMSE_fig(RMSE_cell_line,m_timecourse_data,rbind(model_stat,data_stat))
plot(gg)
if(FALSE) ggsave(plot = gg,filename = "./figures/for_paper/fitted_trajectory_sample_2_errorbar_v2.pdf",width = 6,height = 5)
#if(FALSE) ggsave(plot = gg,filename = "./figures/for_paper/fitted_trajectory_sample_2.pdf",width = 6,height = 5)


# select cell_line sample 3 -- good example
RMSE_cell_line = "184A1"
model_stat = combined_stats %>% filter(cell_line==RMSE_cell_line,type=="model_vs_data",markers=="H3")
data_stat = combined_stats %>% filter(cell_line==RMSE_cell_line,type=="A/B_vs_mean",markers==model_stat$markers)
gg = gen_RMSE_fig(RMSE_cell_line,m_timecourse_data,rbind(model_stat,data_stat))
plot(gg)
if(FALSE) ggsave(plot = gg,filename = "./figures/for_paper/fitted_trajectory_sample_3_errorbar_v2.pdf",width = 6,height = 5)
#if(FALSE) ggsave(plot = gg,filename = "./figures/for_paper/fitted_trajectory_sample_3.pdf",width = 6,height = 5)

# select cell_line sample 4 -- good example
RMSE_cell_line = "BT549"
model_stat = combined_stats %>% filter(cell_line==RMSE_cell_line,type=="model_vs_data",markers=="MAPKAPK2")
data_stat = combined_stats %>% filter(cell_line==RMSE_cell_line,type=="A/B_vs_mean",markers==model_stat$markers)
gg = gen_RMSE_fig(RMSE_cell_line,m_timecourse_data,rbind(model_stat,data_stat))
plot(gg)
if(FALSE) ggsave(plot = gg,filename = "./figures/for_paper/fitted_trajectory_sample_4_errorbar_v2.pdf",width = 6,height = 5)
#if(FALSE) ggsave(plot = gg,filename = "./figures/for_paper/fitted_trajectory_sample_4.pdf",width = 6,height = 5)



# select cell_line sample 4 -- good example
RMSE_cell_line = "HCC1187"
model_stat = combined_stats %>% filter(cell_line==RMSE_cell_line,type=="model_vs_data",markers=="p38")
data_stat = combined_stats %>% filter(cell_line==RMSE_cell_line,type=="A/B_vs_mean",markers==model_stat$markers)
gg = gen_RMSE_fig(RMSE_cell_line,m_timecourse_data,rbind(model_stat,data_stat))
plot(gg)
if(FALSE) ggsave(plot = gg,filename = "./figures/for_paper/fitted_trajectory_sample_5_errorbar_v2.pdf",width = 6,height = 5)




### 3.2 Based on R2 ####
# compare: model vs data correlaton and time_course_A vs time_course_B correlation


# order cell lines by A_VS_B RMSE
tmp = filter(combined_stats,type=="A_vs_B")
tmp = ddply(tmp,.(cell_line),summarise,median_r2 =median(r2) )
tmp = tmp[order(tmp$median_r2),]
combined_stats$cell_line = factor(combined_stats$cell_line,levels = tmp$cell_line)


ggplot(combined_stats,aes(x=cell_line,y=r2,col=type,fill=type)) +
	geom_boxplot()  + coord_flip() + theme_bw() +
	ggtitle("Compare model and experimental reproducibility",subtitle = "model fitting error versus biological reproduction error") +
	xlab("Cell lines") + ylab("squared correlation coefficient (r2)")

if(FALSE) ggsave("./supp_info/results_for_march_2019/figures/R2_modeling_vs_biological.pdf",width = 8,height = 7.5)


ggplot(combined_stats,aes(x=r2,col=type,fill=type)) +
	geom_density(alpha=0.5)  +
	xlab("squared correlation coefficient (r2)") + theme_bw()
if(FALSE) ggsave("./supp_info/results_for_march_2019/figures/R2_modeling_vs_biological_density.pdf",width = 8,height = 3)


## 3.3 Based on R ####
# compare: model vs data correlaton and time_course_A vs time_course_B correlation


# order cell lines by A_VS_B RMSE
tmp = filter(combined_stats,type=="A_vs_B")
tmp = ddply(tmp,.(cell_line),summarise,median_r =median(r) )
tmp = tmp[order(tmp$median_r),]
combined_stats$cell_line = factor(combined_stats$cell_line,levels = tmp$cell_line)


ggplot(combined_stats,aes(x=cell_line,y=r,col=type,fill=type)) +
	geom_boxplot()  + coord_flip() + theme_bw() +
	ggtitle("Compare model and experimental reproducibility",subtitle = "model fitting error versus biological reproduction error") +
	xlab("Cell lines") + ylab("correlation coefficient ")

if(FALSE) ggsave("./supp_info/results_for_march_2019/figures/corr_modeling_vs_biological.pdf",width = 8,height = 7.5)


ggplot(combined_stats,aes(x=r,col=type,fill=type)) +
	geom_density(alpha=0.5)  +
	xlab("correlation coefficient") + theme_bw()
if(FALSE) ggsave("./supp_info/results_for_march_2019/figures/corr_modeling_vs_biological_density.pdf",width = 8,height = 3)

### 3.4  R2 (data_SD > noise) ####

# noise in the data --> null distribution
plot(hist(data_rmse_r2_by_prots$rmse))
noise_095 = quantile(data_rmse_r2_by_prots$rmse,0.95)

# find celllines/signals to be filtered out
#remove_cond = model_rmse_r2_by_prots[which(model_rmse_r2_by_prots$data_sd < noise_095),]
remove_cond = model_rmse_r2_by_prots[which(model_rmse_r2_by_prots$data_sd < noise_095),c("cell_line","markers")]
remove_cond = paste0(remove_cond$cell_line,"_",remove_cond$markers)

keep_data = which(! paste0(data_rmse_r2_by_prots$cell_line,"_",data_rmse_r2_by_prots$markers) %in% remove_cond)
keep_model = which(! paste0(model_rmse_r2_by_prots$cell_line,"_",model_rmse_r2_by_prots$markers) %in% remove_cond)

tmp = ddply(model_rmse_r2_by_prots[keep_data,],.(cell_line),summarise,median_r2 =median(r2) )
tmp = tmp[order(tmp$median_r2),]
model_rmse_r2_by_prots$cell_line = factor(model_rmse_r2_by_prots$cell_line,levels = tmp$cell_line)

cmd_dat = rbind(model_rmse_r2_by_prots[keep_model,-6],data_rmse_r2_by_prots[keep_data,])


ggplot(cmd_dat,aes(x=cell_line,y=r2,col=type,fill=type)) +
	geom_boxplot()  + coord_flip() + theme_bw() +
	ggtitle("Correlation analysis of changing proteins only",subtitle = "correlation between model trajectories and data per protein") +
	xlab("Cell lines") + ylab("correlation coefficient")

if(FALSE) ggsave("./supp_info/results_for_march_2019/figures/R2_modeling_boxplot_filtered_noise.pdf",width = 8,height = 7.5)

# ggplot(data_rmse_r2_by_prots[keep_data,],aes(x=cell_line,y=r2)) +
# 	geom_boxplot()  + coord_flip() + theme_bw() +
# 	ggtitle("Correlation analysis of chaning proteins only",subtitle = "correlation between model trajectories and data per protein") +
# 	xlab("Cell lines") + ylab("correlation coefficient")
#
# ggplot(model_rmse_r2_by_prots[keep_model,],aes(x=cell_line,y=r2)) +
# 	geom_boxplot()  + coord_flip() + theme_bw() +
# 	ggtitle("Correlation analysis of chaning proteins only",subtitle = "correlation between model trajectories and data per protein") +
# 	xlab("Cell lines") + ylab("correlation coefficient")


########### 4. Model Statistics #########################################

## 4.1 RMSE  #####
tmp = ddply(model_rmse_r2_by_prots,.(cell_line),summarise,median_rmse =median(rmse) )
tmp = tmp[order(tmp$median_rmse),]
model_rmse_r2_by_prots$cell_line = factor(model_rmse_r2_by_prots$cell_line,levels = tmp$cell_line)

ggplot(model_rmse_r2_by_prots,aes(x=cell_line,y=rmse)) +
	geom_boxplot()  + coord_flip() + theme_bw() +
	ggtitle("RMSE of fitted model",subtitle = "difference between model trajectories and data per protein") +
	xlab("Cell lines") + ylab("Root mean square error (RMSE)")

if(FALSE) ggsave("./supp_info/results_for_march_2019/figures/RMSE_modeling_boxplot.pdf",width = 8,height = 7.5)



### 4.2 R2  #####

tmp = ddply(model_rmse_r2_by_prots,.(cell_line),summarise,median_r2 =median(r2) )
tmp = tmp[order(tmp$median_r2),]
model_rmse_r2_by_prots$cell_line = factor(model_rmse_r2_by_prots$cell_line,levels = tmp$cell_line)

ggplot(model_rmse_r2_by_prots,aes(x=cell_line,y=r2)) +
	geom_boxplot()  + coord_flip() + theme_bw() +
	ggtitle("Correlation analysis",subtitle = "correlation between model trajectories and data per protein") +
	xlab("Cell lines") + ylab("correlation coefficient")

if(FALSE) ggsave("./supp_info/results_for_march_2019/figures/R2_modeling_boxplot.pdf",width = 8,height = 7.5)


### 4.3  R2 (data_SD > noise) ####

# noise in the data --> null distribution
plot(hist(data_rmse_r2_by_prots$rmse))
noise_095 = quantile(data_rmse_r2_by_prots$rmse,0.95)

tmp = ddply(filter(model_rmse_r2_by_prots,data_sd > noise_095),.(cell_line),summarise,median_r2 =median(r2) )
tmp = tmp[order(tmp$median_r2),]
model_rmse_r2_by_prots$cell_line = factor(model_rmse_r2_by_prots$cell_line,levels = tmp$cell_line)

ggplot(filter(model_rmse_r2_by_prots,data_sd > noise_095),aes(x=cell_line,y=r2)) +
	geom_boxplot()  + coord_flip() + theme_bw() +
	ggtitle("Correlation analysis of changing proteins only",subtitle = "correlation between model trajectories and data per protein") +
	xlab("Cell lines") + ylab("correlation coefficient")

if(FALSE) ggsave("./supp_info/results_for_march_2019/figures/R2_modeling_boxplot_filtered_noise.pdf",width = 8,height = 7.5)
