# Goodness of fit again: compare model error to biological reproducibility


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
best_parameters = ddply(sampling_df_qc,.(cell_line),function(df) df[which.min(df$fobj),])

get_annot_df_from_cnolist = function(signal_list,data_list){
	exps = rep(paste("exp",nrow(signal_list[[1]]),sep=""), ncol(signal_list[[1]])*length(signal_list))
	times = rep(paste("time",1:length(signal_list),sep=""), each=ncol(signal_list[[1]])*nrow(signal_list[[1]]))
	prots = rep(rep(colnames(signal_list[[1]]), each=nrow(signal_list[[1]])),length(signal_list))
	ggdata = data.frame(sim=unlist(signal_list),dat=unlist(data_list),exps,times,prots)
}

sim_vs_data_all = adply(best_parameters,.expand = F,.margins = 1,function(x){
	#x = r2_par_test[1,]

	M = calibrated_models[[x$cell_line]]
	#M$plotFit()
	M$ode_parameters$parValues[M$ode_parameters$index_opt_pars] = as.numeric(x[1,M$ode_parameters$parNames[M$ode_parameters$index_opt_pars]])
	M$ode_parameters$x0Values[M$ode_parameters$index_opt_x0] = as.numeric(x[1,gsub("^1","x0",M$ode_parameters$x0Names[M$ode_parameters$index_opt_x0])])
	M$updateInitialValueEstimates(ode_parameters = M$ode_parameters,initialConditions = M$initialConditions)

	ggdata = get_annot_df_from_cnolist(signal_list =  M$simulate()$signals, data_list = M$exps$signals)
	ggdata$cell_line = x$cell_line
	return(ggdata)
},.progress = progress_text())



# sim vs mean data

model_rmse_r2_by_prots = ddply(sim_vs_data_all,.(cell_line,prots),summarize,
						 rmse = sqrt(sum((sim-dat)^2)/length(sim)),
						 r2 = cor(sim,dat)^2,
						 r = cor(sim,dat),
						 data_sd = sd(dat))





######## 3. COMPARE MODEL AND DATA STATISTICS ##################################


model_rmse_r2_by_prots$type = "model_vs_data"
data_rmse_r2_by_prots$type = "A/B_vs_mean"
data_rmse_r2_by_prots$markers = as.character(data_rmse_r2_by_prots$markers)
data_rmse_r2_by_prots$markers = gsub("p-","",data_rmse_r2_by_prots$markers)

colnames(model_rmse_r2_by_prots)[[2]] = "markers"
model_rmse_r2_by_prots$markers = as.character(model_rmse_r2_by_prots$markers)


combined_stats = rbind(data_rmse_r2_by_prots,model_rmse_r2_by_prots[,-6])


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
	xlab("Cell lines") + ylab("Root mean square error (RMSE)")

if(FALSE) ggsave("./supp_info/results_for_march_2019/figures/RMSE_modeling_vs_biological_mean_exp.pdf",width = 8,height = 7.5)


ggplot(combined_stats,aes(x=rmse,col=type,fill=type)) +
	geom_density(alpha=0.5)  +
	 xlab("Root mean square error (RMSE)") + theme_bw()
if(FALSE) ggsave("./supp_info/results_for_march_2019/figures/RMSE_modeling_vs_biological_density.pdf",width = 8,height = 3)


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
