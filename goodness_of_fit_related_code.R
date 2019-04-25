# goodness of fit for cell line models



library(multiCellNOpt)
library(plyr)
library(dplyr)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
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



# 2. add R2 and other Stats ----------------------------------------------------

r2_par_test = best_parameters

index_par = 2:140


r2_par_test_stats = adply(r2_par_test,.margins = 1,function(x){
	 # x = r2_par_test[1,]
		#x = filter(best_parameters,cell_line =="DU4475")

	M = calibrated_models[[x$cell_line]]
	#M$plotFit()
	M$ode_parameters$parValues[M$ode_parameters$index_opt_pars] = as.numeric(x[1,M$ode_parameters$parNames[M$ode_parameters$index_opt_pars]])
	M$ode_parameters$x0Values[M$ode_parameters$index_opt_x0] = as.numeric(x[1,gsub("^1","x0",M$ode_parameters$x0Names[M$ode_parameters$index_opt_x0])])
	M$updateInitialValueEstimates(ode_parameters = M$ode_parameters,initialConditions = M$initialConditions)

	stats = M$getStatistics()

	data_vec = unlist(M$exps$signals)
	Sim = M$simulate()$signals
	sim_vec = unlist(Sim)
	SStot = sum((data_vec-mean(data_vec))^2)
	SSres= sum((data_vec-sim_vec)^2)

	R2 = 1  - SSres/SStot
	# print(R2)
	corr_coeff = cor(sim_vec,data_vec)^2
	# print(corr_coeff)


	return(cbind(stats,R2=R2, corr2 = corr_coeff ))
},.progress = progress_text())



### Global stats plots ----------------------------------------------------------
# these are influenced by Simson's paradox

# Overall corr2 density plot ( care Simson's paradox is playing big role)
ggplot(r2_par_test_stats) + geom_density(aes(corr2),fill="firebrick") +
	xlab("Correlation coeff.") +
	ggtitle("Correlation coefficients for all celllines")+
	theme_bw() + coord_cartesian(xlim=c(0,1))

# Overall corr2 box plot ( care Simson's paradox is playing big role)
ggplot(r2_par_test_stats) + geom_boxplot(aes(factor("cell_lines"),corr2)) +
	geom_point(aes(1,corr2)) +
	ylab("Correlation coeff.") + xlab("") +
	ggtitle("Correlation coefficients for all celllines")+
	theme_bw() + coord_cartesian(ylim=c(0,1))

# Overall RMSE box plot
ggplot(r2_par_test_stats) + geom_boxplot(aes(factor("cell_lines"),RMSE)) +
	geom_point(aes(1,RMSE)) +
	ylab("RMSE.") + xlab("") +
	ggtitle("RMSE for all celllines")+
	theme_bw() + coord_cartesian(ylim=c(0,0.5))


ggplot(r2_par_test_stats,aes(RMSE,corr2)) + geom_point() +
	geom_text_repel(aes(label=cell_line)) +
	theme_bw() + coord_cartesian(xlim = c(0,0.12),ylim = c(0.5,1))


#######  Compute values for each cell line


# get_annot_df_from_cnolist - - helper funciton
# arrange the CNOlist$signals into a long data.frame format using annotation

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


rmse_r2_by_prots = ddply(sim_vs_data_all,.(cell_line,prots),summarize,
						 rmse = sqrt(sum((sim-dat)^2)/length(sim)),
						 r2 = cor(sim,dat)^2,
						 data_sd = sd(dat))




mean_range_2_text = function(x,digits=2){

	m = round(mean(x),digits)
	r1 = round(range(x)[[1]],digits)
	r2 = round(range(x)[[2]],digits)

	paste(m," [", r1,", ",r2,"]",collapse="")

}

rmse_r2_by_prots$cell_line = as.character(rmse_r2_by_prots$cell_line)

ranked_df = ddply(rmse_r2_by_prots,.(cell_line),function(df)r2=median(df$r2))
ranked_df = ddply(rmse_r2_by_prots,.(cell_line),function(df)r2=median(df$data_sd))
rmse_r2_by_prots$cell_line = factor(rmse_r2_by_prots$cell_line,levels =ranked_df$cell_line[order(ranked_df$V1,decreasing = T)])




ggplot(rmse_r2_by_prots) +
	#geom_violin(aes(x=cell_line,r2),color ="black", fill="firebrick",scale = "width") +
	geom_boxplot(aes(x=cell_line,r2,weight=(data_sd-min(data_sd))/(max(data_sd)-min(data_sd))),color ="black", fill="firebrick",varwidth = T) +
	#geom_boxplot(aes(x=cell_line,data_sd),width=0.1,fill="white",col="black") +
	#geom_dotplot(aes(x=rep(1,nrow(r2_by_time)),y=r2,fill=times,size=10),binaxis = "y", stackdir = "center") +
	#geom_point(aes(x=cell_line,y=r2),size=3) +
	theme_bw() +ggtitle("r2 by cell lines over proteins") + coord_flip(ylim=c(0,1))


## remove the r2 for those where the data_sd is too small (no changein over time)

rmse_r2_by_prots_rm = rmse_r2_by_prots
rmse_r2_by_prots_rm$r2[rmse_r2_by_prots_rm$data_sd<0.1] = NA

ranked_df = ddply(rmse_r2_by_prots_rm,.(cell_line),function(df)r2=median(df$r2,na.rm = T))
rmse_r2_by_prots_rm$cell_line = factor(rmse_r2_by_prots_rm$cell_line,levels =ranked_df$cell_line[order(ranked_df$V1,decreasing = F)])


ggplot(rmse_r2_by_prots_rm) +
	#geom_violin(aes(x=cell_line,r2),color ="black", fill="firebrick",scale = "width") +
	geom_boxplot(aes(x=cell_line,r2),color ="black", fill="grey90") +
	#geom_boxplot(aes(x=cell_line,data_sd),width=0.1,fill="white",col="black") +
	#geom_dotplot(aes(x=rep(1,nrow(r2_by_time)),y=r2,fill=times,size=10),binaxis = "y", stackdir = "center") +
	#geom_point(aes(x=cell_line,y=r2),size=3) +
	theme_bw() +ggtitle("r2 by cell lines over proteins") + coord_flip(ylim=c(0,1))



ggplot(rmse_r2_by_prots) + geom_boxplot(aes(x=cell_line,y=data_sd))

ggplot(rmse_r2_by_prots) + geom_boxplot(aes(x=cell_line,rmse),width=0.3,fill="grey90") +
	#geom_dotplot(aes(x=rep(1,nrow(r2_by_time)),y=r2,fill=times,size=10),binaxis = "y", stackdir = "center") +
	#geom_point(aes(x=cell_line,y=rmse),size=3) +
	theme_bw() +ggtitle("rmse by proteins")  + coord_flip(ylim=c(0,0.25))


ggplot(rmse_r2_by_prots) + geom_point(aes(x=data_sd,r2)) + facet_wrap(~cell_line)
	#geom_dotplot(aes(x=rep(1,nrow(r2_by_time)),y=r2,fill=times,size=10),binaxis = "y", stackdir = "center") +
	#geom_point(aes(x=cell_line,y=rmse),size=3) +
	theme_bw() +ggtitle("rmse by proteins")  + coord_flip(ylim=c(0,0.25))




ggplot(r2_by_prots) + geom_boxplot(aes(x=1,r2),width=0.3,fill="grey90") +
	#geom_dotplot(aes(x=rep(1,nrow(r2_by_time)),y=r2,fill=prots,size=10),binaxis = "y", stackdir = "center") +
	geom_point(aes(x=rep(1,nrow(r2_by_prots)),y=r2,col=prots),size=3) +
	theme_bw() +ggtitle("r2 by proteins") + coord_cartesian(ylim=c(0,1)) +
	geom_text(aes(x=1.10,y=0.8,label=paste("r2=",mean_range_2_text(r2_by_prots$r2) )),hjust = 0)


ggplot() + geom_boxplot(data=r2_by_prots,aes(x=1,r2),width=0.3,fill="grey90") +
	geom_point(data=r2_by_prots, aes(x=rep(1,nrow(r2_by_prots)),y=r2,col=prots),size=3) +
	geom_text(aes(x=1.10,y=0.8,label=paste("r2=",mean_range_2_text(r2_by_prots$r2) )),hjust = 0,size=6)+
	geom_text_repel(data=r2_by_prots, aes(x=rep(1,nrow(r2_by_prots)),y=r2,label=prots)) +

	geom_boxplot(data=rmse_by_prots,aes(x=2,rmse),width=0.3,fill="grey90") +
	geom_point(data=rmse_by_prots, aes(x=rep(2,nrow(rmse_by_prots)),y=rmse,col=prots),size=3) +
	theme_bw() +ggtitle("r2 and RMSE by proteins") + coord_cartesian(ylim=c(0,1),xlim = c(0.5,2.5)) +
	geom_text(aes(x=1.5,y=0.25,label=paste("RMSE=",mean_range_2_text(rmse_by_prots$rmse) )),hjust = 0,size=6)


ggplot() +
	geom_point( aes(x=rmse_by_prots$rmse,y=r2_by_prots$r2,col=r2_by_prots$prots),size=3)+
	theme_bw() +ggtitle("r2 vs RMSE by proteins") + coord_cartesian(ylim=c(0,1),xlim = c(0,0.2)) +
	geom_text_repel( aes(x=rmse_by_prots$rmse,y=r2_by_prots$r2,label=r2_by_prots$prots))


rmse_r2_dat = merge(rmse_by_prots, r2_by_prots)
rmse_r2_dat = melt(rmse_r2_dat,id.vars = "prots")


ggplot(rmse_r2_dat) + geom_boxplot(aes(x=0,value,fill=variable),width=0.3) +
	geom_point(aes(x=0,y=value),size=3) + expand_limits(y=0,x=c(-0.5,0.5)) +
	geom_text_repel(aes(x=0,y=value,label=prots))  + guides(fill=FALSE)+
	facet_wrap(~variable,scales = "free_y") + theme_bw()



ggplot(rmse_r2_dat) +
	geom_point( aes(x=1,y=value),size=3)+
	theme_bw() +ggtitle("r2 vs RMSE by proteins") + coord_cartesian(ylim=c(0,1),xlim = c(0,0.2)) +
	geom_text_repel( aes(x=rmse_by_prots$rmse,y=r2_by_prots$r2,label=r2_by_prots$prots))


ggplot(ggdata,aes(x=sim,y=dat)) + geom_point()
ggplot(ggdata,aes(x=sim,y=dat)) + geom_point(aes(col=times),size=.5)
ggplot(ggdata,aes(x=sim,y=dat)) + geom_point(aes(col=prots),size=.5)
ggplot(ggdata,aes(x=sim,y=dat)) + geom_point(aes(col=exps),size=.5)

ggplot(ggdata,aes(x=sim,y=dat)) + geom_point(size=.5) + facet_wrap(~prots)


range_data = ddply(ggdata,.(prots),function(x)range(x$dat))
range_data$prots = factor(range_data$prots,levels = range_data$prots[order(range_data$V2-range_data$V1)])
range_data$prots = factor(range_data$prots,levels = range_data$prots[order(range_data$V1)])
ggplot(range_data) + geom_errorbar(aes(x=prots,ymin=V1,ymax=V2)) + coord_flip() + ylab("range of values")
