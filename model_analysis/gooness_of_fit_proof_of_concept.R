# goodness of fit test
# this is a conceptual work if correaltion between measured and simulated values is useful or not



library(multiCellNOpt)
library(plyr)
library(dplyr)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

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

# 2. add R2

r2_par_test = filter(sampling_df_qc,cell_line=="184A1")
r2_par_test$cond = "sampled"
# randomised parameters with 10% uniform perturbation
index_par = 2:140

r2_par_test_rand = r2_par_test

index_k = grep("_k_", colnames(r2_par_test))
index_tau = grep("tau_", colnames(r2_par_test))
#r2_par_test_rand[,index_par] = r2_par_test[,index_par] * matrix(runif(prod(dim(r2_par_test[,index_par])),0.5,1.5), nrow = nrow(r2_par_test[,index_par]),ncol = ncol(r2_par_test[,index_par]))
r2_par_test_rand[,index_k] =  matrix(runif(prod(dim(r2_par_test[,index_k])),0,6), nrow = nrow(r2_par_test[,index_k]),ncol = ncol(r2_par_test[,index_k]))
r2_par_test_rand[,index_tau] =  matrix(runif(prod(dim(r2_par_test[,index_tau])),0,5), nrow = nrow(r2_par_test[,index_tau]),ncol = ncol(r2_par_test[,index_tau]))

r2_par_test_rand$cond = 'perturbed'

r2_par_test = rbind(r2_par_test,r2_par_test_rand)

r2_par_test_stats = adply(r2_par_test,.margins = 1,function(x){
	#x = r2_par_test[1,]

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

	return(cbind(stats,R2=R2, corr2 = corr_coeff))
},.progress = progress_text())


ggplot(r2_par_test_stats) + geom_density(aes(corr2,fill=cond),alpha=.5)
ggplot(r2_par_test_stats) + geom_density(aes(RMSE,fill=cond),alpha=.5)

ggplot(r2_par_test_stats) + geom_point(aes(RMSE,corr2,col=cond),size=.5)


### Cehck correlation based on different groupings
x = r2_par_test[1,]

M = calibrated_models[[x$cell_line]]
#M$plotFit()
M$ode_parameters$parValues[M$ode_parameters$index_opt_pars] = as.numeric(x[1,M$ode_parameters$parNames[M$ode_parameters$index_opt_pars]])
M$ode_parameters$x0Values[M$ode_parameters$index_opt_x0] = as.numeric(x[1,gsub("^1","x0",M$ode_parameters$x0Names[M$ode_parameters$index_opt_x0])])
M$updateInitialValueEstimates(ode_parameters = M$ode_parameters,initialConditions = M$initialConditions)

data_vec = unlist(M$exps$signals)
Sim = M$simulate()$signals
sim_vec = unlist(Sim)


SStot = sum((data_vec-mean(data_vec))^2)
SSres= sum((data_vec-sim_vec)^2)

R2 = 1  - SSres/SStot
print(R2)
corr_coeff = cor(sim_vec,data_vec)^2
print(corr_coeff)

# plot(data_vec,sim_vec)
exps = rep(paste("exp",1:7,sep=""), ncol(Sim[[1]])*length(Sim))
times = rep(paste("time",1:length(Sim),sep=""), each=ncol(Sim[[1]])*nrow(Sim[[1]]))
prots = rep(rep(colnames(Sim[[1]]), each=nrow(Sim[[1]])),length(Sim))
ggdata = data.frame(sim=sim_vec,dat=data_vec,exps,times,prots)

r2_by_time = ddply(ggdata,.(times),summarize,r2 = cor(sim,dat)^2)
r2_by_prots = ddply(ggdata,.(prots),summarize,r2 = cor(sim,dat)^2)

rmse_by_prots = ddply(ggdata,.(prots),summarize,rmse = sqrt(sum((sim-dat)^2)/length(sim) ))


sqrt(sum((ggdata$sim-ggdata$dat)^2)/length(ggdata$sim))

mean_range_2_text = function(x,digits=2){

	m = round(mean(x),digits)
	r1 = round(range(x)[[1]],digits)
	r2 = round(range(x)[[2]],digits)

	paste(m," [", r1,", ",r2,"]",collapse="")

}

ggplot(r2_by_time) + geom_boxplot(aes(x=1,r2),width=0.3,fill="grey90") +
	#geom_dotplot(aes(x=rep(1,nrow(r2_by_time)),y=r2,fill=times,size=10),binaxis = "y", stackdir = "center") +
	geom_point(aes(x=rep(1,nrow(r2_by_time)),y=r2,col=times),size=3) +
	theme_bw() +ggtitle("r2 by timepoints") + coord_cartesian(ylim=c(0,1),xlim = c(0.5,1.5)) + geom_text(aes(x=1.10,y=0.8,label=mean_range_2_text(r2_by_time$r2)   ))

ggplot(r2_by_prots) + geom_boxplot(aes(x=1,r2),width=0.3,fill="grey90") +
	#geom_dotplot(aes(x=rep(1,nrow(r2_by_time)),y=r2,fill=prots,size=10),binaxis = "y", stackdir = "center") +
	geom_point(aes(x=rep(1,nrow(r2_by_prots)),y=r2,col=prots),size=3) +
	theme_bw() +ggtitle("r2 by proteins") + coord_cartesian(ylim=c(0,1),xlim = c(0.5,1.5)) +
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
