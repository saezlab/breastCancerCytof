# post fitting analysis

library(multiCellNOpt)
library(plyr)
library(dplyr)
calibrated_model_files = list.files("./data/models/pkn_v4_midas_v3_bootstrap/outputs/","*.RDS",full.names = T)
calibrated_model_files = list.files("./data/models/pkn_v4_midas_v2/outputs/","*.RDS",full.names = T)

calibrated_models = lapply(calibrated_model_files,function(f){
	m = readRDS(f)
	# update the class definition to the newest version
	#m2 = logicODEModel$new(oldModel = m)
}
)

names(calibrated_models) = gsub(".RDS","",basename(calibrated_model_files))

## check first the convergence
# was the length of optimisation enough?


# remove 2 cases which were not finished on the cluster
calibrated_models = calibrated_models[!names(calibrated_models) %in% c("SKBR3_bootstrap_1", "T47D_bootstrap_1","ZR7530_bootstrap_1")]



optimalCosts = as.data.frame(lapply(calibrated_models,function(M){
	if(!is.null(M$fitResults)){
	unlist(lapply(M$fitResults,function(fr){fr$fbest}))}
}))


optimalCosts = melt(optimalCosts,value.name = "costFunction",variable.name = "cell_line")
ggplot(optimalCosts,aes(cell_line,costFunction)) + geom_violin(aes(group=cell_line)) + geom_point() + theme_bw()+coord_flip() + ggtitle('cost distribution of 6 restarts, 40 mins long optimisations')

# the 40 mins optimisation was enough for some but has large variance for others.
# I will increase the number of optimisations from 6 to 10



### Compute R2
# M = calibrated_models[[1]]
get_fitting_stats = function(M){

	SimData = M$simulate()$signals
	ExpData = M$exps$signals
	Nexps = nrow(M$exps$signals[[1]])

	names(SimData) = M$exps$timepoints
	names(ExpData) = M$exps$timepoints

	SimData_df = ldply(SimData,function(D){
		D = as.data.frame(D)
		D$exp = 1:Nexps
		return(D)
	},.id = "time")


	ExpData_df = ldply(ExpData,function(D){
		D = as.data.frame(D)
		D$exp = 1:Nexps
		return(D)
	},.id = "time")


	SimData_df_m = melt(SimData_df,id.vars = c("time","exp"),variable.name = "protein",value.name = "sim_value")
	ExpData_df_m = melt(ExpData_df,id.vars = c("time","exp"),variable.name = "protein",value.name = "exp_value")

	data = merge(ExpData_df_m,SimData_df_m)
	#data = filter(data,exp!=1)  # exp1 is artificial experimental data with no variance --> enforce steady state

	# calculate R2 for each protein and for each experiment
	stats_df = ddply(data,.(protein, exp),function(df){
		# df = filter(data,protein=="STAT5",exp==2)
		# total sum of squares
		SStot = sum((df$exp_value - mean(df$exp_value))^2)
		DataVariance = var(df$exp_value)
		SSres = sum((df$exp_value - df$sim_value)^2)
		SSexpl = SStot - SSres
		RMSE = sqrt(1/length(df$exp_value)*SSres)
		R2 = SSexpl/SStot
		data.frame(SStot=SStot,SSres=SSres,SSexpl=SSexpl,R2=R2,RMSE=RMSE,DataVariance=DataVariance)

	})


}

model_fit_stats = ldply(calibrated_models,get_fitting_stats,.id = "cell_line")

RMSE_sd_categories <- data.frame(
	x = rep(c(0.05,0.15,0.25), 3),
	y = rep(c(0.05,0.15,0.25), each = 3),
	z = factor(c(1,0,1,0,0,0,1,0,1)),
)
####

ggplot(model_fit_stats, aes(RMSE,sqrt(DataVariance)))+geom_tile(data=RMSE_sd_categories,aes(x=x,y=y,fill=z)) + geom_point(aes(col=as.factor(cell_line))) + geom_abline(intercept=0,slope=1)+
	facet_wrap(~protein) + theme_bw() + ggtitle("RMSE vs data stdev")+xlab("RMSE") + ylab("sd(data)") + scale_fill_manual(values = c("1"="grey","0"="white"))






ggplot(stats_df, aes(protein,R2,col=as.factor(exp))) + geom_jitter() + facet_wrap(~protein,scales = "free_x") + theme_bw() + ggtitle("R2 for individual proteins/experiments") + ylab("R2")


M$getStatistics()

M$plotFit()

ggplot(stats_df, aes(protein,R2,col=as.factor(exp))) + geom_jitter() + facet_wrap(~protein,scales = "free_x") + theme_bw() + ggtitle("R2 for individual proteins/experiments") + ylab("R2")
ggplot(stats_df, aes(protein,R2,col=as.factor(exp))) + geom_jitter() + facet_wrap(~protein,scales = "free_x") + theme_bw() + ggtitle("R2 for individual proteins/experiments") + ylab("R2")+ coord_cartesian(ylim = c(0,1))


ggplot(stats_df, aes(RMSE,sqrt(DataVariance),col=as.factor(exp))) + geom_point() +
	facet_wrap(~protein) + theme_bw() + ggtitle("RMSE vs dataVariance")+xlab("RMSE") + ylab("var(data)")


ggplot(stats_df, aes(RMSE,R2,col=as.factor(exp))) + geom_point() +
	facet_wrap(~protein) + theme_bw() + ggtitle("RMSE vs R2")+ylab("R2") + xlab("RMSE")


#### Parameter variance --------------------------------------------------------
# this revealed that the parameter from Serum to nodes and EGF to EGFR are really small,
# so the dynamics were driven mostly by the non-steady state initial conditions.




optimalPars = ldply(calibrated_models,function(M){
	# M = calibrated_models[[2]]
	par_matrix = ldply(M$fitResults,function(fr){fr$xbest},.id = NULL)
	colnames(par_matrix) <- c(M$ode_parameters$parNames[M$ode_parameters$index_opt_pars],M$ode_parameters$x0Names)
	par_matrix
},.id="cell_line_name")



optimalPars = melt(optimalPars,id.vars =  "cell_line_name",value.name = "parameter_value",variable.name = "parameter")


ggplot(optimalPars,aes(parameter,parameter_value,col=cell_line_name)) + geom_point()

optimalPars$type = "IC"
optimalPars$type[grep("_k_",optimalPars$parameter)] = "k"
optimalPars$type[grep("tau_",optimalPars$parameter)] = "tau"


ggplot(filter(optimalPars,type=="k",grepl("SERUM|EGF",parameter) ),aes(parameter,parameter_value,col=cell_line_name)) + geom_point() + coord_flip()


# This shows that the inputs are not causeing the dynamics, but the initial conditions... :(


# check this by creating an experiment with no input and see the dynamics:



M = calibrated_models[[1]]
M$plotFit()
M$exps$cues[1,] = 0
M$exps$stimuli[1,] = 0
M$plotFit()



