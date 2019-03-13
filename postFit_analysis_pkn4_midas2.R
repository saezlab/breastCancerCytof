# post fitting analysis

library(multiCellNOpt)
library(plyr)
library(dplyr)
calibrated_model_files = list.files("./data/models/pkn_v4_midas_v2/outputs/","*.RDS",full.names = T)

calibrated_models = lapply(calibrated_model_files,function(f){
	m = readRDS(f)
	# update the class definition to the newest version
	#m2 = logicODEModel$new(oldModel = m)
}
)

names(calibrated_models) = gsub(".RDS","",gsub("early_model_","",basename(calibrated_model_files)))

## check first the convergence
# was the length of optimisation enough?


calibrated_models[[2]]$reportConvergence()

calibrated_models[[2]]$fitResults[[1]]$fbest


optimalCosts = as.data.frame(lapply(calibrated_models,function(M){
	unlist(lapply(M$fitResults,function(fr){fr$fbest}))
}))


optimalCosts = melt(optimalCosts,value.name = "costFunction",variable.name = "cell_line")
ggplot(optimalCosts,aes(cell_line,costFunction)) + geom_violin(aes(group=cell_line)) + geom_point() + theme_bw() + ggtitle('cost distribution of 6 restarts, 40 mins long optimisations')

# the 40 mins optimisation was enough for some but has large variance for others.
# I will increase the number of optimisations from 6 to 10



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



