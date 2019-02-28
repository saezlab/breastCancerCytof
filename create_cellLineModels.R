# prepare models to Fit on the cluster
# imports each cell lines MIDAS and PKN.
# setup the models, export to RDS files
# RDS files are then moved to the cluser and fitted individually.


library(multiCellNOpt)


# archived settings:  -----

# V2 ------
#model_folder = "./data/models/pkn_v2_midas_v2/inputs/"
#PKNfile = "./data/pkn/cancer_cellLines_v2.sif"

# V3 ------
# model_folder = "./data/models/pkn_v3_midas_v2/inputs/"  # 10.Jan.2019
# PKNfile = "./data/pkn/cancer_cellLines_v3.sif"
# MIDAS_files = list.files("data/MIDAS_v2/",pattern = ".csv",full.names = T)
#

# PKN V4, MIDAS V2
# model_folder = "./data/models/pkn_v4_midas_v2/inputs/"  # 10.Jan.2019
# PKNfile = "./data/pkn/cancer_cellLines_v4.sif"
# MIDAS_files = list.files("data/MIDAS_v2/",pattern = ".csv",full.names = T)

# pkn V4, MIDAS V3 -- with bootstrap -- trial with 1
# model_folder = "./data/models/pkn_v4_midas_v3_bootstrap/inputs/"  # 17.Jan.2019
# PKNfile = "./data/pkn/cancer_cellLines_v4.sif"
# MIDAS_files = list.files("data/MIDAS_v3/",pattern = ".csv",full.names = T)
# re_init_model = "./data/models/pkn_v4_midas_v2/outputs/"  # use these models to initialise parameters
# N_bootstrap = 1

# pkn V4, MIDAS V3 -- with bootstrap
# model_folder = "./data/models/pkn_v4_midas_v3_bootstrap_v2/inputs/"  # 17.Jan.2019
# PKNfile = "./data/pkn/cancer_cellLines_v4.sif"
# MIDAS_files = list.files("data/MIDAS_v3/",pattern = ".csv",full.names = T)
# re_init_model = "./data/models/pkn_v4_midas_v2/outputs/"  # use these models to initialise parameters
# N_bootstrap = 50
#
#
# # pkn V4, MIDAS V4 -- no bootstrap
# model_folder = "./data/models/pkn_v4_midas_v3/inputs/"  # 18.Jan.2019
# PKNfile = "./data/pkn/cancer_cellLines_v4.sif"
# MIDAS_files = list.files("data/MIDAS_v3/",pattern = ".csv",full.names = T)
# re_init_model = "./data/models/pkn_v4_midas_v3_bootstrap/outputs/"  # use these models to initialise parameters
# do.bootstrap = FALSE
# N_bootstrap = 1

# pkn V4, MIDAS V4 -- no bootstrap
model_folder = "./data/models/pkn_v4_midas_v4/inputs/"  # 21.Jan.2019
PKNfile = "./data/pkn/cancer_cellLines_v4.sif"
MIDAS_files = list.files("data/MIDAS_v4/",pattern = ".csv",full.names = T)
re_init_model = "./data/models/pkn_v4_midas_v3_bootstrap/outputs/"  # use these models to initialise parameters
do.bootstrap = FALSE
N_bootstrap = 1


if(do.bootstrap) boot_seed_vec = sample(1:1e5,N_bootstrap) else boot_seed_vec  =1

dir.create(model_folder,recursive = T)



calibrated_model_files = list.files(re_init_model,"*.RDS",full.names = T)

calibrated_models = lapply(calibrated_model_files,function(f){
	m = readRDS(f)
}
)
names(calibrated_models) = gsub("_bootstrap_.*","",  gsub(".RDS","",basename(calibrated_model_files)))

# i = 1
for(i in seq_along(MIDAS_files)){

	M = logicODEModel$new(SIFfile = PKNfile ,exps = MIDAS_files[[i]])
	M$name = gsub(".csv","",basename(MIDAS_files[[i]]))

	M$preprocessing(cutNONC = T,compression = T,inhibANDExpansion = F)

	# M$plotModel()
	# j_bootstrap = 1
	for(j_bootstrap in  1:N_bootstrap){
		M$objectiveFunction = M$getDefaultLS(SSpenalty_fac = 0,
											 SScontrolPenalty_fac = 10,  #
											 SSpenalty_for_unobservedStates = TRUE,
											 use_stdev = F,
											 verbose = F,
											 lambda_tau = -1e-5,
											 lambda_k = 5e-5, # set by the reg par study
											 bootstrap = do.bootstrap,
											 boot_seed = boot_seed_vec[[j_bootstrap]])

		M$initODE_parameters(LB_k = 0, LB_tau = 0, UB_k = 6, UB_tau = 5, default_n =2,  opt_n = F, opt_x0 = T)
		M$transfer_function = 6


		# re-init the models with previously estimated parameters
		known_pars = calibrated_models[[M$name]]$ode_parameters$parValues

		stopifnot(all(names(known_pars) == names(M$ode_parameters$parValues)))


		M$ode_parameters$parValues = known_pars
		M$ode_parameters$x0Values = calibrated_models[[M$name]]$ode_parameters$x0Values
		M$initialConditions = calibrated_models[[M$name]]$initialConditions

		modelName = gsub(".csv","",basename(MIDAS_files[[i]]),fixed = T)

		#modelName = paste(modelName,"bootstrap",j_bootstrap,sep = "_")

		saveRDS(M,paste0(model_folder,modelName,".RDS"))
	}
}
save.image("./data/models/pkn_v4_midas_v4/image_create_cellLineModels.RData")
