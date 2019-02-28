# prepare models to Fit on the cluster
# select 5 cell lines and run them with different regularisation parameters
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

# current settings V4
# model_folder = "./data/models/pkn_v4_midas_v2/inputs/"  # 10.Jan.2019

model_folder = "./data/models/pkn_v4_midas_v2_regpar/inputs/"  # 16.Jan.2019
PKNfile = "./data/pkn/cancer_cellLines_v4.sif"
MIDAS_files = list.files("data/MIDAS_v2/",pattern = ".csv",full.names = T)


dir.create(model_folder,recursive = T)

reg_par_vector = 10^seq(0,-7,length.out = 8)

# select a subset of models:

model_indx = c(18,30,52,58,67)  # sample(seq_along(MIDAS_files),5)

# i = 1
for(j in seq_along(reg_par_vector)){
	for(i in model_indx){

		M= logicODEModel$new(SIFfile = PKNfile ,exps = MIDAS_files[[i]])
		M$name = MIDAS_files[[i]]

		M$preprocessing(cutNONC = T,compression = T,inhibANDExpansion = F)
		# M$plotModel()
		M$objectiveFunction = M$getDefaultLS(SSpenalty_fac = 0,
											 SScontrolPenalty_fac = 0,
											 SSpenalty_for_unobservedStates = FALSE,
											 use_stdev = F,
											 verbose = F,
											 lambda_tau = -1e-5,
											 lambda_k = reg_par_vector[[j]])

		M$initODE_parameters(LB_k = 0, LB_tau = 0, UB_k = 6, UB_tau = 5, default_n = 2,  opt_n = F, opt_x0 = T)
		M$transfer_function = 6

		modelName = gsub(".csv","",basename(MIDAS_files[[i]]),fixed = T)
		modelName = paste0(modelName,"_lambdaK_",j)

		saveRDS(M,paste0(model_folder,modelName,".RDS"))
	}
}
