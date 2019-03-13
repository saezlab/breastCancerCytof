# create models

library(multiCellNOpt)


list.files("data/MIDAS/")

cnodata = CNOlist$new("./data/MIDAS_v2//CAL120.csv")
cnodata$plot()

pkn=PKN$new("./data/pkn/cancer_cellLines_v2.sif")
pkn$plotModel()
pkn=PKN$new("./data/pkn/all_marker_early_sharedAND.sif")

pkn$plotModel(cnolist = cnodata$convertToS4())




M= logicODEModel$new(SIFfile = "./data/pkn/cancer_cellLines_v2.sif",exps = "./data/MIDAS_v2/CAL120.csv")

M$preprocessing(cutNONC = T,compression = F,inhibANDExpansion = F)
M$plotModel()
M$plotData()
M$preprocessing(cutNONC = T,compression = T,inhibANDExpansion = F)
M$plotModel()

M$objectiveFunction = M$getDefaultLS(SSpenalty_fac = 10,
													   SScontrolPenalty_fac = 10,
													   SSpenalty_for_unobservedStates = TRUE,
													   use_stdev = F,verbose = F,lambda_tau = -1e-4)

M$initODE_parameters(LB_k = 0, LB_tau = 0, UB_k = 6, UB_tau = 5, opt_n = T, opt_x0 = T)
M$transfer_function = 6

M$fit(maxEval = 1e5,maxTime = 60,local_solver = 'DHC',nRun = 1,nCores = 1,ndiverse = 20,dim_refset = 20)

M$plotFit(measuredNodesOnly = F)
M$initialConditions
M$simulate()
