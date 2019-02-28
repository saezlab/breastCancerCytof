# open models passed by the argument and fits it 


maxTime = 1200
maxEval = 1e5
local_solver = "DHC"
nRun = 6
nCores = 1
ndiverse = 30
dim_refset = 30


args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	stop("At least one argument must be supplied (input file)", call.=FALSE)
} 

library(multiCellNOpt)

model_file = args[[1]]

model = readRDS(model_file)

model$fit(maxEval = maxEval,
		  maxTime = maxTime,
		  local_solver = local_solver,
		  nRun = nRun,
		  nCores = nCores,
		  ndiverse = ndiverse,
		  dim_refset = dim_refset)


saveRDS(model,model_file)

#
