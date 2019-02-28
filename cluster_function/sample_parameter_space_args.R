# open models passed by the argument and uses mfu sampler to explore the optima
# inputs args:
#	args[[1]]: path to a logicODEmodel .RDS file
#	args[[2]]: id passed to random_seed

approx_time_limit = 1200  # sec

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	stop("At least one argument must be supplied (input file)", call.=FALSE)
}

library(multiCellNOpt)
library(MfUSampler)

model_file = args[[1]]
run_id = args[[2]]

model = readRDS(model_file)



x.init = c(M$ode_parameters$parValues[M$ode_parameters$index_opt_pars],M$ode_parameters$x0Values)
x= x.init
LB = c(M$ode_parameters$LB[M$ode_parameters$index_opt_pars],M$ode_parameters$x0_LB)
UB = c(M$ode_parameters$UB[M$ode_parameters$index_opt_pars],M$ode_parameters$x0_UB)

logpost2 = function(x,M, x0=x.init,  W = diag(1e-6, nrow = length(x.init))){

	#log.prior <- -0.5 * t(x - x0) %*% solve(W) %*% (x - x0)

	LB = c(M$ode_parameters$LB[M$ode_parameters$index_opt_pars],M$ode_parameters$x0_LB)
	UB = c(M$ode_parameters$UB[M$ode_parameters$index_opt_pars],M$ode_parameters$x0_UB)

	log.prior <- 0 #-0.5 * t(x - x0) %*% W %*% (x - x0)

	if(any(x<LB) | any(x>UB)) log.prior = log.prior - 1e6

	### by default the objective function returns the SSE. no SD -> 1
	# lets assume there is a 0.1 standard error only.
	SSE = M$objectiveFunction(x)/(0.1)^2

	LLK = -SSE

	log.post = log.prior + LLK

	return(log.post)
}


### my version with time constraint and verbosity

my_MfU.Sample.Run <- function(x, f,timelimit=Inf, uni.sampler = c("slice", "ars", "arms", "unimet")
							  , ..., control = MfU.Control(length(x)), nsmp = 10, verbose=F) {
	uni.sampler <- match.arg(uni.sampler)

	t <- proc.time()[3]
	x.smp <- array(NA, dim = c(nsmp, length(x)))

	for (n in 1:nsmp) {
		if(verbose) print(paste0("iteration: ",n))
		x <- MfU.Sample(x = x, f = f, uni.sampler = uni.sampler, control = control, ...)
		x.smp[n, ] <- x
		if(proc.time()[3] - t > timelimit) break
	}
	t <- proc.time()[3] - t

	class(x.smp) <- c("MfU", class(x.smp))
	attr(x.smp, "t") <- t
	return (x.smp)
}

#set.seed(22456)
#set.seed(556484)
#set.seed(153467)
set.seed(147*run_id)

nsmp <-300
beta.smp2 <- my_MfU.Sample.Run(x.init, logpost2, timelimit=approx_time_limit, nsmp = nsmp,
							   control = MfU.Control(n = length(x),
							   					  slice.w = rep(0.01,length(x)),
							   					  slice.lower = LB,
							   					  slice.upper = UB)
							   , M = M,verbose=TRUE)

outfile = gsub(".RDS","",basename(model_file),fixed = T)

saveRDS(beta.smp2,paste0("mfusampler_rundid_", run_id,"_",outfile,".RDS"))

