getReactionActivity = function( self = NULL, timeSignals = NULL, trfType=as.character(self$transfer_function),ode_parameters=self$ode_parameters){

	require(plyr)
	require(pracma)

	if(!is.null(timeSignals)){
		# merge the timeSignals with the time from the data.
		# we need coinciding simulation and data time.
		tsim = sort(unique(c(timeSignals,unlist(lapply(self$exps,function(x)x$timepoints)))))
	}else tsim = self$exps$timepoints

	trfParams = self$ode_parameters$get_reactionParameters()

	# we plot the transfer functions evaluated at the simulaiton points.

	# change the transfer function for the simulation
	tmp_trf = self$transfer_function
	self$transfer_function = as.double(trfType)
	sim = self$simulateStates()
	self$transfer_function = tmp_trf

	# add cues values to the simulations
	sim = self$fillCueValuesInSimulation(sim=sim)

	# compute the edge activity for each models, because they have different edge-parameters
	trf_timecourse = list()

	for(imodel in 1:length(sim) ){


		signals = ldply(1:length(sim[[imodel]]$signals),function(i, S,T){
			#	browser()
			df = cbind(as.data.frame(S[[i]]),time=T[[i]], exp=paste("exp",1:nrow(S[[i]])))

			return(df)
		},S=sim[[imodel]]$signals,T=sim[[imodel]]$timepoints)


		signals$modelName = sim[[imodel]]$name


		trf_timecourse[[imodel]] = ldply(trfParams[imodel,],.fun = function(parset,signals,trfType){
			#parset = trfParams[[1,1]]
			#if(parset$name=="!ERK12_dm2+FAK=RAS") browser()
			y = 1
			for(j in seq_along(parset$inputNodes)){
				x = signals[,parset$inputNodes[j]]

				# for each reaction we select the transfer function type
				if(trfType=="3"){ # all trf3
					trf = getTransferFunction("3")

				}else if(trfType=="4"){# all trf4
					trf = getTransferFunction("4")
				}else if(trfType=="6"){# activation are tf4, inhibition tf3
					if(parset$edgeSign[j]==1)
						trf = getTransferFunction("4")
					else
						trf = getTransferFunction("3")
				}


				if(parset$edgeSign[j]==1)
					y = y*trf(x,parset$n[j],parset$k[j])
				else
					y = y*(1-trf(x,parset$n[j],parset$k[j]))

			}

			res = data.frame(reacID = parset$reacID, exp=signals[,"exp"], time=signals[,"time"], modelName=signals[,"modelName"], y=y)

			#names(res) = parset$name
			return(res)
		},signals,trfType)


	}
	trf_timecourse = do.call("rbind",trf_timecourse)

	yData=list()
	yData$timecourse = trf_timecourse
	yData$summary = ddply(yData$timecourse,c("reacID","exp","modelName"),summarise, meanActivity = mean(y,na.rm = T), AUC = pracma::trapz(time,y) )

	return(yData)

}




getTransferFunction=function(trfType=c("1","2","3","4")){

	trfType<-match.arg(trfType)
	trf = switch(as.character(trfType),
				 "1" = function(x,n,k)x,
				 "2" = function(x,n,k) x^n/(x^n + k^n),
				 "3" = function(x,n,k){
				 	x = ifelse(x<0, 0, x)
				 	x = ifelse(x>1, 1, x)
				 	if (k==0) {
				 		y = ifelse(x==0, 1.0, x^n/(x^n + k^n)*(1+k^n))
				 	}else
				 		y = x^n/(x^n + k^n)*(1+k^n)
				 	return(y)
				 },
				 "4"=  function(x,n,k){if (k==0) return(x*0.0)
				 	x = ifelse(x<0, 0, x)
				 	x = ifelse(x>1, 1, x)
				 	y = 1-(1-x)^n/((1-x)^n + k^n)*(1+k^n)
				 	return(y)},
				 "default"= stop("no such transfer function type")
	)
}
