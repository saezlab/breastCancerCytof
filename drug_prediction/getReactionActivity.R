# input: logicODEmodel
# output: list, each item on the list is a list that corresponds to a ReacID. It contains
# 	reacid, input and output nodes, reaction parameters
get_reactionParameters_logicODEModel=function(model){

	reacID = model$pkn$network$reacID

	nReactions = length(reacID)

	trfParams = list()
	odeParameters = model$ode_parameters


	mysplit = function(x){strsplit(x, "=")}
	mysplit2 = function(x){strsplit(x, "+", fixed=TRUE)}
	#trfParams  = vector("list",nTransferFunctions)
	trfName=list()



	reacs = mysplit(reacID) # separate the reactions into left and right parts
	tmp <- unlist(mysplit(reacID))
	reacs = t(matrix(unlist(mysplit(reacID)), ncol=length(tmp)/2)) # reordering

	#  i =1
	for (i in 1:dim(reacs)[1]){
		inputs<-unlist(strsplit(reacs[i,1],"+", fixed=TRUE))

		inputSet = vector()
		outputSet = reacs[i,2]
		edgeSign = vector()
		par_k_set = vector()
		par_n_set = vector()
		for (j in seq_along(inputs)){

			v1<-inputs[j]
			edgeSign[j] = 1
			# detect inhibitory edge
			if (length(grep("!", v1))){
				v1 = sub("!", "", v1)
				edgeSign[j] = -1
			}

			inputSet[j] = v1
			par_k_set[j] = odeParameters$parValues[paste0(v1,"_k_",outputSet)]
			par_n_set[j] = odeParameters$parValues[paste0(v1,"_k_",outputSet)]


		}

		trfParams[[i]] = list(name=reacID[[i]],
							  reacID = reacID[[i]],
							  inputNodes=inputSet,
							  endNode=outputSet,
							  n = par_n_set,
							  k = par_k_set,
							  edgeSign=edgeSign)
	}
	trfParams

}


getReactionActivity_logicModel = function( logicODEmodel = NULL, timeSignals = NULL, trfType=as.character(logicODEmodel$transfer_function),ode_parameters=logicODEmodel$ode_parameters){

	require(plyr)
	require(pracma)

	# update time if needed
	if(!is.null(timeSignals)){
		tsim = sort(unique(timeSignals))
	}else tsim = logicODEmodel$exps$timepoints

	trfParams = get_reactionParameters_logicODEModel(logicODEmodel)

	# we plot the transfer functions evaluated at the simulaiton points.

	# change the transfer function for the simulation
	tmp_trf = logicODEmodel$transfer_function
	logicODEmodel$transfer_function = as.double(trfType)
	sim = logicODEmodel$simulateStates(timeSignals = tsim)
	logicODEmodel$transfer_function = tmp_trf

	# add cues values to the simulations
	sim$signals = logicODEmodel$fillCueValuesInSimulation(signals=sim$signals)

	# compute the edge activity for each models, because they have different edge-parameters
	trf_timecourse = list()

	signals = ldply(1:length(sim$signals),function(i, S,T){
		#	browser()
		df = cbind(as.data.frame(S[[i]]),time=T[[i]], exp=paste("exp",1:nrow(S[[i]])))

		return(df)
	},S=sim$signals,T=sim$timepoints)




	signals$modelName = sim$name


	trf_timecourse = ldply(trfParams,.fun = function(parset,signals,trfType){
		#parset = trfParams[[1,1]]
		#if(parset$name=="!ERK12_dm2+FAK=RAS")

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
