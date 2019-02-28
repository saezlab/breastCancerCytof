# breastCancerCytof
modelling breast cancer cytof data with logic ODEs


### README Cell Line Cytof 

pipeline:
	- prepare_data: 			creates MIDAS frmo raw data  
	- create_celllineModel.R: 	create logicODEmodelfrom MIDAS and PKN, saves to RDS files  
	- run optimisation on the RWTH cluster  
	- drug_model_assoc.R: 		check model parameters against drug resistance  
	- export_features.R 		get model parameters, state value to a matrix to processed by Mi Yang  



### Functions and scripts:

	* prepare_data.R: takes the raw data, prepares for CellNOpt and exports to midas files.
	different version of MIDAS files available:
	midas_v2 was not satisfactory: despite the good fit, since we havent included
	the control experiment, the states trajectories were changing without any input signal
	midas_v3: includes a pseudo control experiment (copying the first timepoint across time),
	but we made the mistake to decrease the cue levels to 0.75 for inhibitors.
	midas_v4: stimuli EGF and SERUM is 0.75, but inhibitors are set to 1.

	* analyse_regularisation_parameter.R: scans the regularisation parameter for k
	between 1e-7 and 1 on a logarithmic scale (every magnitude).

	* correlation_ranking_marker_pairs.R: checks the correlation and covariance of measured
	signal values. Ths was to help integrating some of the signals to the PKN. It turned out we have to calculate,
	the correlation for each cell line and each experiment. then aggregate the scores. Correlation was not a good
	measure when the absulute change in the value was small, therefore covariance was tried.

	* create_cellLineModel.R: creates a cell line model from the PKN and MIDAS files,
	 sets up the parameter ranges for the optimisaiton and prepares everything to run a
	 parameter estimation on the cluster.

	* create_cellLineModes_regularisation_scan.R: same as above but includes the
	setup for different regularisation parameters

	* drug_model_assoc.R: find correaltion between model parameters and drug sensitivity.

	* export_features: takes the calibrated models and extract different features from the models,
	these can be just model parameters, directly nonobservable states or interaction strength
