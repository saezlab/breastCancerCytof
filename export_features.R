# export features

library(multiCellNOpt)
library(plyr)
library(dplyr)

# calibrated_model_files = list.files("./data/models/pkn_v2_midas_v2/outputs/","*.RDS",full.names = T)
calibrated_model_files = list.files("./data/models/pkn_v4_midas_v4/outputs/","*.RDS",full.names = T)
feature_folder = "./data/models/pkn_v4_midas_v4/features/"

dir.create(feature_folder,recursive = T)

calibrated_models = lapply(calibrated_model_files,function(f){
	m = readRDS(f)
	# update the class definition to the newest version
	#m2 = logicODEModel$new(oldModel = m)
}
)

names(calibrated_models) = gsub(".RDS","",basename(calibrated_model_files))


### Model features 1: parameters -----------------------------------------------
# we export the estimated model parameters (k and tau) for each cell-lines.

feature_table_1 = data.frame()
feature_table_1 = ldply(calibrated_models,function(M){
	M$ode_parameters$parValues[M$ode_parameters$index_opt_pars]
},.id = "cell_line")

save(feature_table_1,file=paste0(feature_folder,"/feature_table_1_raw_parameters.RData"))



### Model features 1.5: parameter variance -----------------------------------------------
# we export the estimated model parameters variance a coeff of variation for each cell-lines.


sampling_df = readRDS("./data/models/pkn_v4_midas_v4/mfuSampling/00.estim_params_collected_v2.RDS")
sampling_df_qc = sampling_df[sampling_df$rel_fobj<1.1,]

parameter_columns = 2:140

feature_table_par_sd = ddply(sampling_df_qc,.(cell_line),function(M){
	par_sample = M[,parameter_columns]
	apply(par_sample,2,sd)
},.id = "cell_line")

colnames(feature_table_par_sd)[parameter_columns] = paste0("sd_",colnames(feature_table_par_sd)[parameter_columns])

save(feature_table_par_sd,file=paste0(feature_folder,"feature_table_par_sd.RData"))

feature_table_par_CoV = ddply(sampling_df_qc,.(cell_line),function(M){
	par_sample = M[,parameter_columns]
	apply(par_sample,2,function(x)sd(x)/mean(x))
},.id = "cell_line")
colnames(feature_table_par_CoV)[parameter_columns] = paste0("CoV_",colnames(feature_table_par_CoV)[parameter_columns])

save(feature_table_par_CoV,file=paste0(feature_folder,"feature_table_par_CoV.RData"))





### Model features 2: edge activity --------------------------------------------
# compute for each edge the strengh of the interaction by computing the value of the
# tranfer function (uses the edge parameters and input nodes) over time, then
# takes the mean activity over time for each experiment. Finally we take the maximum
# among each experiment.
source("../KinaseLibraryScan/R/getReactionActivity.R")
# model = calibrated_models[[1]]
# stat = "meanActivity"
# tau_threshold = 1e-3
transfer_functions_to_feature = function(model,stat=c("meanActivity","AUC"), tau_threshold = 1e-3){

	stat = match.arg(stat)
	reacData = getReactionActivity_logicModel(self = model)

	#k2color <- colorRamp(RColorBrewer::brewer.pal(7,"Greys"))

	# reacData$summary contains the activity of the reaction averaged over time.
	# average computed by area under curve and mean value over time
	# computes the median over all experiments and timepoint:
	max_stat = ddply(reacData$summary,.(reacID,modelName),function(df){
		#browser()
		max_stat = median(df[,stat])
		})

	pre_features = max_stat$V1
	names(pre_features) = as.character(max_stat$reacID)

	tau_pars = model$ode_parameters$parValues[model$ode_parameters$index_tau]
	reacs = as.character(max_stat$reacID)

	get_downstream_node_name = function(reacID){
		strsplit(reacID,split = "=")[[1]][[2]]
	}
	#i=1
	features = pre_features
	for(i in seq_along(reacs)){
		outnode = get_downstream_node_name(reacs[[i]])
		outnode_tau = tau_pars[paste0("tau_",outnode)]
		if(outnode_tau <= tau_threshold){
			features[reacs[i]] = 0
		}
	}
	return(features)

}

feature_table_2 = ldply(calibrated_models,function(M){
	features = transfer_functions_to_feature(model = M,stat = "meanActivity", tau_threshold = 1e-3)
},.id = "cell_line" )
save(feature_table_2,file=paste0(feature_folder,"/feature_table_2_edge_mean_activity.RData"))


## Feature 2 worked pretty well.


### Model feautures 2.2 : timecourse edge strenght -----------------------------
timecourse_strength_to_feature = function(model){
	reacData = getReactionActivity_logicModel(self = model)
	return(reacData$timecourse)
}

feature_table_2_2 = ldply(calibrated_models,function(M){
	features = timecourse_strength_to_feature(model = M)
},.id = "cell_line" ,.progress = progress_text())


feature_table_2_2_timecourse_edge_strength = dcast(feature_table_2_2,cell_line+exp+time+modelName~reacID,value.var = "y")
# remove the control exp
feature_table_2_2_timecourse_edge_strength = feature_table_2_2_timecourse_edge_strength[feature_table_2_2_timecourse_edge_strength$exp != "exp 1",]

save(feature_table_2_2_timecourse_edge_strength,file=paste0(feature_folder,"/feature_table_2_time_course_activity.RData"))


# import and interpolate for all timepoints
load(file=paste0(feature_folder,"/feature_table_2_time_course_activity.RData"))
# load variable feature_table_2_2_timecourse_edge_strength

head(feature_table_2_2_timecourse_edge_strength)
unique_time = unique(feature_table_2_2_timecourse_edge_strength$time)


feature_cols = grep("=", colnames(feature_table_2_2_timecourse_edge_strength))

feature_table_2_3_full_timecourse_edge_strength = ddply(feature_table_2_2_timecourse_edge_strength,.(cell_line,exp),function(df){
	# df = filter(feature_table_2_3_full_timecourse_edge_strength,cell_line=="184A1",exp=="exp 3")
 	# approx(x=df$time,y = df$`PIP3=AKT_S473`,xout = unique_time)
	interp_matrix = apply(df[,feature_cols],2,function(col_dat){approx(x = df$time,y=col_dat,xout = unique_time,method = "linear",rule = 2)[['y']]})
	cbind(time=unique_time,interp_matrix)
})

# exp 1 — not included
# exp 2 —  EGF+SERUM
# exp 3 — iEGFR
# exp 4 — iMEK
# exp 5 — imTOR
# exp 6 — iPI3K
# exp 7 — iPKC
colnames(feature_table_2_3_full_timecourse_edge_strength)[[2]] = "treatment"
feature_table_2_3_full_timecourse_edge_strength$treatment = as.character(feature_table_2_3_full_timecourse_edge_strength$treatment)
feature_table_2_3_full_timecourse_edge_strength$treatment[feature_table_2_3_full_timecourse_edge_strength$treatment == "exp 2"] = "EGF"
feature_table_2_3_full_timecourse_edge_strength$treatment[feature_table_2_3_full_timecourse_edge_strength$treatment == "exp 3"] = "iEGFR"
feature_table_2_3_full_timecourse_edge_strength$treatment[feature_table_2_3_full_timecourse_edge_strength$treatment == "exp 4"] = "iMEK"
feature_table_2_3_full_timecourse_edge_strength$treatment[feature_table_2_3_full_timecourse_edge_strength$treatment == "exp 5"] = "imTOR"
feature_table_2_3_full_timecourse_edge_strength$treatment[feature_table_2_3_full_timecourse_edge_strength$treatment == "exp 6"] = "iPI3K"
feature_table_2_3_full_timecourse_edge_strength$treatment[feature_table_2_3_full_timecourse_edge_strength$treatment == "exp 7"] = "iPKC"
head(feature_table_2_3_full_timecourse_edge_strength)

saveRDS(feature_table_2_3_full_timecourse_edge_strength,file=paste0(feature_folder,"/feature_table_2_3_full_timecourse_edge_strength.RDS"))


ggplot(filter(feature_table_2_3_full_timecourse_edge_strength,cell_line=="184A1")) + geom_line(aes(time,`ERK12=p90RSK` ,col=treatment))

###  Model features 3 ---------------------------------------------------------
# measured and predicted states at iPI3K, timepoint 13
iexp_iPI3K = 5
predicted_nodes_ipi3k_time13 = ldply(calibrated_models,function(M){
	s = M$simulateStates()
	time_13 = which(s$timepoints==13)
	if(length(time_13)== 0) {
		print(paste("time13_missing ", M$exps$name))
		# then take the closest measured time
		time_13 = which.min(abs(s$timepoints - 13))
	}

	s$signals[[time_13]][iexp_iPI3K,]

},.id = "cell_line")


MIDAS_files = list.files("data/MIDAS_v2/",pattern = ".csv",full.names = T)
cellline_names = gsub(".csv","",basename(MIDAS_files),fixed = T)
MIDAS_files = as.list(MIDAS_files)
names(MIDAS_files) = cellline_names

# i = 1
measured_nodes_ipi3k_time13  =  ldply(MIDAS_files,function(m){
	s = CNOlist$new(m,verbose = F)

	time_13 = which(s$timepoints==13)
	if(length(time_13)== 0) {
		print(paste("time13_missing ", M$exps$name))
		# then take the closest measured time
		time_13 = which.min(abs(s$timepoints - 13))
	}

	s$signals[[time_13]][iexp_iPI3K,]


},.id = "cell_line")


# merge the measured and the predicted states:
# order the same way
measured_nodes_ipi3k_time13 = measured_nodes_ipi3k_time13[order(measured_nodes_ipi3k_time13$cell_line), ]
predicted_nodes_ipi3k_time13 = predicted_nodes_ipi3k_time13[order(predicted_nodes_ipi3k_time13$cell_line),]

# get the nodes from the models that were not measured, for other nodes we use the scaled values:
predicted_nodes_only_iPI3K_time13 = predicted_nodes_ipi3k_time13[,!colnames(predicted_nodes_ipi3k_time13) %in% colnames(measured_nodes_ipi3k_time13)]

all_measured_and_predicted_iPI3K_time13 = cbind(measured_nodes_ipi3k_time13,predicted_nodes_only_iPI3K_time13)

save(all_measured_and_predicted_iPI3K_time13,file="all_measured_and_predicted_iPI3K_time13.RData")



# get all the predicted values plus the measured values for the nodes that are not modelled:
measured_nodes_only_iPI3K_time13 = measured_nodes_ipi3k_time13[,!colnames(measured_nodes_ipi3k_time13) %in% colnames(predicted_nodes_ipi3k_time13)]


measured_and_all_predicted_iPI3K_time13 =  cbind(predicted_nodes_ipi3k_time13,measured_nodes_only_iPI3K_time13)

save(measured_and_all_predicted_iPI3K_time13,file="measured_and_all_predicted_iPI3K_time13.RData")



### Plotting model features on the graph ---------------------------------------
source("./R/plotAnnotatedModel.R")
# my_color_scale = rev(c("#B2182B", RColorBrewer::brewer.pal(7,"RdBu"),"#2166AC"))
my_color_scale = RColorBrewer::brewer.pal(9,"YlGnBu")

num_2_my_color_scale = colorRamp(my_color_scale)
num_2_my_color_scale = colorRamp(c('#AAAAAA','#FF1111'))


reac_colors = rgb(num_2_my_color_scale(feature_table_2[1,2:ncol(feature_table_2)]),maxColorValue = 255)
names(reac_colors) = colnames(feature_table_2)[2:ncol(feature_table_2)]

plotAnnotatedModel(model = calibrated_models[[1]]$pkn$network,
				   CNOlist = calibrated_models[[1]]$exps$convertToS4(),
				   annotEdgeColors = reac_colors,annotEdgeMatching = 'reacID')

pdf("figures/fitted_models_plotFit.pdf",width = 21,height = 3)
l_ply(seq_along(calibrated_models),function(i,M){
	M[[i]]$plotFit(measuredNodesOnly = F)
	title(paste0(i,"-",M[[i]]$exps$name))
},M = calibrated_models,.progress = progress_text())
dev.off()

pdf("figures/fitted_models_plotAnnotatedModel_features_2.pdf",width = 8,height = 8)
l_ply(seq_along(calibrated_models),function(i,M){

	reac_colors = rgb(num_2_my_color_scale(feature_table_2[i,2:ncol(feature_table_2)]),maxColorValue = 255)
	names(reac_colors) = colnames(feature_table_2)[2:ncol(feature_table_2)]

	plotAnnotatedModel(model = calibrated_models[[i]]$pkn$network,
					   CNOlist = calibrated_models[[i]]$exps$convertToS4(),
					   annotEdgeColors = reac_colors,annotEdgeMatching = 'reacID')

	title(paste0(i,"-",gsub(".csv","",M[[i]]$exps$name,fixed = T)))
},M = calibrated_models,.progress = progress_text())
dev.off()


###### plot PKN with gggraph
#
# library(ggraph)
# #> Loading required package: ggplot2
# library(tidygraph)
# library(igraph)
#
# sif = read.table("./data/pkn/cancer_cellLines_v2.sif")
# #sif = unique(sif)
# #write.table(sif,"./data/pkn/cancer_cellLines_v3.sif",quote = F,sep = '\t',col.names = F,row.names = F)
#
# graph = as_tbl_graph(sif[,c(1,3)])
# E(graph)$sign = sif[,2]
# ggraph::create_layout(g = graph,layout = 'graphopt')
# ggraph(graph, layout = 'graphopt') + # 'kk'
# 	geom_edge_fan(aes(start_cap = label_rect(node1.name),
# 					   end_cap = label_rect(node2.name)),
# 				   arrow = arrow(length = unit(4, 'mm'))) +
# 	geom_node_text(aes(label=name))+
# 	theme_graph(foreground = 'steelblue', fg_text_colour = 'white')
#
#
#
# G = igraph::graph_from_data_frame(sif[,c(1,3)])
# igraph::plot.igraph(G)



# computes the AUC of the combined transfer functions described in the reacID over the time.
# then maps this values to color scale
transfer_functions_to_edge_colors = function(model,stat=c("meanActivity","AUC")){
	source("R/getReactionActivity.R")
	stat = match.arg(stat)
	reacData = getReactionActivity_logicModel(self = model,timeSignals = seq(0,600,by=60) )

	k2color <- colorRamp(RColorBrewer::brewer.pal(7,"Greys"))

	edge_value = reacData$summary[reacData$summary$exp=="exp 3",stat]
	names(edge_value) = reacData$summary[reacData$summary$exp=="exp 3","reacID"]

	edge_colors = rgb(k2color(edge_value/max(edge_value)),maxColorValue = 255)
	names(edge_colors) = names(edge_value)
	return(list(edge_colors = edge_colors,
				edge_values = edge_value))
}


node_activation_to_node_colors = function(model, type=c("control_to_low","low_to_high")){

	type = match.arg(type)

	S = model$fillCueValuesInSimulation()

	if(type=="control_to_low"){
		# effect = Signal_low_time2 - Signal_control_time2
		effect = S[[2]][2,] - S[[2]][1,]
	}else if(type=="low_to_high"){
		# effect = Signal_high_time2 - Signal_low_time2
		effect = S[[2]][3,] - S[[2]][2,]
	}
	# the effect is \in [-1, 1], therefore
	# effect/2 + 0.5 is \in [0,1] and 0.5 corresponds to the unscaled_effect == 0


	effect2color <- colorRamp(effect2node_colors)

	node_colors = rgb(effect2color(0.5+0.5*effect),maxColorValue = 255)
	names(node_colors) = names(effect)
	return(list(node_colors = node_colors,
				node_values = effect))
}

transfer_functions_to_bitstring = function(model,reacID_values,threshold=0.01, tau_threshold=1e-8){


	# weight the edge_values with the downstream node's tau parameter
	tau_pars = model$ode_parameters$parValues[model$ode_parameters$index_tau]
	reacs = names(reacID_values)

	get_downstream_node_name = function(reacID){
		strsplit(reacID,split = "=")[[1]][[2]]
	}
	#i=1
	for(i in seq_along(reacs)){
		outnode = get_downstream_node_name(reacs[[i]])
		outnode_tau = tau_pars[paste0("tau_",outnode)]
		if(outnode_tau <= tau_threshold){
			reacID_values[reacs[i]] = 0
		}
	}

	bitString = as.logical(reacID_values>threshold)
	names(bitString) = names(reacID_values)
	return(bitString)
}


source("../KinaseLibraryScan/R/plotAnnotatedModel.R")
plotAnnotatedModel()


