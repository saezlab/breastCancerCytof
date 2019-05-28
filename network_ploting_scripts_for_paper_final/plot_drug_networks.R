### plot Drug prediction in Cytoscape:
# This is the final pipeline -- hopefully
# we push a network annotated by the performance of each drug + all drug to Cytoscape.

library(RCy3)
library(tidyverse)
library(tidygraph)

source("./network_ploting_scripts_for_paper_final//plot_drug_networks_helper_functions.R")

feature_files = list.files("./data/data_for_graphs/features/","*.rds",full.names = T)
feature_files = list.files("./data/data_for_graphs/features_may3/","*.rds",full.names = T)
pkn_used = "data/pkn/cancer_cellLines_v4.sif"
midas_used  = "./data/MIDAS_v4/184A1.csv"


# get info about nodes and interactions modelled/measured:
exp_info = collect_exp_info(sif_file = pkn_used, midas_file = midas_used)
pknNodes = exp_info$pknNodes
dataNodes = exp_info$dataNodes
modelled_node_names = exp_info$modelled_node_names
modelled_reac_names = exp_info$modelled_reac_names


#### Drug sensitivity per drug --------------------------------------------------

edge_featureMatrix = list()
node_featureMatrix = list()


for(N_feature in seq_along(feature_files)){

	featureMatrix = readRDS(feature_files[[N_feature]])
	#drugName = gsub("feature_BH_25fdr_","",gsub("\ _[0-9]+.rds","",basename(feature_files[[N_feature]])))
	drugName = gsub("feature_BH_25fdr_","",gsub(".rds","",basename(feature_files[[N_feature]])))

	# fix feature names
	featureMatrix$variable = fix_feature_variable_names(variable = featureMatrix$variable,
														modelled_node_names,
														modelled_reac_names)

	# make the features unique:
	# MArco: use CV_raw if both CV_row and MEDIAN is available:
	# so we rank MEDIAN as 2, all other stats as 1 and for each variable we take the highest ranked stats:
	featureMatrix <- featureMatrix %>% as_tibble() %>% mutate(stat_num=ifelse(statistic == "Median",2,1))  %>%
		group_by(variable) %>% arrange(stat_num) %>% slice(1)

	# augment features with modelled nodes/edges -- > for graph plotting
	featureMatrix = merge_features_and_model(featureMatrix,
											 modelled_node_names,
											 modelled_reac_names)

	edge_featureMatrix[[N_feature]] = get_drug_edges(featureMatrix, drugName)
	node_featureMatrix[[N_feature]] = get_drug_nodes(featureMatrix, drugName, dataNodes)
}


node_to_export <- node_featureMatrix %>% reduce(rbind)
edge_to_export <- edge_featureMatrix %>% reduce(rbind)


#### Global drug sensitivity  --------------------------------------------------

feature_data = readRDS("data/data_for_graphs/featureAnova_PandN.rds")

featureMatrix <- feature_data %>%ungroup() %>% filter(statistic %in% c("Edge_interpolated", "CV_raw")) %>%
	filter(sign != "FALSE")
drugName = "all_drug"

# fix feature names
featureMatrix$variable = fix_feature_variable_names(variable = featureMatrix$variable,
													modelled_node_names,
													modelled_reac_names)

# naming changes between Marco's table, so let's rename to ours:
featureMatrix <- featureMatrix %>% select(variable,avarage) %>%
	rename(performance = avarage)

# augment features with modelled nodes/edges -- > for graph plotting
featureMatrix = merge_features_and_model(featureMatrix,
										 modelled_node_names,
										 modelled_reac_names)

edge_to_export_global = get_drug_edges(featureMatrix, drugName)
node_to_export_global = get_drug_nodes(featureMatrix, drugName, dataNodes)




#### Mix Global drug sensitivity  and individual drugs ------------------------


node_performance <- node_to_export %>%
	bind_rows(node_to_export_global) %>% # filter(drug=="all_drug", variable=="and1")
	select(drug,performance,variable,nodetype) %>%
	mutate(drug= gsub("-","_",drug)) %>%
	spread(drug,performance,fill = 0) %>%
	rename(name=variable)


edge_performance <- edge_to_export %>%
	bind_rows(edge_to_export_global) %>%
	mutate(parameter = paste0(from,"_",to))%>%
	select(drug,performance,from,to,parameter) %>%
	mutate(drug= gsub("-","_",drug)) %>%
	spread(drug,performance, fill=0)


# renaming nodes: shorten the node names:
target_names <- c("SMAD23", "MAPKAPK2", "MEK12", "ERK12", "MEK12_S221", "p70S6K", "MSK12")
names(target_names) <-  c("SMAD", "MK2", "MEK", "ERK", "MEK_S221", "S6K","MSK")

# i=1
for(i in seq_along(target_names)){
	node_performance$name[node_performance$name==target_names[[i]]] = names(target_names)[[i]]
	edge_performance$from = gsub(target_names[[i]],names(target_names)[[i]], edge_performance$from,fixed = T)
	edge_performance$to = gsub(target_names[[i]],names(target_names)[[i]], edge_performance$to,fixed = T)
	edge_performance$parameter = gsub(target_names[[i]],names(target_names)[[i]], edge_performance$parameter,fixed = T)
}



G  <-  tbl_graph(nodes = data.frame(node_performance),edges = data.frame(edge_performance), directed = TRUE)

createNetworkFromIgraph(G,title = "DRUG_FEATURES",collection = "DRUG_features")
setVisualStyle("CC_v_drug")
layoutNetwork(layout.name = "apply preferred")
#layoutNetwork(layout.name = "hierarchical")

#cols = RColorBrewer::brewer.pal(2,"Reds")
cols = gplots::col2hex(c("navy", "grey80","firebrick3"))
# light a bit the navyblue:
cols[[1]]= "#0000DD"
#cols = c("#0099FF", "#EEEEEE","#FF4040")
#edge_cols = RColorBrewer::brewer.pal(3,"Set1")
drugs = colnames(node_performance)[3:ncol(node_performance)]

# performance distribution
all_performance = as.numeric(c(edge_to_export$performance, node_to_export$performance,edge_to_export_global$performance, node_to_export_global$performance))

all_drug_performance = as.numeric(c(node_performance$all_drug,edge_performance$all_drug))


# performance mapping :
performance_breaks = c(min(all_performance,na.rm = T),
					   0.9*min(all_drug_performance,na.rm = T),
					   0,
					   0.9*max(all_drug_performance,na.rm = T),
					   max(all_performance,na.rm = T))
performance_colors = c(cols[[1]],cols,cols[[3]])
performance_edgeWidth = c(30,8,30)

# For each cell-line edge strength:
# i_drug = 1
for(i_drug in 1:length(drugs)){
	drug_name = drugs[[i_drug]]
	print(drug_name)

	var_name = paste0(drug_name)
	if(var_name=="all_drug"){
		performance_breaks = c(-max(abs(all_drug_performance),na.rm = T),
							   -0.9*max(abs(all_drug_performance),na.rm = T),
							   0,
							   0.9*max(abs(all_drug_performance),na.rm = T),
							   max(abs(all_drug_performance),na.rm = T))
		performance_breaks_line = performance_breaks[c(1,3,5)]
	}else{
		performance_breaks = c(-max(abs(all_performance),na.rm = T),
							   #0.9*min(all_performance,na.rm = T),
							   0,
							   #0.9*max(all_performance,na.rm = T),
							   max(abs(all_performance),na.rm = T))
		performance_breaks_line = performance_breaks[c(1,2,3)]
	}
	setEdgeColorMapping(var_name, performance_breaks, performance_colors, style.name = "CC_v_drug"); Sys.sleep(0.3)
	setNodeColorMapping(var_name, performance_breaks, performance_colors, style.name = "CC_v_drug"); Sys.sleep(0.3)

	setEdgeLineWidthMapping(var_name,performance_breaks_line,widths = performance_edgeWidth,style.name = "CC_v_drug"); Sys.sleep(0.3)
	#setNodeBorderColorMapping("node_type",table.column.values = c("non-measured","measured","input"),colors = edge_cols,style.name = "CC_v1",mapping.type = "d"); Sys.sleep(0.1)

	layoutNetwork(layout.name = "apply preferred"); Sys.sleep(1)

	exportImage(filename = file.path(getwd(),"supp_info/results_for_april_2019/drug_figs_v6",paste0(var_name,".pdf")),type = "PDF")
}




