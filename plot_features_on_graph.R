# plot the performance of the features from Marco on a graph


library(CellNOptR)
library(RColorBrewer)
source("./R/plotAnnotatedModel.R")


model =  readSIF("data/pkn/cancer_cellLines_v4.sif")
cno = CNOlist("./data/MIDAS_v4/184A1.csv")

cno_cut = cutCNOlist(cno,model)
model = preprocessing(data=cno_cut,model=model,compression = T,expansion = F)


feature_files = list.files("./data/data_for_graphs/features/","*.rds",full.names = T)


pdf("./figures/combined_features_on_graphs.pdf")
for(N_feature in 1:11){
	feature = readRDS(feature_files[[N_feature]])
	# fix feature names --------------------------------------------------------------
	# node names in the model:
	#print(model$namesSpecies)

	# to annotate nodes and edges on the graph we use the following convention:
	# nodes are called <nodeName>
	# edges are called <nodeName1>_k_<nodeName2>

	# 1. repalce feature names to node and edge names on the graph:
	feature$variable = as.character(feature$variable)
	feature$variable = gsub("p.","",feature$variable,fixed = T)
	feature$variable = gsub("Akt.Ser473.","AKT_S473",feature$variable,fixed = T)
	feature$variable = gsub("AKT.Thr308.","AKT_T308",feature$variable,fixed = T)
	feature$variable = gsub("GSK3b","GSK3B",feature$variable,fixed = T)
	feature$variable = gsub("Ki.67","Ki-67",feature$variable,fixed = T)
	feature$variable = gsub("MKK3.MKK6","MKK36",feature$variable,fixed = T)
	feature$variable = gsub("X.GSK3B.PI3K.PIP3","!GSK3B+PI3K.PIP3",feature$variable,fixed = T)

	feature$variable = gsub("ERK[12]*","ERK12",feature$variable,fixed = F)
	feature$variable = gsub("MEK[12]*","MEK12",feature$variable,fixed = F)
	# should be fixed for all possible features...


	# now add column type: which tells us if it is an edge or node related feature
	feature$type = "node"  # set all to node then change for edges
	feature$type[grep(".",feature$variable,fixed = T)] = "edge"  # find the . between variable names -> edge

	feature$variable = gsub(".","=",feature$variable,fixed = T)  # replace the . with _k_ characters for the edges

	# check if all nodes are corectly included:
	# print nodes, which are neither in the model nor in the data.
	modelNodes = model$namesSpecies
	dataNodes = colnames(cno@signals[[1]])
	#
	problematic_nodes = feature$variable[feature$type=="node" & !feature$variable %in% union(modelNodes,dataNodes)]  # should be empty otherwise check which names to change
	problematic_edges = feature$variable[feature$type=="edge" & !feature$variable %in% model$reacID]  # should be empty otherwise check which names to change

	if(length(problematic_nodes) > 0 | length(problematic_edges)>0) stop(paste("edge or node naming issue with variable:",problematic_nodes,problematic_edges))

	# 2. Scale the performance scores to map to colors -----------------------------
	# we will the colorRamp() function to map numerical values to colors.
	# important notes:
	# colorRamp maps [0,1] interval to colors, so we have to scale the features
	# we use a diverse color scale (performance score has a sign) so we transform 0 to 0.5 and extreme value to 0 and 1.

	feature$scaled_performance = 0.5 +  feature$performance/(2*max(abs(feature$performance)))

	# check the scaling:
	#plot(feature$performance,feature$scaled_performance)

	# 3. convert feature performance to colors -------------------------------------
	# select diverse color scale:
	# RColorBrewer::display.brewer.all()
	# my_color_scale = RColorBrewer::brewer.pal(9,"Spectral")
	my_color_scale = RColorBrewer::brewer.pal(9,"RdYlBu")


	num_2_my_color_scale = colorRamp(my_color_scale)

	feature$colors = rgb(num_2_my_color_scale(feature$scaled_performance),maxColorValue = 255)

	reac_colors = rep("#BBBBBB",length(model$reacID))
	names(reac_colors) = model$reacID

	reac_colors[feature$variable[feature$type=="edge"]] = feature$colors[feature$type=="edge"]

	node_colors = feature$colors[feature$type=="node"]
	names(node_colors) = feature$variable[feature$type=="node"]



	plotAnnotatedModel(model = model,
					   #CNOlist = cno,
					   #graphvizParams = list(edgecolor="grey80"),
					   annotEdgeColors = reac_colors, annotEdgeMatching = 'reacID',annotNodeColors = node_colors
	)
	title(basename(feature_files[[N_feature]]))
}

dev.off()
