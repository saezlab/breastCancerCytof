
# Gets the nodes and reacid's used in the CNORode model, PKN and MIDAS files.
collect_exp_info <- function(sif_file, midas_file){

	model =  CellNOptR::readSIF(sif_file)
	cno = CellNOptR::CNOlist(midas_file)

	cno_cut = CellNOptR::cutCNOlist(cno,model)
	model = CellNOptR::preprocessing(data=cno_cut,model=model,compression = T,expansion = F)

	# get the names of nodes/edges that appear in the model or in data:
	pknNodes = model$namesSpecies
	dataNodes = colnames(cno@signals[[1]])
	modelled_node_names = union(pknNodes,dataNodes)
	modelled_reac_names = model$reacID

	return(list(modelled_node_names = modelled_node_names,
				modelled_reac_names = modelled_reac_names,
				dataNodes = dataNodes,
				pknNodes = pknNodes))

}


# this function should be adjusted by hand
# we take the variable names found in Marco's feature matrix and convert to names
# the modelled network
fix_feature_variable_names <- function(variable,modelled_node_names,modelled_reac_names){
	# fix feature names --------------------------------------------------------------
	# node names in the model:
	#print(model$namesSpecies)
	# browser()

	# to annotate nodes and edges on the graph we use the following convention:
	# nodes are called <nodeName>
	# edges are called <nodeName1>_k_<nodeName2>

	# 1. repalce feature names to node and edge names on the graph:
	variable = as.character(variable)
	variable = gsub("p.","",variable,fixed = T)
	variable = gsub("Akt.Ser473.","AKT_S473",variable,fixed = T)
	variable = gsub("AKT.Thr308.","AKT_T308",variable,fixed = T)
	variable = gsub("GSK3b","GSK3B",variable,fixed = T)
	variable = gsub("Ki.67","Ki-67",variable,fixed = T)
	variable = gsub("^S6K$","p70S6K",variable,fixed = F)
	variable = gsub("MKK3.MKK6","MKK36",variable,fixed = T)
	variable = gsub("X.GSK3B.PI3K.PIP3","!GSK3B+PI3K.PIP3",variable,fixed = T)
	variable = gsub("X.ERK12.EGFR.RAS","!ERK12+EGFR.RAS",variable,fixed = T)
	variable = gsub("X.EGFR_FB.EGF.EGFR","!EGFR_FB+EGF.EGFR",variable,fixed = T)
	variable = gsub("X.AKT.RAS.MEK12_S221","!AKT+RAS.MEK12_S221",variable,fixed = T)

	variable = gsub("^ERK[12]*$","ERK12",variable,fixed = F)
	variable = gsub("^MEK[12]*$","MEK12_S221",variable,fixed = F)
	# should be fixed for all possible features...

	# now add column type: which tells us if it is an edge or node related feature
	variable_type = rep("node",length(variable))  # set all to node then change for edges
	variable_type[grep(".",variable,fixed = T)] = "edge"  # find the . between variable names -> edge

	variable = gsub(".","=",variable,fixed = T)  # replace the . with _k_ characters for the edges


	# check if all nodes are corectly included:
	# print nodes, which are neither in the model nor in the data.

	problematic_nodes = variable[variable_type=="node" & !variable %in% modelled_node_names]  # should be empty otherwise check which names to change
	problematic_edges = variable[variable_type=="edge" & !variable %in% modelled_reac_names]  # should be empty otherwise check which names to change

	if(length(problematic_nodes) > 0 | length(problematic_edges)>0) stop(paste("edge or node naming issue with variable:",problematic_nodes,problematic_edges))

	return(variable)
}


# takes a reacID encoding an AND gate interaction and splits it to seperate interactions
resolve_complex_edge <- function(ANDreacID,AND_id){

	sides = strsplit(ANDreacID,"=")[[1]]
	left_terms = strsplit(sides[[1]],"+",fixed = T)[[1]]

	and_node = paste0("and",AND_id)

	new_reacs = paste0(left_terms,"=",and_node)
	new_reacs = c(new_reacs,paste0(and_node,"=",sides[[2]]))

}

# synchronise features an modelled edges/nodes
merge_features_and_model <- function(featureMatrix,modelled_node_names,modelled_reac_names){

	# merge in the nodes that we measire in the network
	featureMatrix = merge(featureMatrix,data.frame(nodename=modelled_node_names),by.x = "variable",by.y = "nodename",all = T)

	N_and_gates = length(grep("+",modelled_reac_names, fixed = T))
	if(N_and_gates>0){
		featureMatrix = merge(featureMatrix,data.frame(nodename=paste0("and",1:N_and_gates)),by.x = "variable",by.y = "nodename",all = T)
	}
	# merge in the features the modelled network: we need for visualisation
	featureMatrix = merge(featureMatrix,data.frame(reacID = modelled_reac_names),by.x = "variable",by.y = "reacID",all = T)

}


# takes the feature Matrix -> convert complex AND gates to individual edges and
# and creates a suitable output for tidy graph
# inputs:
# featureMatrix: data.frame, req. columns: variable, performance
get_drug_edges <- function(featureMatrix,drugName){

	if(!all(c("variable","performance") %in% colnames(featureMatrix))) stop("feature matrix missing required columns")

	# find reacID's encoding complex edges
	index_complex = grep("+",featureMatrix$variable,fixed = T)
	count_ANDs = 0
	reacs_to_add = data.frame()
	for(ind in index_complex){
		count_ANDs = count_ANDs + 1
		reacs = resolve_complex_edge(featureMatrix$variable[ind], AND_id=count_ANDs)

		# replicate the statistics of the edge and then insert the new names
		tmp_df = do.call("rbind",replicate(3,featureMatrix[ind,],simplify = F))
		tmp_df$variable = reacs
		# remove old
		reacs_to_add = rbind(reacs_to_add,tmp_df)
	}
	# remove the complex
	featureMatrix = featureMatrix[-index_complex,]
	# add all new:
	featureMatrix = rbind(featureMatrix,reacs_to_add)

	# separate nodes and edges:
	featureMatrix$type = "node"  # set all to node then change for edges
	featureMatrix$type[grep("=",featureMatrix$variable,fixed = T)] = "edge"  # find the . between variable names -> edge

	# create edge table:
	edge_to_export = featureMatrix %>% as_tibble() %>%
		filter(type=="edge") %>%
		separate(variable,c("from","to"),"=", remove = FALSE) %>%  # edgeid -> from/to nodes
		mutate(drug=drugName) %>%  # add drug name (determined from filename)
		mutate(sign = ifelse(grepl("!",from),"-","+")) %>%  # find inhibitory edged
		mutate(from = gsub("!","",from)) %>%  # inhibition is stored in sign column
		#mutate(stat = ifelse(is.na(stat),"-",stat)) %>%  # remove NA for later issues
		mutate(performance = ifelse(is.na(performance),0,performance)) %>% # remove NA for later issues
		#mutate(abs_performance = abs(performance)) %>% # absolute performance for edge weight
		#mutate(p_val = ifelse(is.na(p_val),1,p_val)) %>%
		select(-variable,-type)

}

# takes the feature Matrix -> gets node table for tidygraph
# inputs:
# featureMatrix: data.frame, req. columns: variable, performance
get_drug_nodes <- function(featureMatrix,drugName,dataNodes){

	if(!all(c("variable","performance") %in% colnames(featureMatrix))) stop("feature matrix missing required columns")

	# separate nodes and edges:
	featureMatrix$type = "node"  # set all to node then change for edges
	featureMatrix$type[grep("=",featureMatrix$variable,fixed = T)] = "edge"  # find the . between variable names -> edge

	node_to_export= featureMatrix %>% as_tibble() %>%
		filter(type=="node") %>%
		mutate(nodetype = ifelse(variable %in% dataNodes,"measured",  # it is either measured, non-measured or AND gate
								 ifelse(grepl("^and[0-9]*$",variable),"and_gate","non-measured"))) %>%
		mutate(drug=drugName) %>% # add drug name
		#mutate(stat = ifelse(is.na(stat),"-",stat)) %>%
		mutate(performance = ifelse(is.na(performance),0,performance)) %>%
		mutate(performance = ifelse(nodetype!="measured",0,performance)) %>%  # performance of non-measured nodes are set to missing --> appears in white
		#mutate(p_val = ifelse(is.na(p_val),1,p_val)) %>%
		select(-type)

}
