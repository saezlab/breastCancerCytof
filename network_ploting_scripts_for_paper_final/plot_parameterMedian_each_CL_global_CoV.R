# updates

# - 24 th July:
# we got back to the best fitting parameters: the drug prediction is re-ran on the
# new GDSC data where we used the best parameters from the MFU sampling
# - 28 of May:
#    instead of the best fitting parameter vector we adapt the pipeline for the
#  best optimised parameter vector. This might have worse fit for some models.
#	par_stats.RDS got updated today with the new stats
# - no filtering for accurate parameters?


# Summary script to produce network figures:

library(tidyverse)
library(tidygraph)
library(RCy3)

par_stats = readRDS("./data/models/pkn_v4_midas_v4/par_stats.RDS") %>% as_tibble()

get_edge_data <- function(par_stats, stat = c("CL_median","CL_mean","CL_best","CL_optim")){

	stat = match.arg(stat)
	# 1. arrange the cell_line data into column --> not tidy
	# 1.2 edges
	edges <- par_stats %>%
		filter(grepl("_k_",parameter)) %>%
		separate(parameter,into = c("from","to"),"_k_", remove = FALSE)

	# CL_median for edges:
	new_names = unique(edges$cell_line)
	names(new_names) = paste0(stat,"_",new_names)
	edges_stat <- edges %>%
		select(cell_line,parameter,from,to,stat) %>%
		spread(cell_line,stat) %>%
		rename(!!new_names)

}

get_node_data <- function(par_stats,stat = c("CL_median","CL_mean","CL_best","CL_optim") ){

	stat = match.arg(stat)

	# 2. add information on measurements and inputs
	cno = CellNOptR::CNOlist("./data/MIDAS_v4/184A1.csv")
	dataNodes = colnames(cno@signals[[1]])

	# 1.1 nodes
	nodes <-  par_stats %>%
		filter(grepl("tau_",parameter)) %>%  # remove x0: keep only tau and k parameters
		mutate(name = gsub("tau_","",parameter))

	# CL_median table for nodes:
	new_names = unique(nodes$cell_line)
	names(new_names) = paste0(stat,"_",new_names)

	nodes_stat <- nodes %>%
		select(cell_line,stat,name) %>%
		spread(cell_line,stat) %>%
		rename(!!new_names)


	nodes_stat <- nodes_stat %>%
		mutate(nodetype = ifelse(nodes_stat$name %in% dataNodes,"measured","non-measured"))

	nodes_stat <- nodes_stat %>% full_join(y = data.frame(name=c("EGF","SERUM"),nodetype=c("input","input"),stringsAsFactors=F),by = c("name","nodetype"))
}

# From par_stats generate the CL_median data for each edge/node in a table, where
# each column is for a different cell_line.
#nodes_stat = get_node_data(par_stats,stat = "CL_mean") # until May 25
#nodes_stat = get_node_data(par_stats,stat = "CL_optim") # befoire re-running macau
#edges_stat = get_edge_data(par_stats,stat = "CL_optim")

nodes_stat = get_node_data(par_stats,stat = "CL_best")
edges_stat = get_edge_data(par_stats,stat = "CL_best")

# # renaming nodes: shorten the node names:
target_names <- c("SMAD23", "MAPKAPK2", "MEK12", "ERK12", "MEK12_S221", "p70S6K", "MSK12")
names(target_names) <-  c("SMAD", "MK2", "MEK", "ERK", "MEK_S221", "S6K","MSK")

# i=1
for(i in seq_along(target_names)){
	nodes_stat$name[nodes_stat$name==target_names[[i]]] = names(target_names)[[i]]
	edges_stat$from = gsub(target_names[[i]],names(target_names)[[i]], edges_stat$from,fixed = T)
	edges_stat$to = gsub(target_names[[i]],names(target_names)[[i]], edges_stat$to,fixed = T)
	edges_stat$parameter = gsub(target_names[[i]],names(target_names)[[i]], edges_stat$parameter,fixed = T)
}



# Add Global variables
# renaming nodes: shorten the node names:
par_stats$parameter = as.character(	par_stats$parameter)
for(i in seq_along(target_names)){
	par_stats$parameter = gsub(target_names[[i]],names(target_names)[[i]], par_stats$parameter,fixed = T)
}

par_stats %>% mutate( type = ifelse(grepl("_k_",parameter),"k",	ifelse(grepl("tau_",parameter),"tau","x0"))) %>%
	ggplot(.) + geom_point(aes(parameter,CL_cov)) + facet_wrap(~type,scales = "free")

par_stats %>% mutate( type = ifelse(grepl("_k_",parameter),"k",	ifelse(grepl("tau_",parameter),"tau","x0"))) %>%
	ggplot(.) + geom_boxplot(aes(parameter,CL_cov)) + facet_wrap(~type,scales = "free_x") + geom_hline(yintercept = 0.5)



# up to may 25 2019
# global_stats <- par_stats %>%
# 	mutate(accurate = CL_cov < 0.5) %>%
# 	group_by(parameter) %>%
# 	filter(accurate,.preserve = T) %>%
# 	summarise(CoV_par = sd(CL_mean)/mean(CL_mean), n_celline = sum(accurate)) %>%
# 	mutate(CoV_par = ifelse(is.na(CoV_par),0,CoV_par))   # if the mean was 0, then COV is NA, but then it is accurate.

# we are no longer filtering, and we use the optimised parameter vector
# global_stats <- par_stats %>%
# 	mutate(accurate = TRUE) %>%
# 	group_by(parameter) %>%
# 	#filter(accurate,.preserve = T) %>%
# 	summarise(CoV_par = sd(CL_optim)/mean(CL_optim), n_celline = sum(accurate)) %>%
# 	mutate(CoV_par = ifelse(is.na(CoV_par),0,CoV_par))   # if the mean was 0, then COV is NA, but then it is accurate.

# Since Macau rerun: we are no longer filtering, and we use the Best parameter vector
global_stats <- par_stats %>%
	mutate(accurate = TRUE) %>%
	group_by(parameter) %>%
	#filter(accurate,.preserve = T) %>%
	summarise(CoV_par = sd(CL_best)/mean(CL_best),
			  n_celline = sum(accurate),
			  CL_mean_of_bests = mean(CL_best)) %>%
	mutate(CoV_par = ifelse(is.na(CoV_par),0,CoV_par))   # if the mean was 0, then COV is NA, but then it is accurate.



global_edges <- global_stats %>%
	filter(grepl("_k_",parameter)) %>%
	separate(parameter,into = c("from","to"),"_k_", remove = FALSE)


global_nodes <- global_stats %>%
	filter(grepl("tau_",parameter)) %>%  # keep only tau
	mutate(name = gsub("tau_","",parameter)) %>% select(-parameter)



## Merge 2 tables for same layout
all_nodes <- nodes_stat %>% full_join(global_nodes,by = "name")
all_edges <- edges_stat %>% full_join(global_edges,by = c("parameter" ,"from","to"))


# Set the inputs: EGF and SERUM values to -1 --> coloring them separately
# all_nodes_tmp = all_nodes
all_nodes[is.na(all_nodes)] = -1
# 1. create an igraph object
G  <-  tbl_graph(nodes = data.frame(all_nodes),edges = data.frame(all_edges), directed = TRUE)


########################################################################################
########################################################################################
##### !!!  START CYTOSCAPE

createNetworkFromIgraph(G,title = "all plot")
### Export indivdual networks --------------------------------------------------
visual_style <- "CC_v_drug" #"CC_v1"
setVisualStyle(visual_style)
layoutNetwork(layout.name = "apply preferred")
#layoutNetwork(layout.name = "hierarchical")

#cols = RColorBrewer::brewer.pal(2,"Reds")
# the first color is used to color input nodes, the other two for building a gradient
cols = gplots::col2hex(c("chartreuse3","grey90","firebrick3"))
#edge_cols = RColorBrewer::brewer.pal(3,"Set1")
cell_lines = unique(par_stats$cell_line)


# For each cell-line edge strength:
# i_celline = 1
bar = progress::progress_bar$new(total = length(cell_lines))
for(i_celline in 1:length(cell_lines)){
	bar$tick()
	cl_name = cell_lines[[i_celline]]
#	print(cl_name)

	var_name = paste0("CL_best","_",cl_name)
	# colors: first is a distinct to color input nodes, second and thirds used to create a gradient
	setEdgeColorMapping(var_name,c(-1,-0.001,5),cols,style.name = visual_style); Sys.sleep(0.3)
	setNodeColorMapping(var_name,c(-1,-0.001,5),cols,style.name = visual_style); Sys.sleep(0.3)
	setEdgeLineWidthMapping(var_name,c(0,5),widths = c(8,30),style.name = visual_style); Sys.sleep(0.3)
	#setNodeBorderColorMapping("nodetype",table.column.values = c("non-measured","measured","input"),colors = edge_cols,style.name = "CC_v1",mapping.type = "d"); Sys.sleep(0.1)

	layoutNetwork(layout.name = "apply preferred"); Sys.sleep(1)

	exportImage(filename = file.path(getwd(),"./data/models/pkn_v4_midas_v4/networks/best_par_based",paste0(var_name,".pdf")),type = "PDF")
}
# we remove the -1 for the coloring to generate the LEGEND (manually)
setEdgeColorMapping(var_name,c(0,5),cols[-1],style.name = visual_style); Sys.sleep(0.3)
setNodeColorMapping(var_name,c(0,5),cols[-1],style.name = visual_style); Sys.sleep(0.3)

### EXPORT LEGEND MANUALLY!!!!!

## global parameter COV_par (covariance of parameters), CL_mean_of_bests (mean parameter across cellslings)
###  export COV_par network ----------------------------------------------------
stat_name = "CoV_par"
cols = gplots::col2hex(c("grey90","navyblue"))
cols = gplots::col2hex(c("chartreuse3","grey90","orange"))
#cols[[2]] = "#0099FF"

var_name = stat_name
setEdgeColorMapping(var_name,
					c(-1,quantile(c(all_nodes$CoV_par,all_edges$CoV_par),probs = c(0,0.1,0.9,1),na.rm = T,names = F)),
					c(cols[1:2],cols[2:3],cols[[3]]),style.name = visual_style); Sys.sleep(0.1)
setNodeColorMapping(var_name,
					c(-1,quantile(c(all_nodes$CoV_par,all_edges$CoV_par),probs = c(0,0.1,0.9,1),na.rm = T,names = F)),
					c(cols[1:2],cols[2:3],cols[[3]]),style.name = visual_style); Sys.sleep(0.1)
setEdgeLineWidthMapping(var_name,
						quantile(c(all_nodes$CoV_par,all_edges$CoV_par),probs = c(0,0.1,0.9,1),na.rm = T,names = F),
						widths = c(8,8,30,30),style.name = visual_style); Sys.sleep(0.1)


layoutNetwork(layout.name = "apply preferred"); Sys.sleep(1)

exportImage(filename = file.path(getwd(),"./data/models/pkn_v4_midas_v4/networks/best_par_based/",paste0(var_name,"_orange.pdf")),type = "PDF")

# Correcting legend: removing coloring hacks
setEdgeColorMapping(var_name,
					c(quantile(c(all_nodes$CoV_par,all_edges$CoV_par),probs = c(0,0.1,0.9,1),na.rm = T,names = F)),
					c(cols[2],cols[2:3],cols[[3]]),style.name = visual_style); Sys.sleep(0.1)

setNodeColorMapping(var_name,
					c(quantile(c(all_nodes$CoV_par,all_edges$CoV_par),probs = c(0,0.1,0.9,1),na.rm = T,names = F)),
					c(cols[2],cols[2:3],cols[[3]]),style.name = visual_style); Sys.sleep(0.1)

# EXPORT LEGEND NOW !!!

stat_name = "CL_mean_of_bests"
# cols = gplots::col2hex(c("chartreuse3","grey90","orange"))
cols = gplots::col2hex(c("chartreuse3","grey90","firebrick3"))
#cols[[2]] = "#0099FF"

var_name = stat_name
setEdgeColorMapping(var_name,
					c(-1,quantile(c(all_nodes$CL_mean_of_bests,all_edges$CL_mean_of_bests),probs = c(0,0.1,0.9,1),na.rm = T,names = F)),
					c(cols[1:2],cols[2:3],cols[[3]]),style.name = visual_style); Sys.sleep(0.1)

setNodeColorMapping(var_name,
					c(-1,quantile(c(all_nodes$CL_mean_of_bests,all_edges$CL_mean_of_bests),probs = c(0,0.1,0.9,1),na.rm = T,names = F)),
					c(cols[1:2],cols[2:3],cols[[3]]),style.name = visual_style); Sys.sleep(0.1)

setEdgeLineWidthMapping(var_name,
						quantile(c(all_nodes$CL_mean_of_bests,all_edges$CL_mean_of_bests),probs = c(0,0.1,0.9,1),na.rm = T,names = F),
						widths = c(8,8,30,30),style.name = visual_style); Sys.sleep(0.1)
# setEdgeColorMapping(var_name,c(-1,-0.001,5),cols,style.name = visual_style); Sys.sleep(0.3)
# setNodeColorMapping(var_name,c(-1,-0.001,5),cols,style.name = visual_style); Sys.sleep(0.3)
# setEdgeLineWidthMapping(var_name,c(0,5),widths = c(8,30),style.name = visual_style); Sys.sleep(0.3)



layoutNetwork(layout.name = "apply preferred"); Sys.sleep(1)

exportImage(filename = file.path(getwd(),"./data/models/pkn_v4_midas_v4/networks/best_par_based/",paste0(var_name,".pdf")),type = "PDF")

#exportImage(filename = file.path(getwd(),"supp_info/results_for_april_2019/network_figs",paste0(var_name,"_plain.pdf")),type = "PDF")



# Correcting legend: removing coloring hacks
setEdgeColorMapping(var_name,
					c(quantile(c(all_nodes$CoV_par,all_edges$CL_mean_of_bests),probs = c(0,0.1,0.9,1),na.rm = T,names = F)),
					c(cols[2],cols[2:3],cols[[3]]),style.name = visual_style); Sys.sleep(0.1)

setNodeColorMapping(var_name,
					c(quantile(c(all_nodes$CoV_par,all_edges$CL_mean_of_bests),probs = c(0,0.1,0.9,1),na.rm = T,names = F)),
					c(cols[2],cols[2:3],cols[[3]]),style.name = visual_style); Sys.sleep(0.1)

# EXPORT LEGEND NOW !!!
