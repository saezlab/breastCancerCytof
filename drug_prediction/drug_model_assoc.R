# compare the models similarity with the drug response.
# Let someone more talented in drug prediction do this....
library(reshape2)
calibrated_model_files = list.files("./data/models/pkn_v4_midas_v2/outputs/","*.RDS",full.names = T)

calibrated_models = lapply(calibrated_model_files,function(f){
	m = readRDS(f)
	# update the class definition to the newest version
	#m2 = logicODEModel$new(oldModel = m)
}
)

names(calibrated_models) = gsub(".RDS","",gsub("early_model_","",basename(calibrated_model_files)))



drug_response = readxl::read_excel("data/v17.3_fitted_dose_response.xlsx",sheet = 1)



### Find cell lines in GDSC  ------------------
# Fix naming
# no nanmes in the model contains space or dash"-"
unique(drug_response$CELL_LINE_NAME)
drug_response$CELL_LINE_NAME_new = gsub("-","",drug_response$CELL_LINE_NAME)
drug_response$CELL_LINE_NAME_new = gsub(" ","",drug_response$CELL_LINE_NAME_new)


drug_response  = filter(drug_response, CELL_LINE_NAME_new %in% names(calibrated_models))

tmp = filter(drug_response, DRUG_NAME=="Rapamycin")
tmp = filter(drug_response, CELL_LINE_NAME_new=="MFM223")


tmp = tmp[order(tmp$CELL_LINE_NAME_new),]


ggplot(filter(drug_response,DRUG_NAME %in% sample(unique(drug_response$DRUG_NAME),25))) + geom_point(aes(CELL_LINE_NAME_new,LN_IC50)) + coord_flip() + facet_wrap(~DRUG_NAME)




#### Try to cluster cell lines based on drug respone
# we need to get rid of NA's




drug_response_df = dcast(as.data.frame(drug_response), DRUG_NAME~CELL_LINE_NAME_new, value.var = "LN_IC50", fun.aggregate = function(x){max(x)} )

rownames(drug_response_df) = drug_response_df$DRUG_NAME
drug_response_df$DRUG_NAME = NULL
missing_IC50_for_Drug = sort(colSums(t(drug_response_df)==-Inf))
drugs_to_remove = names(missing_IC50_for_Drug[missing_IC50_for_Drug > 10])
drug_response_df = drug_response_df[!rownames(drug_response_df) %in% drugs_to_remove,]



sort(colSums(drug_response_df==-Inf))
missing_IC50_for_CellLines = sort(colSums(drug_response_df==-Inf))
celllines_to_remove = names(missing_IC50_for_CellLines[missing_IC50_for_CellLines > 10])
drug_response_df = drug_response_df[,!colnames(drug_response_df) %in% celllines_to_remove]




image(drug_response_df==-Inf )


# inpute the -INF data by mean over other cell lines\
# irow = 2
for(irow in 1:nrow(drug_response_df)){

	if(sum(drug_response_df[irow,]==-Inf)>0){
		drug_response_df[irow,drug_response_df[irow,]==-Inf] = NA
		drug_response_df[irow,is.na(drug_response_df[irow,])] = 	mean(as.numeric(drug_response_df[irow,]),na.rm=TRUE)
	}


}
pheatmap::pheatmap(drug_response_df,clustering_distance_rows = "manhattan",clustering_distance_cols ="manhattan",fontsize_row = 4 )



#### Find extreme drug responses:


# scale by range for each drug
drug_response_df_range = adply(drug_response_df,1,range)
hist(drug_response_df_range$V2-drug_response_df_range$V1, main=" histogram of IC50 range per cell lines", xlab="IC50 range (max-min)")

# around 25 extreme cases where the range is larger than 6 order of magnitude.
extremes = which((drug_response_df_range$V2-drug_response_df_range$V1) > 6) # 6
drug_response_df[extremes,]
library(RColorBrewer)
extreme_drug_resp_df = drug_response_df[extremes,]

extreme_cell_lines = colnames(extreme_drug_resp_df)

pheatmap::pheatmap(extreme_drug_resp_df,clustering_distance_rows = "correlation",clustering_distance_cols ="manhattan",
				   color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(7),breaks = seq(-13, 13,length.out = 8),fontsize_row = 6)



extreme_response_drugnames = rownames(extreme_drug_resp_df)



drug_target = unique(filter(drug_response, DRUG_NAME %in% extreme_response_drugnames)[,c("DRUG_NAME","PUTATIVE_TARGET")])
drug_target = as.data.frame(drug_target)
drug_target = drug_target[order(drug_target$PUTATIVE_TARGET),]

indx.match = match(extreme_response_drugnames,drug_target$DRUG_NAME)
drug_target[indx.match,]



pheatmap::pheatmap(extreme_drug_resp_df,clustering_distance_rows = "correlation",clustering_distance_cols ="manhattan",
				   color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(7),breaks = seq(-13, 13,length.out = 8),
				   labels_row = drug_target[indx.match,"PUTATIVE_TARGET"],fontsize_row = 6)



### Model features 1: parameters -----------------------------------------------
library(multiCellNOpt)
library(plyr)
library(dplyr)

calibrated_model_folder = c("./data/models/pkn_v4_midas_v2/outputs/")

# imports and names the models
import_models = function(folder){

	calibrated_model_files = list.files(folder,"*.RDS",full.names = T)

	calibrated_models = lapply(calibrated_model_files,function(f){
		m = readRDS(f)
	}
	)

	names(calibrated_models) = gsub(".RDS","",basename(calibrated_model_files))
	return(calibrated_models)
}


models = import_models(calibrated_model_folder)
models = models[extreme_cell_lines]

# get estimated model parameters --------------------

# all
model_parameters = ldply(models, function(m)m$ode_parameters$parValues[m$ode_parameters$index_opt_pars])
rownames(model_parameters) = model_parameters$.id
model_parameters$.id = NULL

pheatmap::pheatmap(t(model_parameters),clustering_distance_rows = "manhattan",clustering_distance_cols ="manhattan",
				   color = colorRampPalette(brewer.pal(n = 7, name ="Reds"))(7),breaks = seq(0, 10,length.out = 8))

# k
model_parameters_k = ldply(models, function(m)m$ode_parameters$parValues[m$ode_parameters$index_k])
rownames(model_parameters_k) = model_parameters_k$.id
model_parameters_k$.id = NULL

pheatmap::pheatmap(t(model_parameters_k),clustering_distance_rows = "manhattan",clustering_distance_cols ="manhattan",
				   color = colorRampPalette(brewer.pal(n = 7, name ="Reds"))(7),breaks = seq(0, 10,length.out = 8))
# tau
model_parameters_tau = ldply(models, function(m)m$ode_parameters$parValues[m$ode_parameters$index_tau])
rownames(model_parameters_tau) = model_parameters_tau$.id
model_parameters_tau$.id = NULL

pheatmap::pheatmap(t(model_parameters_tau),clustering_distance_rows = "manhattan",clustering_distance_cols ="manhattan",
				   color = colorRampPalette(brewer.pal(n = 7, name ="Reds"))(7),breaks = seq(0, 10,length.out = 8))



### Focus on MEK/ERK targeting drugs and model parameters

ind = grep("MEK|ERK",drug_target[indx.match,"PUTATIVE_TARGET"])

mek_erk_drug = extreme_drug_resp_df[ind,]

pheatmap::pheatmap(mek_erk_drug,clustering_distance_rows = "correlation",clustering_distance_cols ="manhattan",
				   color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(7),breaks = seq(-10, 10,length.out = 8),
				   labels_row = drug_target[indx.match,"PUTATIVE_TARGET"][ind])

pheatmap::pheatmap(mek_erk_drug,clustering_distance_rows = "correlation",clustering_distance_cols ="manhattan",
				   color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(7),breaks = seq(-10, 10,length.out = 8))

# k
model_parameters_k = ldply(models, function(m)m$ode_parameters$parValues[m$ode_parameters$index_k])
rownames(model_parameters_k) = model_parameters_k$.id
model_parameters_k$.id = NULL




ind2 = grep("MEK|ERK",colnames(model_parameters_k))
model_parameters_k_mekerk = model_parameters_k[,ind2]

hist(log10(unlist(model_parameters_k_mekerk)),20)

model_parameters_k_mekerk = log10(model_parameters_k_mekerk + 1e-5)




pheatmap::pheatmap(t(model_parameters_k_mekerk),clustering_distance_rows = "euclidean",clustering_distance_cols ="euclidean",
				   color = c("#bbbbbb",colorRampPalette(brewer.pal(n = 7, name ="Reds"))(7)))



### Focus on AKT/Pi3K targeting drugs and model parameters

ind = grep("AKT|PI3K|PDK",drug_target[indx.match,"PUTATIVE_TARGET"])

akt_pi3K_pdk_drug = extreme_drug_resp_df[ind,]

pheatmap::pheatmap(akt_pi3K_pdk_drug,clustering_distance_rows = "correlation",clustering_distance_cols ="manhattan",
				   color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(7),breaks = seq(-10, 10,length.out = 8),
				   labels_row = drug_target[indx.match,"PUTATIVE_TARGET"][ind])

pheatmap::pheatmap(akt_pi3K_pdk_drug,clustering_distance_rows = "euclidean",clustering_distance_cols ="manhattan",
				   color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(7),breaks = seq(-10, 10,length.out = 8))

# k
model_parameters_k = ldply(models, function(m)m$ode_parameters$parValues[m$ode_parameters$index_k])
rownames(model_parameters_k) = model_parameters_k$.id
model_parameters_k$.id = NULL




ind2 = grep("AKT|PI3K|PIP3|PDK",colnames(model_parameters_k))
model_parameters_k_akt_pi3k_pdk = model_parameters_k[,ind2]

hist(log10(unlist(model_parameters_k_akt_pi3k_pdk)),20)

model_parameters_k_akt_pi3k_pdk = log10(model_parameters_k_akt_pi3k_pdk + 1e-7)




pheatmap::pheatmap(t(model_parameters_k_akt_pi3k_pdk),clustering_distance_rows = "euclidean",clustering_distance_cols ="euclidean",
				   color = c("#bbbbbb",colorRampPalette(brewer.pal(n = 7, name ="Reds"))(7)))




### Focus on AKT/Pi3K targeting drugs and model parameters

ind = grep("MTOR",drug_target[indx.match,"PUTATIVE_TARGET"])

mtor_drug = extreme_drug_resp_df[ind,]

pheatmap::pheatmap(mtor_drug,clustering_distance_rows = "correlation",clustering_distance_cols ="manhattan",
				   color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(7),breaks = seq(-10, 10,length.out = 8),
				   labels_row = drug_target[indx.match,"PUTATIVE_TARGET"][ind])

pheatmap::pheatmap(mtor_drug,clustering_distance_rows = "euclidean",clustering_distance_cols ="manhattan",
				   color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(7),breaks = seq(-10, 10,length.out = 8))

# k
model_parameters_k = ldply(models, function(m)m$ode_parameters$parValues[m$ode_parameters$index_k])
rownames(model_parameters_k) = model_parameters_k$.id
model_parameters_k$.id = NULL




ind2 = grep("mTOR",colnames(model_parameters_k))
model_parameters_k_mtor = model_parameters_k[,ind2]

hist(log10(unlist(model_parameters_k_mtor)),20)

model_parameters_k_mtor = log10(model_parameters_k_mtor + 1e-5)




pheatmap::pheatmap(t(model_parameters_k_mtor),clustering_distance_rows = "euclidean",clustering_distance_cols ="euclidean",
				   color = c("#bbbbbb",colorRampPalette(brewer.pal(n = 7, name ="Reds"))(7)))
