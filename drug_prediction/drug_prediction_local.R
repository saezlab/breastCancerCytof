# find significant features:
library(plyr)
library(dplyr)
## User inputs ------
IC50_cellline_matrix_file = "./data/DATA_GDSC/GDSC_IC50_breast"
IC50_matrix = read.csv(IC50_cellline_matrix_file)

feature_data_file  = "./data/models/pkn_v4_midas_v4/features/feature_table_1_raw_parameters.RData"
load(feature_data_file)

feature_data = feature_table_1[2:129]
measured_cell_lines = as.character(feature_table_1$cell_line)
rownames(feature_data) = measured_cell_lines


# get only measured cell lines from drugs and get only drugged cell lines from measurements
drug_cell_lines = colnames(IC50_matrix)[2:49]

common_cell_lines = intersect(measured_cell_lines, drug_cell_lines)
common_IC50_matrix = IC50_matrix[,c("X",common_cell_lines)]
common_features = feature_data[common_cell_lines,]


interesting_drugs = c("Paclitaxel","Bortezomib","CP466722","Phenformin","AS605240","Genentech Cpd 10","Masitinib", "PIK-93","SB52334","Camptothecin","Vinblastine","Docetaxel","TW 37")
interesting_drugs_ids = which(common_IC50_matrix$X %in% interesting_drugs)


## User inputs END ------


## Build a sparse model for each drug multiple times from the feaures to predict IC50

#' analyse drug
#' @param cell_lines character vector of cell line names (length C)
#' @param features C*F matrix, for each cell line (in a row) the F number of features
#' @param drug_response length C drug response IC05 of the cell lines
analyse_drug = function(cell_lines,features,drug_response){

	stopifnot((length(cell_lines) == nrow(features)) & (length(cell_lines) == length(drug_response)))


	numeric_id = which(!is.na(drug_response))
	cell_lines = cell_lines[numeric_id]
	features = features[numeric_id,]
	drug_response = drug_response[numeric_id]
	drug_response_random = sample(drug_response)

	# standardize tau and K parameters
	features = as.matrix(features)
	ix_tau<-grep("tau_", colnames(features))
	ix_k<-grep("_k_", colnames(features))

	mean_k<-mean(features[,ix_k])
	sd_k<-sd(features[,ix_k])

	mean_tau<-mean(features[,ix_tau])
	sd_tau<-sd(features[,ix_tau])

	features_bk<-features

	features[,ix_tau]<-(features[,ix_tau]-mean_tau)/sd_tau
	features[,ix_k]<-(features[,ix_k]-mean_k)/sd_k


	# we do a systematic leave one out prediction
	# theIteration = 6
	library(caret)
	#res = plyr::ldply(1:length(drug_response), function(theIteration){
	res = plyr::ldply(1:4, function(theIteration){
		print(paste0("iter: ",theIteration))

		training_data = drug_response[-theIteration]
		training_data_rnadom = drug_response_random[-theIteration]
		test_data = drug_response[theIteration]
		training_features = features[-theIteration,]
		test_features = features[theIteration,]

		# do leave one out cross validation on the remaining data (using caret)
		trC = caret::trainControl(method="LOOCV")
		fitM = caret::train(training_features, training_data, trControl = trC, method="glmnet", tuneLength = 20,)
		fitM_rm = train(training_features, training_data_rnadom, trControl = trC, method="glmnet", tuneLength = 20)

		# save best rmse
		rmse_best<-min(fitM$results$RMSE)
		rmse_best_rm<-min(fitM_rm$results$RMSE)


		cat(theIteration, ": alpha=", round(fitM$bestTune$alpha,2),
			", lambda=", round(fitM$bestTune$lambda,2),
			", RMSE=", round(rmse_best, 2),
			", RMSE random=", round(rmse_best_rm, 2),
			"\n", sep="")

		# save estimated parameters
		a<-as.matrix(coef(fitM$finalModel, s=fitM$bestTune$lambda))
		#a<-a[allLinks.names,]
		a<-c(a, RMSE=rmse_best, RMSE_rm=rmse_best_rm)
		return(a)
	})
	return(res)

}




library("doParallel")
cl <- makeCluster(3)  ## creates 48 virtual clusters
registerDoParallel(cl)


idrug=interesting_drugs_ids[[1]]

par_res = llply(interesting_drugs_ids[1:4], function(idrug){

	drug_IC50 = common_IC50_matrix[idrug,]
	drug_IC50 = reshape2::melt(drug_IC50,id.vars = "X",variable.name = "cell_line",value.name = "IC50")
	rownames(drug_IC50) = drug_IC50$cell_line
	drug_IC50$cell_line = NULL

	analyse_drug(cell_lines = common_cell_lines, features = common_features, drug_response = drug_IC50[,"IC50"])


},.parallel = TRUE,.paropts = list(.packages=c("caret","reshape2"),.export=c("common_IC50_matrix","common_cell_lines","common_features","analyse_drug")))



for(idrug in 1:nrow(common_IC50_matrix)){


	drug_IC50 = common_IC50_matrix[idrug,]
	drug_IC50 = melt(drug_IC50,id.vars = "X",variable.name = "cell_line",value.name = "IC50")
	rownames(drug_IC50) = drug_IC50$cell_line
	drug_IC50$cell_line = NULL

	analyse_drug(cell_lines = common_cell_lines, features = common_features, drug_response = drug_IC50[,"IC50"])
}





