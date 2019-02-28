# Associate the feature matrix (inputs) with drug response by fitting multiple elastic net models.

# inputs:
# args[[1]]: IC50 matrix
# args[[2]]: feature matrix
# args[[3]]: drug_id

args = commandArgs(trailingOnly=TRUE)

# test with custom argument:
# args = list()
# args[[1]] = "./data/models/pkn_v4_midas_v4/features_vs_drug/inputs/common_IC50_matrix.RDS"
# args[[2]] = "./data/models/pkn_v4_midas_v4/features_vs_drug/inputs/common_rawPar_features.RDS"
# args[[3]] = 49

if (length(args)<3) {
	stop("At least three argument must be supplied (IC50, feature matrix file and drug_id)", call.=FALSE)
}

library(plyr)
library(dplyr)
library(caret)
library(reshape2)


IC50_file = args[[1]]
feature_data_file = args[[2]]
drug_id = as.numeric(args[[3]])


common_IC50_matrix = readRDS(IC50_file)
common_features = readRDS(feature_data_file)


if(drug_id>nrow(common_IC50_matrix)) stop("invalida drug_id. Must be in range 1-256.")
if(!class(common_IC50_matrix[,1]) %in% c("factor","character")) stop("IC50 matrix should be a dataframe, first column with drug names.")

feature_name = gsub(".RDS","",basename(feature_data_file),fixed = T)

if(drug_id==0) print("drug_id==0 Testing...")
## Build a sparse model for each drug multiple times from the feaures to predict IC50

#' analyse drug
#' @param cell_lines character vector of cell line names (length C)
#' @param features C*F matrix, for each cell line (in a row) the F number of features
#' @param drug_response length C drug response IC05 of the cell lines
analyse_drug = function(features,drug_response){

	stopifnot( nrow(features) == length(drug_response))

	numeric_id = which(!is.na(drug_response))
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

	# res = plyr::ldply(1:length(drug_response), function(theIteration){
	res = plyr::ldply(1:1, function(theIteration){
		print(paste0("iter: ",theIteration))

		training_data = drug_response[-theIteration]
		training_data_rnadom = drug_response_random[-theIteration]
		test_data = drug_response[theIteration]
		training_features = features[-theIteration,]
		test_features = features[theIteration,]

		# do leave one out cross validation on the remaining data (using caret)
		trC = caret::trainControl(method="LOOCV")
		fitM = caret::train(training_features, training_data, trControl = trC, method="glmnet", tuneLength = 20)
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

if(drug_id!=0){

	drug_IC50 = common_IC50_matrix[drug_id,]
	colnames(drug_IC50)[[1]]="drug"
	drug_IC50 = reshape2::melt(drug_IC50,id.vars = "drug",variable.name = "cell_line",value.name = "IC50")
	rownames(drug_IC50) = drug_IC50$cell_line
	drug_IC50$cell_line = NULL


	drug_res = analyse_drug(features = common_features, drug_response = drug_IC50[,"IC50"])

}else{
	drug_res = data.frame(a="test",b="success")
}

if(!dir.exists("./outputs")) dir.create("./outputs")
saveRDS(drug_res, paste0("./outputs/predict_drug_id_",drug_id,"_",feature_name,".RDS"))

