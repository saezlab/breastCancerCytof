# compare model fits between versions

library(multiCellNOpt)
library(plyr)
library(dplyr)

calibrated_model_folder = c("./data/models/pkn_v4_midas_v2/outputs/",
							"./data/models/pkn_v2_midas_v2/outputs/")

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


modelsA = import_models(calibrated_model_folder[[1]])
modelsB = import_models(calibrated_model_folder[[2]])

modelsA_stats = ldply(modelsA,function(m)m$getStatistics())
modelsB_stats = ldply(modelsB,function(m)m$getStatistics())


modelsA_stats$model_version = "pkn_v4_midas_v2"
modelsB_stats$model_version = "pkn_v2_midas_v2"

models_stats = rbind(modelsA_stats,modelsB_stats)

ggplot(models_stats,aes(model_version,RMSE)) +
	geom_violin(aes(fill=model_version)) +
	geom_boxplot(width=0.3) +
	geom_jitter(col="darkblue",width = 0.1) +
	ylim(0,0.15) +
	theme_bw() +
	ggtitle("Comparing fitting performance by RMSE") +
	xlab("model version")

