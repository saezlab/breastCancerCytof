# plot the CNOlist files for each cell lines

library(multiCellNOpt)
library(plyr)
library(dplyr)

# calibrated_model_folder = "./data/models/pkn_v2_midas_v2/outputs/"
# calibrated_model_folder = "./data/models/pkn_v3_midas_v2/outputs/"
# calibrated_model_folder = "./data/models/pkn_v4_midas_v2/outputs/"
#calibrated_model_folder = "./data/models/pkn_v4_midas_v3_bootstrap/outputs/"
calibrated_model_folder = "./data/models/pkn_v4_midas_v4/outputs/"


calibrated_model_files = list.files(calibrated_model_folder,"*.RDS",full.names = T)

calibrated_models = lapply(calibrated_model_files,function(f){
	m = readRDS(f)
	# update the class definition to the newest version
	#m2 = logicODEModel$new(oldModel = m)
}
)

names(calibrated_models) = gsub(".RDS","",basename(calibrated_model_files))

# pdf("figures/fitted_models_plotFit_pkn_v3_midas_v2.pdf",width = 21,height = 3)
# pdf("figures/fitted_models_plotFit_pkn_v4_midas_v2.pdf",width = 21,height = 3)

pdf("figures/fitted_models_convergence_pkn_v4_midas_v4.pdf",width = 7,height = 7)
l_ply(seq_along(calibrated_models),function(i,M){   #
	#i=1
	#M = calibrated_models
	plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
		 xaxt='n', yaxt='n', xlab='', ylab='')
	text(1,4,paste0(i,"-",M[[i]]$exps$name), pos=4)

	if(is.null(M[[i]]$fitResults)) return()

	M[[i]]$reportConvergence()

},M = calibrated_models,.progress = progress_text())
dev.off()


pdf("figures/fitted_models_plotFit_pkn_v4_midas_v4.pdf",width = 21,height = 3)
l_ply(seq_along(calibrated_models),function(i,M){
	M[[i]]$plotFit(measuredNodesOnly = F)
	title(paste0(i,"-",M[[i]]$exps$name))
},M = calibrated_models,.progress = progress_text())
dev.off()

pdf("figures/fitted_models_plotFit_obsonly_pkn_v4_midas_v4.pdf",width = 21,height = 3)
l_ply(seq_along(calibrated_models),function(i,M){
	M[[i]]$plotFit(measuredNodesOnly = T)
	title(paste0(i,"-",M[[i]]$exps$name))
},M = calibrated_models,.progress = progress_text())
dev.off()
