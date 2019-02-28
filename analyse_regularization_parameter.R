# evaluate the regularisation parameter scanning results for 5 cell lines.
# we launched estimatation problems with regulariosation parameters for parameter k (lambda_k) between 1e-7 and 1.
# Here we plot the results and draw the conclusion that the regularisation should be somewhere
# between 1e-4 and 1e-5.




library(multiCellNOpt)
library(plyr)
library(dplyr)

# import calibrated models
calibrated_model_folder = "./data/models/pkn_v4_midas_v2_regpar/outputs/"
calibrated_model_files = list.files(calibrated_model_folder,"*.RDS",full.names = T)
calibrated_models = lapply(calibrated_model_files,function(f){
	m = readRDS(f)
	# update the class definition to the newest version
	#m2 = logicODEModel$new(oldModel = m)
}
)
names(calibrated_models) = gsub(".RDS","",basename(calibrated_model_files))

# remove 2 cases which were not finished on the cluster
calibrated_models = calibrated_models[!names(calibrated_models) %in% c("HCC2157_lambdaK_3", "HCC1143_lambdaK_4")]


# get model parameters and statistics on models
model_stats = ldply(calibrated_models,function(M){
	# M = calibrated_models[[1]]
	opt_pars = M$ode_parameters$parValues[M$ode_parameters$index_k]
	opt_stats = M$getStatistics()
	cbind(t(as.data.frame(opt_pars)),opt_stats)

},.id="model",.inform = T)

model_stats$cell_line = gsub("_.*","",model_stats$model)

# select k parameters, lambda_k influence only those

index_parameters = grep("_k_",colnames(model_stats))

# add statistics on k parameters
model_stats$sumPars = rowSums(model_stats[,index_parameters])

model_stats$numNonZero  = rowSums(model_stats[,index_parameters]>1e-2)

library(ggrepel)
# plot sum of parameters vs SSE model fit
ggplot(model_stats,aes(SSE, sumPars, col=cell_line)) +
	geom_point()+
	geom_line() +
	geom_text_repel(aes(label=lambda_k)) +
	theme_bw() + ggtitle('Regularisation: fitting cost vs parameters') +
	xlab("Sum of squared error (SSE)") + ylab("sum of k parameters")
dir.create("./figures/regularisation",showWarnings = F)
if(FALSE) ggsave("./figures/regularisation/SSE_vs_SSPars_pkn_v4_midas_v2.pdf")



# plot number of non-zero parameters (>0.01) vs SSE model fit
ggplot(model_stats,aes(SSE, numNonZero, col=cell_line)) +
	geom_point()+
	geom_line() +
	geom_text_repel(aes(label=lambda_k)) +
	theme_bw() + ggtitle('Regularisation: fitting cost vs parameters') +
	xlab("Sum of squared error (SSE)") + ylab("sum of non-zero k parameters")
dir.create("./figures/regularisation",showWarnings = F)
if(FALSE) ggsave("./figures/regularisation/SSE_vs_non_zero_k_pars_pkn_v4_midas_v2.pdf")


library(gridExtra)
# how much fit do we sacrafice?
gg1= ggplot(model_stats,aes(lambda_k, SSE, col=cell_line)) +
	geom_point()+
	geom_line() +
	geom_text_repel(aes(label=lambda_k)) +
	scale_x_log10() +
	theme_bw() + ggtitle('Regularisation influence on fitting cost') +
	xlab("regularisation parameter") + ylab("Sum of squared error (SSE)")


gg2 = ggplot(model_stats,aes(lambda_k, numNonZero, col=cell_line)) +
	geom_point()+
	geom_line() +
	geom_text_repel(aes(label=lambda_k)) +
	scale_x_log10() +
	theme_bw() + ggtitle('Regularisation influence on parameters') +
	xlab("regularisation parameter") + ylab("sum of non-zero k parameter")

grid.arrange(gg1,gg2,nrow = 2, ncol = 1)

if(FALSE) ggsave("./figures/regularisation/regularisation_on_SSE_and_nonzeroPars_v4_midas_v2.pdf")



# scale number of parameters and sum of squares to overlap
model_stats = model_stats %>% group_by(cell_line) %>% mutate(rel_SSE = (SSE/min(SSE))) %>% ungroup()
model_stats = model_stats %>% group_by(cell_line) %>% mutate(rel_numNonZero = (numNonZero/max(numNonZero))) %>% ungroup()

model_stats = as.data.frame(model_stats)


ggplot(model_stats,aes(rel_SSE, rel_numNonZero, col=cell_line)) + geom_point() + geom_line() + geom_text_repel(aes(label=lambda_k))


mean_model_stats = model_stats %>% group_by(lambda_k) %>% summarise(mean_rel_SSE = mean(rel_SSE),
												 sd_rel_SSE = sd(rel_SSE),
												 mean_rel_numNonZero = mean(rel_numNonZero),
												 sd_rel_numNonZero = sd(rel_numNonZero))

gg1= ggplot(mean_model_stats,aes(lambda_k, mean_rel_SSE)) +
	geom_errorbar(aes(ymin=mean_rel_SSE-sd_rel_SSE,ymax=mean_rel_SSE+sd_rel_SSE))+
	geom_line(linetype="dashed") +
	geom_point() +
	scale_x_log10() + coord_cartesian( ylim = c(0.9, 1.5))+
	theme_bw() + ggtitle('Regularisation influence on fitting cost') +
	xlab("regularisation parameter") + ylab("averaged relative sum of squared error")


gg2= ggplot(mean_model_stats,aes(lambda_k, mean_rel_numNonZero)) +
	geom_errorbar(aes(ymin=mean_rel_numNonZero-sd_rel_numNonZero,ymax=mean_rel_numNonZero+sd_rel_numNonZero))+
	geom_line(linetype="dashed") +
	geom_point() +
	scale_x_log10() +
	theme_bw() + ggtitle('Regularisation influence on fitting cost') +
	xlab("regularisation parameter") + ylab("averaged relative number of nonzero parameters")

grid.arrange(gg1,gg2,nrow = 2, ncol = 1)
if(FALSE) ggsave("./figures/regularisation/regularisation_on_averaged_SSE_and_nonZeroPars_v4_midas_v2.pdf")


### SUMMARY #####
# The 5 cell lines overall shows that regularisation matters aroun 1e-5 and has almost
# no effect below this value.
# Model fitting error gets pretty high at 1e-2 and above (more than 20%, which
# is huge considering that many measurements are flat)
#
# The optimal regularisation parameter is around 1e-4 - 1e-5. At this range we can see a clear effect
# on the number of nonzero parameters, but the fitting error bias is less than 10%



