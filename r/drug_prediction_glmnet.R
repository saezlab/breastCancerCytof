# testing features by linear models
#
# --- Getting IC50 values - cell line name issues -----------------------------
# Mi gave a file that contains a _subset_ of cell lines (data/GDSC_CytoF_iPI3K_13)
# then I tried to get the cell-line drug reponse from the GDSC panel
# database: data/v17.3_fitted_dose_response.xlsx
# the problem is cell-line names.
# Marco used cell line names and alternative cell line names, but sometimes, these
# are not in the GDSC panel.
# --> we should contact Marco about this naming issue, until then I use what i got from Mi.

library(glmnet)
drug_reponse_IC50 = read.csv("data/GDSC_CytoF_iPI3K_13")
cell_line_names = colnames(drug_reponse_IC50)[2:ncol(drug_reponse_IC50)]

load("features/all_measured_and_predicted_iPI3K_time13.RData")
load("features/feature_table_2.RData")
features = all_measured_and_predicted_iPI3K_time13
features = feature_table_2


rownames(features)= features$cell_line
features$cell_line = NULL
feature_matrix = t(features)

# order as it is in the IC50 data
feature_matrix = feature_matrix[,cell_line_names]
feature_matrix = t(feature_matrix)
#feature_matrix = feature_matrix[,!colnames(feature_matrix)%in% c("EGF","SERUM")]  # remove NA columns

### for each drug build a linear model based on the features:
# idrug = 100
idrug = which(drug_reponse_IC50$X=="TW 37") # all_measured_and_predicted_iPI3K_time13
idrug = which(drug_reponse_IC50$X=="GW 441756") # feature.table.2
idrug = which(drug_reponse_IC50$X=="Camptothecin") # feature.table.2

for(idrug in 1: nrow(drug_reponse_IC50)){
	IC50 = t(drug_reponse_IC50[idrug,-1])

	if(!all(rownames(IC50) == rownames(feature_matrix))) stop('feature matrix and IC50 cell line names are not passing to each other')

	remove_NA_cellline = which(!is.na(IC50))
	IC50_clean = IC50[remove_NA_cellline]
	feature_matrix_clean = feature_matrix[remove_NA_cellline,]

	m.cv = cv.glmnet(y=IC50_clean, x=feature_matrix_clean,standardize.response=TRUE)

	plot(m.cv)
	coef(m.cv$glmnet.fit)[,which(m.cv$glmnet.fit$lambda == m.cv$lambda.min)]
	plot(predict.cv.glmnet(m.cv,newx = feature_matrix_clean,s = "lambda.min"),IC50_clean)
	print(m)
	coef(m)
}



############ plot correlation between drug and feature ########################


cell_line_data = melt(drug_reponse_IC50,id.vars = "X",variable.name = "cell_line_name",value.name = "IC50")

colnames(cell_line_data)[[1]] = "Drug"

feature_matrix_df = melt(feature_matrix,varnames = c("cell_line_name","feature"),value.name = "feature_value")


M = merge(cell_line_data,feature_matrix_df)

feature_drug_corr = ddply(M,.(Drug,feature),function(x){
	#x = dplyr::filter(M,Drug=="Erlotinib",feature=="IdU")
	c = cor.test(x$feature_value,x$IC50,use = "complete.obs")
	data.frame(t.stat = c$statistic,
			   df = c$parameter,
			   p.value = c$p.value,
			   corr.estim = c$estimate,
			   CI.low.95 = c$conf.int[[1]],
			   CI.high.95 = c$conf.int[[2]])
},.progress = progress_text()
)
feature_drug_corr = feature_drug_corr[order(feature_drug_corr$corr.estim),]
feature_drug_corr = feature_drug_corr[order(feature_drug_corr$p.value),]

feature_drug_corr$signif = ifelse(feature_drug_corr$p.value<0.05,TRUE,FALSE)
library(ggrepel)
ggplot(dplyr::filter(feature_drug_corr,Drug=="GW 441756"),aes(corr.estim,-log10(p.value))) + geom_point() +
	geom_text_repel(aes(label=ifelse(signif,as.character(feature),"")),color="red")+
	xlab("corr (IC50 vs reaction strenght)") + theme_bw() + ggtitle("DRUG: GW 441756")



ggplot(dplyr::filter(M,Drug=="GW 441756",feature=="MAPKAPK2=CREB")) + geom_point(aes(feature_value,IC50))
ggplot(dplyr::filter(M,Drug=="GW 441756")) + geom_point(aes(feature_value,IC50)) + facet_wrap(~feature)
ggplot(dplyr::filter(M,Drug=="Camptothecin")) + geom_point(aes(feature_value,IC50)) + facet_wrap(~feature)

