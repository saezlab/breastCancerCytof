# load elastic net predictons, compute test statistics and generate figures.

library(stringr)
library(ggplot2)
library(ggrepel)
library(plyr)
library(dplyr)
library(reshape2)


### User inputs ------

featureStr = "rawPars"
featureStr = "interStr"
featureStr = "interStr_iEGFR_13"

features = readRDS("./data/models/pkn_v4_midas_v4/features_vs_drug/inputs/common_rawPar_features.RDS")
features = readRDS("./data/models/pkn_v4_midas_v4/features_vs_drug/inputs/common_interStr_features.RDS")
features = readRDS("./data/models/pkn_v4_midas_v4/features_vs_drug/inputs/common_interStr_iEGFR_13_features.RDS")

figOutFolder = "./figures/drug_prediction/rawPars"
figOutFolder = "./figures/drug_prediction/interStr"
figOutFolder = "./figures/drug_prediction/interStr_iEGFR_13"

drug_files = list.files("./data/models/pkn_v4_midas_v4/features_vs_drug/outputs/",".*rawPar",full.names = T)
drug_files = list.files("./data/models/pkn_v4_midas_v4/features_vs_drug/outputs/",".*interStr_features",full.names = T)
drug_files = list.files("./data/models/pkn_v4_midas_v4/features_vs_drug/outputs/",".*iEGFR",full.names = T)

IC50 = readRDS("./data/models/pkn_v4_midas_v4/features_vs_drug/inputs/common_IC50_matrix.RDS")

plotFigs = T

#### ------
dir.create(figOutFolder)

show_features_on_heatmap = function(features,fileName=NULL){
	##
	if(!is.null(fileName)){
		pdf(fileName,width = 7.3,height = 4.5)}
	pheatmap::pheatmap(features,cluster_cols = F,
					   cluster_rows = F,colorRampPalette(RColorBrewer::brewer.pal(n = 7, name ="Reds"))(100),
					   fontsize_col = 4, fontsize_row = 4)
	if(!is.null(fileName)){	dev.off()}

}


if(plotFigs){
	show_features_on_heatmap(features)
	show_features_on_heatmap(features,fileName =file.path(figOutFolder,paste0(featureStr,"_feature_heatmap.pdf")))
}


# for hypothesis testing :
# bonferroni corrected p-value
p_limit = 0.05/nrow(IC50)

feature_names = colnames(features)

# f= drug_files[[2]]
# load the elastic net model results:
# drug_assoc:  df, columns: coeff vector for each feature; RMSE, RMSE_rm, drug_id for each model (1 cell line left out in each)
drug_assoc = ldply(drug_files,function(f){
	drug_id = stringr::str_extract(basename(f),"[0-9]+")
	df = readRDS(f)
	df$drug_id = drug_id
	return(df)
})

drug_assoc$drug_id = as.numeric(drug_assoc$drug_id)
drug_assoc = drug_assoc[order(drug_assoc$drug_id),]
colnames(drug_assoc)[grep("V[0-9]+",colnames(drug_assoc))] = c("intercept", feature_names)
drug_assoc$drug_name = as.character(IC50$X)[drug_assoc$drug_id]

# drug_assoc_test: data.frame, checking significance between RMSE and RMSE_rn by 3 test
drug_assoc_test = ddply(drug_assoc,.(drug_id),function(res){
	# res= filter(drug_assoc,drug_id=="4")

	RMSE = res$RMSE
	RMSE_rm = res$RMSE_rm

	t.res<-t.test(RMSE, RMSE_rm, paired=T, alternative = "less")
	t2.res<-t.test(RMSE, RMSE_rm, paired=F, alternative = "less")
	w.res<-wilcox.test(RMSE, RMSE_rm, paired=T, alternative = "less")

	#print(t.res)
	#print(w.res)
	p.vals = data.frame(t.paired=t.res$p.value,
						t.unpaired=t2.res$p.value,
						wilcox=w.res$p.value,
						RMSE_mean = mean(RMSE),
						t.paired.mean.diff = t.res$estimate)
})

drug_assoc_test$drug_name =  as.character(IC50$X)[drug_assoc_test$drug_id]
drug_assoc_test = drug_assoc_test[order(drug_assoc_test$wilcox),]

if(FALSE) write.csv(drug_assoc_test,file = paste0("./data/models/pkn_v4_midas_v4/features_vs_drug/predictable_drugs_by_",featureStr,".csv"),quote = F,row.names = F)


### PLOT 1: p-values for different statistics between RMSE and RMSE_rn: model fittid to random perturbation of IC50.
if(plotFigs){
	ggdata = cbind(drug_assoc_test,drug_name=IC50$X)
	ggdata = melt(ggdata,id.vars = c("drug_name","drug_id","RMSE_mean","t.paired.mean.diff"),variable.name = "test",value.name = "p.value")
	ggdata$label = as.character(ggdata$drug_name)
	ggdata$label[ggdata$p.value>0.05/265] = ""


	ggplot(filter(ggdata,test=="t.paired"),aes(t.paired.mean.diff,-log10(p.value),col=nchar(label)>1)) +
		geom_point()+geom_text_repel(aes(label=label)) +
		theme_bw() +
		ggtitle("Drug prediction, comparison to random models (2-sample, paired T-test)") +
		xlab("RMSE difference") +
		scale_color_manual(values = c("TRUE"="red","FALSE"="gray"),name="signif.")

	ggsave(file.path(figOutFolder,paste0(featureStr,"_feature_RMSE_random_vs_original_drug_vulcano.pdf")),width = 9,height = 7)
}

if(plotFigs){
	ggplot(ggdata,aes(test,p.value,col=RMSE_mean)) + geom_point()+
		scale_y_log10() +
		#scale_color_gradient(low="red",high="gray",name="RMSE")+
		scale_color_gradientn(colors=heat.colors(7, alpha = 1),name="RMSE")+
		geom_text_repel(aes(label=label)) +
		ggtitle("Drug prediction, comparison to random models") +
		xlab("2-sample tests") +
		ylab("p-value") +
		theme_bw()

	ggsave(file.path(figOutFolder,paste0(featureStr,"_feature_RMSE_random_vs_original_3test_pvalues.pdf")),width = 9,height = 7)
}

## Choose the paired Wilcoxon rank-sum test

# strict: all test:
# signif_rows = which(drug_assoc_test$t.paired<p_limit &drug_assoc_test$t.unpaired<p_limit & drug_assoc_test$wilcox<p_limit)
signif_rows = which(drug_assoc_test$wilcox<p_limit)
signif_drug_ids = drug_assoc_test$drug_id[signif_rows]
signif_drug_assoc = drug_assoc[drug_assoc$drug_id %in% signif_drug_ids,]


if(plotFigs){
	ggdata_signif_drug_assoc_violin = signif_drug_assoc
	ggdata_signif_drug_assoc_violin = melt(ggdata_signif_drug_assoc_violin,id.vars = c("drug_id","drug_name"),measure.vars = c("RMSE","RMSE_rm"),variable.name = "data_type",value.name = "RMSE")
	ggdata_signif_drug_assoc_violin$data_type = as.character(ggdata_signif_drug_assoc_violin$data_type)
	ggdata_signif_drug_assoc_violin$data_type[ggdata_signif_drug_assoc_violin$data_type=="RMSE"] = "original"
	ggdata_signif_drug_assoc_violin$data_type[ggdata_signif_drug_assoc_violin$data_type=="RMSE_rm"] = "random"

	ggplot(ggdata_signif_drug_assoc_violin, aes(x=data_type,y=RMSE,col=data_type)) +
		geom_violin() +
		geom_dotplot(aes(fill=data_type),binaxis='y', stackdir='center', dotsize=1) +
		scale_color_discrete(name="data",aesthetics = c("colour","fill")) +
		#geom_jitter(alpha=.5,height = 0	) +
		facet_wrap(~drug_name,scales = "free_y") +
		ggtitle("compare IC50 predictions based on random/real data ") +
		theme_bw() + xlab("")

	ggsave(file.path(figOutFolder,paste0(featureStr,"_feature_RMSE_random_vs_original_violin.pdf")),width = 12, height = 9)
}

### Parameter tests ************
# for the models with significant improvement, which are the predicting parameters
library(gridExtra)
if(plotFigs){ pdf(file.path(figOutFolder,paste0(featureStr,"_feature_importance.pdf")),width = 10,height = 9)}
signif_features_stats = ddply(signif_drug_assoc,.(drug_name),function(df,plotFlag){
	# df = filter(signif_drug_assoc,drug_name == "PIK-93")

	non_zero_features = which(colSums(as.matrix(df[,feature_names]))!=0)
	drug_name = df$drug_name[[1]]


	dfm=melt(df,measure.vars = feature_names,variable.name = "feature",value.name = "value")

	# calculate how many times the feature's coeff is non-zero
	NNZ_features = summarise_at(df,feature_names,function(x)sum(abs(x)>0))
	NNZ_df = data.frame(feature = colnames(NNZ_features),NNZ = as.numeric(NNZ_features[1,]))
	NNZ_df = NNZ_df[order(NNZ_df$NNZ,decreasing = T),]
	NNZ_df$feature = factor(NNZ_df$feature,levels = NNZ_df$feature)

	# set the order of the features according to mean values
	dfm_feature_sorted = ddply(dfm,.(feature),function(f){mean(f$value)})
	dfm_feature_sorted  = dfm_feature_sorted[order(dfm_feature_sorted$V1),]
	dfm = dfm[order(dfm$feature),]
	dfm$feature = factor(dfm$feature,levels=unique(dfm_feature_sorted$feature))



	feature_coeff_stats = ddply(dfm,.(feature),summarise,mean_value=mean(value),median_value=median(value),mad_value=mad(value))

	# plot results
	if(plotFlag){
		gg1 = ggplot(filter(NNZ_df,NNZ>0),aes(feature,NNZ)) + geom_col() +
			theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
			ggtitle(paste0(drug_name,", feature occurance")) + ylab('occurance as predictor')

		#, ! feature %in% c("EGFR=STAT5","EGFR=STAT1","EGFR=PKC")
		gg2 = ggplot(filter(dfm,feature %in% names(non_zero_features)),aes(feature,value))+
			geom_boxplot() +
			theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))+
			ggtitle(paste0(drug_name, ", feature coefficients")) + ylab('coeff. value')
		print(gridExtra::grid.arrange(gg1,gg2,nrow=2))
	}

	return(merge(NNZ_df,feature_coeff_stats))

},plotFlag=plotFigs)
if(plotFigs) {dev.off()}

# check the targets of the significant drugs

drug_descriptions = read.csv("./data/DATA_GDSC/DRUG_ANALYSIS_SET_update_improved.csv")

drug_targets= signif_drug_assoc[!duplicated(signif_drug_assoc$drug_name),c("drug_id","drug_name")]
drug_targets$drug_name_alternative = drug_targets$drug_name

# checking/checking some drugname by hand ...
drug_targets$drug_name[which(!drug_targets$drug_name %in% drug_descriptions$Drug.Name)]
drug_targets$drug_name_alternative[drug_targets$drug_name=="Mitomycin C"] = "Mitomycin-C"
drug_targets$drug_name_alternative[drug_targets$drug_name=="Zibotentan, ZD4054"] = "Zibotentan"
drug_targets$drug_name_alternative[drug_targets$drug_name=="GDC0941 (rescreen)"] = "Pictilisib"
drug_targets$drug_name_alternative[drug_targets$drug_name=="BIRB 0796"] = "Doramapimod"
drug_targets$drug_name_alternative[drug_targets$drug_name=="piperlongumine"] = "Piperlongumine"
drug_targets$drug_name_alternative[drug_targets$drug_name=="PLX4720 (rescreen)"] = "PLX-4720"
drug_targets$drug_name_alternative[drug_targets$drug_name=="AZD-2281"] = "Olaparib"
drug_targets$drug_name_alternative[drug_targets$drug_name=="Bleomycin (50 uM)"] = "Bleomycin"

# remove drugs for which target is not reported (why not found?)
drug_targets = drug_targets[which(! drug_targets$drug_name_alternative=="XL-880"),]

# irow = 14
indx_match = c()
for(irow in 1:nrow(drug_targets)){
	# check first the name (exactly, it is important!)
	idx = grep(paste0("^",drug_targets[irow,"drug_name_alternative"],"$"),drug_descriptions$Drug.Name)

	# if name is not matching, check Synonymes
	if(length(idx)==0) idx = grep(drug_targets[irow,"drug_name_alternative"],as.character(drug_descriptions$Synonyms))

	if(length(idx)==0)  {stop(paste(drug_targets[irow,"drug_name_alternative"]," not found"))}

	if(length(idx)>1) {
		# print(drug_descriptions[idx,])
		# stop if target/pathway is not the same, i.e. unique returns longer than 1.
		stopifnot(nrow(unique(drug_descriptions[idx,c("Target" , "Target.Pathway")]))==1)
		idx = idx[1]
	}

	indx_match = c(indx_match,idx)

}

drug_targets = cbind(drug_targets,drug_descriptions[indx_match,c("Target","Target.Pathway")])



if(FALSE) saveRDS(drug_targets,file=paste0("./data/models/pkn_v4_midas_v4/features_vs_drug/signif_drugs_by_",featureStr,".RDS"))

#####   Check individual drugs by their targets and features
if(featureStr=="interStr_iEGFR_13"){
	feature_drug_all = merge(signif_features_stats,drug_targets)
	feature_drug_all = merge(feature_drug_all,drug_assoc_test)

	filter(feature_drug_all,drug_name == "SB-715992")
	filter(feature_drug_all,drug_name == "QL-XII-47")
	filter(feature_drug_all,drug_name == "Genentech Cpd 10")
	filter(feature_drug_all,drug_name == "YM201636")
	filter(feature_drug_all,drug_name == "TG101348")

	# sort features by mean
	temp = ddply(feature_drug_all,.(feature),function(f){mean(f$median_value)})
	temp  = temp[order(temp$V1),]

	feature_drug_all$feature = factor(feature_drug_all$feature,levels=unique(temp$feature))
	ggplot(filter(feature_drug_all,abs(median_value)<1e6),aes(feature,median_value)) + geom_boxplot() +
		theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))+
		ggtitle("overall effect of edge features on IC50 (by coeff.)") + ylab('coeff. value')

	ggplot(filter(feature_drug_all,abs(median_value)<1e6),aes(feature,as.numeric(median_value))) + geom_violin(draw_quantiles = .5, scale = "width") +
		theme_bw() +theme(axis.text.x = element_text(angle = 90, hjust = 1))+
		ggtitle("overall effect of edge features on IC50 (by coeff.)") + ylab('coeff. value')


	M = readRDS("./data/models/pkn_v4_midas_v4/outputs/184A1.RDS")
	source("R/plotAnnotatedModel.R")

	# idrug = 6
	binariseCols = TRUE
	pdf(file.path(figOutFolder,paste0(featureStr,"_feature_importance_on_graph_binary.pdf")),width = 10, height = 10)
	for(idrug in 1:16){

		drug = drug_assoc_test$drug_name[[idrug]]

		coeffs = filter(feature_drug_all,drug_name == drug)


		edge_weight = coeffs$mean_value
		names(edge_weight) = coeffs$feature

		# if edge does not show up in 50 % then remove
		edge_weight[coeffs$NNZ < max(coeffs$NNZ)/2] = 0

		# remove edges where feature is anyway 0
		edge_weight[abs(edge_weight)>1e6] = 0


		# convert to binary
		if(binariseCols){
		edge_weight[edge_weight>0] = 1
		edge_weight[edge_weight<0] = -1
		}

		pal = colorRampPalette(rev(RColorBrewer::brewer.pal(7,"RdBu")))(101)

		edge_col = map2color(edge_weight,pal,limits=c(-max(abs(edge_weight)),max(abs(edge_weight))))

		reacID_colors = edge_col[match(M$pkn$network$reacID, names(edge_col))]

		plotAnnotatedModel(model = M$pkn$network,CNOlist = M$exps$convertToS4(),annotEdgeColors = reacID_colors,annotEdgeMatching = "reacID")

		target = as.character(drug_targets[drug_targets$drug_name==drug,c("Target")])
		target_pw = as.character(drug_targets[drug_targets$drug_name==drug,c("Target.Pathway")])
		title(main = paste0(drug," targets ",target," (",target_pw,")"))
	}

	dev.off()


}
