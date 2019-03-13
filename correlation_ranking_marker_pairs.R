# plot correlation between signals


library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)

regenerateFigures = FALSE
plotABmergingIssues = FALSE

celline_data_raw = readRDS("./data/Median_allsamples_nocontrols_withcellcount.rds")
# merge in basically the cellcount only

colnames(celline_data_raw)[which(colnames(celline_data_raw)=="dmt$cellcount")] ="cell_count"

head(celline_data_raw)

celline_data = celline_data_raw

# rename markers ---------------------------------------------------------------
markers = colnames(celline_data)[5:38]

markers[which(markers=="p-Akt(Ser473)")] = "p-AKT_S473"
markers[which(markers=="p-MKK3-MKK6")] = "p-MKK36"
markers[which(markers=="p-AKT(Thr308)")] = "p-AKT_T308"
markers[which(markers=="p-S6K")] = "p-p70S6K"
markers[which(markers=="p-ERK")] = "p-ERK12"
markers[which(markers=="p-MEK")] = "p-MEK12_S221"
markers[which(markers=="p-GSK3b")] = "p-GSK3B"

celline_data$treatment = as.character(celline_data$treatment)
celline_data$treatment = gsub("iMEK","iMEK12",celline_data$treatment,fixed = T)


colnames(celline_data)[5:38] = markers
# time courses by cell lines ----------------------------------------------
markers = colnames(celline_data)[5:38]
cell_lines = as.character(unique(celline_data$cell_line))

gg_markers = gsub("-","_",markers)
colnames(celline_data)[5:38] = gg_markers

i=1
j=2

for(i  in 1:(length(gg_markers)-1) ){
	for(j in (i+1):length(gg_markers)){
		pdf(paste0("./figures/correlationPlots/marker_correlation_",gg_markers[[i]],"_",gg_markers[[j]],".pdf"),width = 13,height = 10)

		gg = ggplot(celline_data,aes_string(x=gg_markers[[i]],y=gg_markers[[j]],col="treatment")) +geom_line(alpha=.2) +facet_wrap(~cell_line,scales = "free")+ geom_point() + theme_bw()
		print(gg)
		dev.off()
	}
}



# correlation analysis:

celline_data_clean = filter(celline_data,treatment!= "full")

celline_cor_cov = ddply(celline_data_clean,.(cell_line),function(df){
	#df = filter(celline_data_clean,cell_line=="MX1")
	cor_res = cor(df[,gg_markers])
	cov_res = cov(df[,gg_markers])

	# remove upper triangle of the matrix:
	# cor_res_long = subset(as.data.frame(as.table(cor_res)),
	# 	   match(Var1,gg_markers) > match(Var2, gg_markers))
	# cov_res_long = subset(as.data.frame(as.table(cov_res)),
	# 					  match(Var1,gg_markers) > match(Var2, gg_markers))

	# DONT  remove upper triangle of the matrix:
	cor_res_long = as.data.frame(as.table(cor_res))
	cov_res_long = as.data.frame(as.table(cov_res))


	cor_res_long = cor_res_long[order(cor_res_long$Freq,decreasing = T),]
	cov_res_long = cov_res_long[order(cor_res_long$Freq,decreasing = T),]
	colnames(cor_res_long)[[3]] = "corr"
	colnames(cov_res_long)[[3]] = "cov"
	merge(cor_res_long,cov_res_long)

})



# add ranking: in which percent of the top correlation it is
celline_cor_cov$corr_percent_rank = percent_rank(celline_cor_cov$corr)
celline_cor_cov$cov_percent_rank = percent_rank(celline_cor_cov$cov)

celline_cor_cov = celline_cor_cov[order(celline_cor_cov$corr_percent_rank,decreasing = T),]
celline_cor_cov = celline_cor_cov[order(celline_cor_cov$cov_percent_rank,decreasing = T),]
head(celline_cor_cov)

write.csv(celline_cor_cov,file = "data/celline_cor_cov.csv")

# compute the ranking for each cell line then avera over cell lines:
# compiute  ranking of pairs for each cell line
celline_cor_cov_ranks = ddply(celline_cor_cov,.(cell_line),function(df){
	#df = filter(celline_cor_cov,cell_line=="MX1")
	df$corr_percent_rank = percent_rank(df$corr)
	df$cov_percent_rank = percent_rank(df$cov)
	df

	})
# avera the rank:
mean_ranked_cor_cov = ddply(celline_cor_cov_ranks,.(Var1,Var2),function(df){
	#df = filter(celline_cor_cov,Var1=="p_H3",Var2=="p_S6")
	colMeans(df[,4:7])
})

mean_ranked_cor_cov = mean_ranked_cor_cov[order(mean_ranked_cor_cov$cov_percent_rank,decreasing = T),]

corrplot::corrplot(cor_res)
hist(cor_res_long$Freq)
plot(df$p_p90RSK,df$p_ERK12)
plot(df$p_MKK4,df$p_H3)
p_MKK4         p_H3

