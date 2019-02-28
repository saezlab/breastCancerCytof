# plot the CNOlist files for each cell lines

library(multiCellNOpt)

#MIDAS_files = list.files("data/MIDAS_v2/",pattern = ".csv",full.names = T)
MIDAS_files = list.files("data/MIDAS_v3/",pattern = ".csv",full.names = T)


pdf("./figures/MIDAS_v3_all_cellLines.pdf",width = 20,height = 5)
for(i in seq_along(MIDAS_files)){
	cno = CNOlist$new(MIDAS_files[[i]])
	cno$plot()
	title(gsub(".csv","",basename(MIDAS_files[[i]]),fixed = T))
}
dev.off()
