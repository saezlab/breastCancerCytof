# prepare data
# Author: A. Gabor.
# date: 14.09.2018
# import data from Maco:
# there are 2 datasets:
#   - Median_allsamples_nocontrols_withcellcount.rds  I got  a bit later, but
#		contains cell counts
#	- Arcsinh_Median_Normalized.rds: contains all experiments, but no cell counts.
#
# preprocessing steps:
# 1. import
# 2. visualisation of the raw data: ploting each cell lines with all signals
# 3. simple scaling the data: we scale the data over all cell lines, for each marker.
#		this way signaling might be comparable between cell-lines.
#		Instead of the whole range of the signals, we use the 99% interquatiles and set the
#		1% to 0 or 1.
# 4. combine time-course A and B (2 biological replica) by interpolating the values
# 		at the union of the time-points of the experiments. Then we take the mean of
#		the real measurements/interpolated measurements. -- This should be checked in some
#		cases, because the 2 time courses might be significantly different.
# 5. export to MIDAS files. We build CNOlist structures and export the data to MIDAS file.


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
# time courses by cell lines ----------------------------------------------
markers = colnames(celline_data)[5:38]
cell_lines = as.character(unique(celline_data$cell_line))

celline_data_melted = melt(celline_data,measure.vars = markers,
						   variable.name = "markers",value.name = "signal")


# plot the time course for individual  cell-line markers:
# I am checking if we see signals
if(regenerateFigures){
	pdf("figures/cellLine_timecourse_raw.pdf",width = 14,height = 7)
	for(selected_cell_line in cell_lines){
		#selected_cell_line = "MDAMB157" "MX1"
		gg_data = filter(celline_data_melted,cell_line==selected_cell_line)
		gg_data$grouping = paste(gg_data$time_course, gg_data$treatment,sep = "_")
		gg = ggplot(gg_data,aes(as.numeric(time),signal,col=treatment)) +
			geom_point(aes(shape=time_course)) + geom_line(aes(group=grouping)) +
			facet_wrap(~markers) + ggtitle(selected_cell_line) + theme_bw() +
			xlab("time")

		print(gg)
	}
	dev.off()
}

# time course all cell lines ---------------------------------------------------
# to see the range of markers

if(regenerateFigures){
	gg_data = celline_data_melted
	gg_data$grouping = paste(gg_data$cell_line,gg_data$time_course, gg_data$treatment,sep = "_")

	pdf("figures/all_cellLine_dist_timecourse.pdf",width = 14,height = 7)
	# density plot
	gg = ggplot(gg_data) +
		geom_violin(aes(x=treatment, y=signal,fill=treatment),alpha=.3) +
		facet_wrap(~markers) + ggtitle("signal distributions: all time, all cell lines ") + theme_bw()
	print(gg)

	# this shows, that some markers might have outliers, like IdU, p-STAT, p-AKT
	# time course
	gg = ggplot(gg_data,aes(as.numeric(time),signal,col=treatment)) +
		geom_point(aes(shape=time_course)) + geom_line(aes(group=grouping)) +
		facet_wrap(~markers) + ggtitle("time-course: all cell lines") + theme_bw() +
		xlab("time")
	print(gg)

	dev.off()
}
# scale data for CellNOptR -----------------------------------------------------
#
# version 1:
# take the range of the signals over all cell-lines and conditions
# scaling is based on the .99 interquantile range.

celline_data_melted_ext = ddply(celline_data_melted, .(markers), function(df){
	#df = filter(celline_data_melted,markers=="IdU")

	q99 = quantile(df$signal,probs = c(0.005,0.995))
	df$scaled_signal = (df$signal - q99[[1]])/(q99[[2]]-q99[[1]])
	df$scaled_signal[df$scaled_signal>1] = 1
	df$scaled_signal[df$scaled_signal<0] = 0
	return(df)
})

if(regenerateFigures){
	gg_data = celline_data_melted_ext
	gg_data$grouping = paste(gg_data$cell_line,gg_data$time_course, gg_data$treatment,sep = "_")

	pdf("figures/all_cellLine_dist_timecourse_scaled.pdf",width = 14,height = 7)
	# density plot
	gg = ggplot(gg_data) +
		geom_violin(aes(x=treatment, y=scaled_signal,fill=treatment),alpha=.3) +
		facet_wrap(~markers) + ggtitle("signal distributions: all time, all cell lines ") + theme_bw()
	print(gg)

	# this shows, that some markers might have outliers, like IdU, p-STAT, p-AKT
	# time course
	gg = ggplot(gg_data,aes(as.numeric(time),scaled_signal,col=treatment)) +
		geom_point(aes(shape=time_course)) + geom_line(aes(group=grouping)) +
		facet_wrap(~markers) + ggtitle("time-course: all cell lines") + theme_bw() +
		xlab("time")
	print(gg)

	dev.off()
}

# create CNO object -----------------------------------------------------------
# we setup a CNOlist object from CellNOptR for each cell-line and
# then we can export the data to MIDAS files

library(CellNOptR)
cell_lines = as.character(unique(celline_data_melted_ext$cell_line))

# do naming manipulation if neccessary !
reporters = as.character(unique(celline_data_melted_ext$markers))

treatments = c("EGF","iEGFR", "iMEK", "imTOR","iPI3K","iPKC")
library(progress)
pb = progress_bar$new(total=length(cell_lines))

for(iCell_line in 1:length(cell_lines)){
	pb$tick()
	#cno_data = filter(celline_data_melted_ext,cell_line == "T47D", treatment != "full")
	cno_data = filter(celline_data_melted_ext,cell_line == cell_lines[[iCell_line]], treatment != "full")
	cno_data[cno_data$time=="0short"] = "0"
	cno_data$time = as.numeric(as.character(cno_data$time))

	timepoints = sort(unique(cno_data$time))

	cno = list()
	cno$namesCues = c("EGF","SERUM", "EGFR", "MEK", "mTOR", "PI3K", "PKC")
	cno$namesStimuli = c("EGF","SERUM")
	cno$namesInhibitors = c("EGFR", "MEK", "mTOR", "PI3K", "PKC")
	cno$namesSignals = reporters
	cno$timeSignals = timepoints

	cno$valueCues = matrix(0,nrow = 6,ncol = 7)
	colnames(cno$valueCues) = cno$namesCues
	cno$valueCues[2:6,3:7] = diag(5)
	cno$valueCues[,c("EGF","SERUM")] = 1

	cno$valueInhibitors = cno$valueCues[,cno$namesInhibitors]

	cno$valueStimuli = cno$valueCues[,cno$namesStimuli,drop=FALSE]


	### time-course generation  ----------------------------------------------------
	# There are 2 biological replica A and B. Measured in different time-point.
	# Marco's initial idea is to use linear interpolation across time in A and B and
	# compute missing signal values for the union of measurement times.
	# then take the simple mean of the measured and interpolated values at each and
	# every time point.
	#
	# Option B: we could fit a polynomial on all the measured point and predict the
	# measured values. This would have some smoothing effect if the time course in
	# A and B are not very nice.

	req_time = timepoints
	cno_data_combined = ddply(cno_data,.(markers,treatment),function(df){

		#df = filter(cno_data,treatment=="EGF",markers=="IdU")
		#df = filter(cno_data,treatment=="iPI3K",markers=="p-S6")

		df = df[order(df$time),]
		timecourse_A = filter(df,time_course=="A")
		timecourse_B = filter(df,time_course=="B")

		y_interp_A = approx(timecourse_A$time,timecourse_A$scaled_signal,req_time, method = "linear",rule=2)
		y_interp_B = approx(timecourse_B$time,timecourse_B$scaled_signal,req_time, method = "linear",rule=2)

		# plot(timecourse_A$time,timecourse_A$scaled_signal, col="blue",type="b",ylim=c(0,1),xlim=c(0,80))
		# lines(timecourse_B$time,timecourse_B$scaled_signal, col="red" ,type="b")
		# points(y_interp_A$x,y_interp_A$y,col="blue")
		# points(y_interp_B$x,y_interp_B$y,col="red")

		meanSignal = (y_interp_A$y + y_interp_B$y)/2
		sdSignal = abs(y_interp_A$y - y_interp_B$y)/sqrt(2)
		dfnew = data.frame(time=req_time,combined_signal=meanSignal,sdSignal=sdSignal)

	})

	# quality control: this only works in the loop of T47D:
	if(plotABmergingIssues){
		# plot the distance of the measirements versus the mean: the distance should be
		# small comparing to the mean, otherwise it indicates that the two tme courses are
		# different
		cno_data_combined$relSD = cno_data_combined$sdSignal/cno_data_combined$combined_signal
		cno_data_combined[order(cno_data_combined$relSD,decreasing = T),]
		ggplot(cno_data_combined) + geom_point(aes(combined_signal,sdSignal,col=markers))

		ggplot(filter(celline_data_melted_ext,cell_line == "T47D", treatment == "iPI3K", markers=="p-S6"),aes(as.numeric(as.character(time)),signal,col=time_course)) +
			geom_point() + geom_line(aes(group=time_course))
		ggplot(filter(celline_data_melted_ext,cell_line == "T47D", treatment == "iPI3K", markers=="p-S6K"),aes(as.numeric(as.character(time)),signal,col=time_course)) +
			geom_point() + geom_line(aes(group=time_course))

		ggplot(filter(celline_data_melted_ext,cell_line == "T47D", treatment == "iPI3K", markers%in%c("p-S6", "p-S6K","p-p90RSK")),aes(as.numeric(as.character(time)),signal,col=time_course)) +
			geom_point() + geom_line(aes(group=time_course)) + facet_wrap(~markers)

		ggplot(filter(celline_data_melted_ext,cell_line == "T47D", treatment == "iMEK", markers=="p-S6"),aes(as.numeric(as.character(time)),signal,col=time_course)) +
			geom_point() + geom_line(aes(group=time_course))

		ggplot(filter(celline_data_melted_ext,cell_line == "T47D", treatment == "iEGFR", markers=="p-CREB"),aes(as.numeric(as.character(time)),signal,col=time_course)) +
			geom_point() + geom_line(aes(group=time_course))

		ggplot(filter(celline_data_melted_ext,cell_line == "T47D", treatment == "iPI3K", markers=="p-STAT1"),aes(as.numeric(as.character(time)),signal,col=time_course)) +
			geom_point() + geom_line(aes(group=time_course))
	}


	valueSignals = vector("list",length = length(timepoints))
	names(valueSignals) = timepoints

	for(time_ind in 1:length(timepoints)){
		#time_ind = 2
		act_time = timepoints[[time_ind]]

		cno_data_t = filter(cno_data_combined,time==act_time)
		formated_S = dcast(filter(cno_data_t,time==act_time),time+treatments~markers,value.var = "combined_signal",fill = NA)

		valueSignals[[time_ind]] = as.matrix(formated_S[,markers])
	}
	cno$valueSignals = valueSignals


	writeMIDAS(CNOlist(cno),filename = paste0('data/MIDAS/',cell_lines[[iCell_line]],".csv"))
	#plot(CNOlist(paste0('data/MIDAS/',cell_lines[[iCell_line]],".csv")))
}



## ============================= END of EXPORT =================================


## Control experiments ---------------------------------------------------------

# separate technical controls --------------------------------------------------
celline_data_raw = readRDS("./data/Arcsinh_Median_Normalized.rds")
celline_data_technicalControl = filter(celline_data_raw,time_course == "control")
celline_data_noTech = filter(celline_data_raw,time_course != "control")

# separate biological controls -------------------------------------------------
# when the associated_cell_line is not NA, then it is a control: the cell line
# is HCC70 in this case
unique(celline_data_noTech[!is.na(celline_data_noTech$associated_cell_line),]$cell_line)
celline_data_bioControl = filter(celline_data_noTech,!is.na(associated_cell_line))
celline_data = filter(celline_data_noTech,is.na(associated_cell_line))

markers = colnames(celline_data)[38:71]

# show summary of data ---------------------------------------------------------
print(summary(celline_data,maxsum = 67))


# remove singlets treatment (whats that?)
celline_data = filter(celline_data,treatment!="singlets")

# the following cell-lines have more data then others - why? - ask MARCO
exceptions_cell_lines = c("MDAMB157","MDAMB436","HCC2157","MCF10F")


celline_data_technicalControl = filter(celline_data_raw,time_course == "control")
# control cell Lines
control_CLs = unique(celline_data_technicalControl$cell_line)



### HCC70 as technical control -------------------------------------------------
# HCC70 is used for technical and biological controls
# when time_course == control  --> represents a technical control
# technical control:
#	- only EGF stimulation
#	- time 0 and 10
# 	- associated_cell_line: shows that with which cell-line the technical control was used with

filter(celline_data_technicalControl,cell_line=="HCC70",time_course=="control")

# confirm EGF stimulation via p-ERK activity in technical controls:
ggplot(
	filter(celline_data_technicalControl,cell_line=="HCC70",time_course=="control"),
	aes(time,`p-ERK`,col=associated_cell_line)) +
	geom_point() +
	geom_line(aes(group=associated_cell_line)) +
	ggtitle("p-ERK in technical controls (HCC70) associated to each cell lines (EGF stim.)")+
	theme_bw() + ylab("p-ERK signal")
# on each plate we see an increased p-ERK from T0 to T10, which confirms that EGF was probably used.


### MDAMB453 as technical control ----------------------------------------------
# MDAMB453 is used for technical
# when time_course == control  --> represents a technical control
# technical control:
#	- only EGF stimulation
#	- time 0, 7 and 15 mins
# 	- associated_cell_line: shows that with which cell-line the technical control was used with

filter(celline_data_technicalControl,cell_line=="MDAMB453",time_course=="control")

# confirm EGF stimulation via p-ERK activity in technical controls:
ggplot(
	filter(celline_data_technicalControl,cell_line=="MDAMB453",time_course=="control"),
	aes(as.numeric(time),`p-ERK`,col=associated_cell_line)) +
	geom_point() +
	geom_line(aes(group=associated_cell_line)) +
	ggtitle("p-ERK in technical controls (MDAMB453) associated to each cell lines (EGF stim.)")+
	theme_bw() + ylab("p-ERK signal")
# on each plate we see an increased p-ERK from T0 to T10, which confirms that EGF was probably used.
# !!! p-ERK measurement drops at T7 on TWO plates strongly and on many plates slightly, what is the interpretation of this ?


