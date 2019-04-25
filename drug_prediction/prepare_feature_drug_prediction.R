# prepare cluster computationL


# prepare features from raw parameters -----------------------------------------

feature_data_file  = "./data/models/pkn_v4_midas_v4/features/feature_table_1_raw_parameters.RData"
load(feature_data_file)
output_folder = "./data/models/pkn_v4_midas_v4/features_vs_drug/inputs"
dir.create(output_folder,recursive = T)

IC50_cellline_matrix_file = "./data/DATA_GDSC/GDSC_IC50_breast"
IC50_matrix = read.csv(IC50_cellline_matrix_file)

feature_data = feature_table_1[2:129]
measured_cell_lines = as.character(feature_table_1$cell_line)
rownames(feature_data) = measured_cell_lines

# get only measured cell lines from drugs and get only drugged cell lines from measurements
drug_cell_lines = colnames(IC50_matrix)[2:49]

common_cell_lines = intersect(measured_cell_lines, drug_cell_lines)
common_IC50_matrix = IC50_matrix[,c("X",common_cell_lines)]
common_features = feature_data[common_cell_lines,]

saveRDS(common_cell_lines,file = file.path(output_folder,"common_cell_lines.RDS"))
saveRDS(common_IC50_matrix,file = file.path(output_folder,"common_IC50_matrix.RDS"))
saveRDS(common_features,file = file.path(output_folder,"common_rawPar_features.RDS"))



# prepare features from median interaction strenght -----------------------------------------

feature_data_file  = "./data/models/pkn_v4_midas_v4/features/feature_table_2_edge_mean_activity.RData"
load(feature_data_file)
output_folder = "./data/models/pkn_v4_midas_v4/features_vs_drug/inputs"
#dir.create(output_folder,recursive = T)

IC50_cellline_matrix_file = "./data/DATA_GDSC/GDSC_IC50_breast"
IC50_matrix = read.csv(IC50_cellline_matrix_file)

feature_data = feature_table_2[2:85]
measured_cell_lines = as.character(feature_table_2$cell_line)
rownames(feature_data) = measured_cell_lines

# get only measured cell lines from drugs and get only drugged cell lines from measurements
drug_cell_lines = colnames(IC50_matrix)[2:49]

common_cell_lines = intersect(measured_cell_lines, drug_cell_lines)
common_IC50_matrix = IC50_matrix[,c("X",common_cell_lines)]
common_features = feature_data[common_cell_lines,]

#saveRDS(common_cell_lines,file = file.path(output_folder,"common_cell_lines.RDS"))
#saveRDS(common_IC50_matrix,file = file.path(output_folder,"common_IC50_matrix.RDS"))
saveRDS(common_features,file = file.path(output_folder,"common_interStr_features.RDS"))



# prepare features from timecourse interaction strength -----------------------------------------

feature_data_file  = "./data/models/pkn_v4_midas_v4/features/feature_table_2_time_course_activity.RData"
load(feature_data_file)
output_folder = "./data/models/pkn_v4_midas_v4/features_vs_drug/inputs"
#dir.create(output_folder,recursive = T)

IC50_cellline_matrix_file = "./data/DATA_GDSC/GDSC_IC50_breast"
IC50_matrix = read.csv(IC50_cellline_matrix_file)
# exp 1 — not included
# exp 2 —  EGF+SERUM
# exp 3 — iEGFR
# exp 4 — iMEK
# exp 5 — imTOR
# exp 6 — iPI3K
# exp 7 — iPKC
# The next line casuse error in the end, bc for 3 cell lines there is no time 13, but 14.
# sub_data = filter(feature_table_2_2_timecourse_edge_strength,exp=="exp 3",time %in% c(13.0))
# so use features from 13/14. table() command shows that either one or the orther exist, never both
sub_data = filter(feature_table_2_2_timecourse_edge_strength,exp=="exp 3",time %in% c(13.0,14.0) )
# table(sub_data2[,c("cell_line","time")])






feature_data =sub_data [,c(5:88)]
measured_cell_lines = as.character(sub_data$cell_line)
rownames(feature_data) = measured_cell_lines

# get only measured cell lines from drugs and get only drugged cell lines from measurements
drug_cell_lines = colnames(IC50_matrix)[2:49]

common_cell_lines = intersect(measured_cell_lines, drug_cell_lines)
common_IC50_matrix = IC50_matrix[,c("X",common_cell_lines)]
common_features = feature_data[common_cell_lines,]

#saveRDS(common_cell_lines,file = file.path(output_folder,"common_cell_lines.RDS"))
#saveRDS(common_IC50_matrix,file = file.path(output_folder,"common_IC50_matrix.RDS"))
saveRDS(common_features,file = file.path(output_folder,"common_interStr_iEGFR_13_features.RDS"))

