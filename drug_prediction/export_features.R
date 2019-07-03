# export features

# revision: 01.07.2019
# - unify the timepoints of the predictions
# - use the best parameters that was obtained by the MFu sampling after the initial optimisation

library(multiCellNOpt)
library(tidyverse)

# calibrated_model_files = list.files("./data/models/pkn_v2_midas_v2/outputs/","*.RDS",full.names = T)
calibrated_model_files = list.files("./data/models/pkn_v4_midas_v4/outputs/","*.RDS",full.names = T)
feature_folder = "./data/models/pkn_v4_midas_v4/features/"

dir.create(feature_folder,recursive = T)

calibrated_models = lapply(calibrated_model_files,function(f){
	m = readRDS(f)
	# update the class definition to the newest version
	#m2 = logicODEModel$new(oldModel = m)
}
)

names(calibrated_models) = gsub(".RDS","",basename(calibrated_model_files))


### Model features 1: parameters -----------------------------------------------
# we export the estimated model parameters (k and tau) for each cell-lines.


# first export the parameters found in the optimisation:
feature_table_1 <-  calibrated_models %>% map_dfr(function(M){
	pars = data.frame(as.list(M$ode_parameters$parValues[M$ode_parameters$index_opt_pars]))
}) %>% bind_cols(tibble(cell_line=names(calibrated_models)))


write_rds(feature_table_1, path = paste0(feature_folder,"/ft_1_raw_parameters_optimised.rds"))


### Model features 1.2: best parameters found during MFU_sampling---------------
# import sampling parameters
mfu_parameters <- read_rds('./data/models/pkn_v4_midas_v4/parameter_samplings_QC.RDS') %>% as_tibble()

# find the parameters with smallest fobj for each cell-line.
feature_table_1_2 = mfu_parameters %>%
	group_by(cell_line) %>% top_n(1,-fobj) %>%
	select(cell_line,everything(),-fobj,-rel_fobj,-sampling_id,-run_id,-starts_with("x0"))

write_rds(feature_table_1_2, path = paste0(feature_folder,"/ft_2_raw_parameters_mfu.rds"))





### Model feautures 2.2 : timecourse edge strenght -----------------------------
source("./drug_prediction/getReactionActivity.R")
base_time = c(0,5.5,7,9,13,17,23,30,40,60)


# update the calibrated models with the best parameters found
mfu_parameters <- read_rds('./data/models/pkn_v4_midas_v4/parameter_samplings_QC.RDS') %>% as_tibble()

# find the parameters with smallest fobj for each cell-line.
feature_table_1_2 = mfu_parameters %>%
	group_by(cell_line) %>% top_n(1,-fobj) %>%
	select(cell_line,everything(),-fobj,-rel_fobj,-sampling_id,-run_id,-starts_with("x0"))


timecourse_strength_to_feature = function(model){
	reacData = getReactionActivity_logicModel(logicODEmodel = model,timeSignals = base_time)
	return(reacData$timecourse)
}

feature_table_2_2 = ldply(calibrated_models,function(M){

	mfu_min_fobj <- mfu_parameters %>% filter(cell_line == M$name) %>% pull(fobj) %>% min()
	if (M$getStatistics()$fObj > mfu_min_fobj){

		best_vec = mfu_parameters %>% filter(cell_line == M$name) %>%  # find the cell-line
			arrange(fobj) %>% slice(1) %>% # sleect the row with smallest objective function value
			select(everything(),-cell_line, -fobj,-rel_fobj,-sampling_id,-run_id)

		best_vec_pars <- best_vec %>% select(-starts_with("x0"))
		best_vec_x0 <- best_vec %>% select(starts_with("x0"))

		stopifnot(all(M$ode_parameters$parNames[M$ode_parameters$index_opt_pars]==names(best_vec_pars)))

		M$ode_parameters$parValues[M$ode_parameters$index_opt_pars] = best_vec_pars %>% flatten_dbl()

		M$ode_parameters$x0Values[M$ode_parameters$index_opt_x0] = best_vec_x0 %>% flatten_dbl()

		M$updateInitialValueEstimates(ode_parameters = M$ode_parameters,
									  initialConditions = M$initialConditions)
	}

	features = timecourse_strength_to_feature(model = M)
},.id = "cell_line" ,.progress = progress_text())

feature_table_2_2 <- feature_table_2_2 %>% as_tibble() %>% spread(reacID,y) %>%
	filter(exp !="exp 1") %>% mutate(exp = as.character(exp)) %>%
	mutate(exp = case_when(
		exp == "exp 2" ~ "EGF",
		exp == "exp 3" ~ "iEGFR",
		exp == "exp 4" ~ "iMEK",
		exp == "exp 5" ~ "imTOR",
		exp == "exp 6" ~ "iPI3K",
		exp == "exp 7" ~ "iPKC")) %>% select(-modelName)

write_rds(feature_table_2_2, path = paste0(feature_folder,"/ft_3_edge_strength.rds"))




