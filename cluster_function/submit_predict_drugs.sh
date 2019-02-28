#!/bin/bash
echo 'Submitting drug predictions from features'

IC50_file="./inputs/common_IC50_matrix.RDS"
#feature_file="./inputs/common_rawPar_features.RDS"
feature_file="./inputs/common_interStr_features.RDS"

# for drug_ID in 0; do
for drug_ID in {1..265}; do

	jobname=interStr_$drug_ID

	bsub -J $jobname -P jrc_combine -M 8192 -W 01:00 -o ${jobname}_out.txt -e ${jobname}_err.txt Rscript predict_drug_args.R $IC50_file $feature_file $drug_ID
    echo "$jobname submitted"
    sleep 1
done

