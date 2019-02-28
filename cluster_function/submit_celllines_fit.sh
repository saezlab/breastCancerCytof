#!/bin/bash
echo 'Running model fitting'


for filename in ./pkn_v2_midas_v2/inputs/*.RDS; do
	
	jobname=cellline_${filename##*/}
	
	bsub -J $jobname -M 8192 -W 05:00 -o ${jobname}_out.txt -e ${jobname}_err.txt Rscript fit_model_args.R $filename
    echo "$filename submitted"
    sleep 0.5
done

