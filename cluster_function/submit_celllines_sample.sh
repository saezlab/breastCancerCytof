#!/bin/bash
echo 'Running model fitting'


for filename in ./pkn_v4_midas_v4/outputs/*.RDS; do

	jobname=cellline_${filename##*/}
	run_id=1

	bsub -J $jobname -M 8192 -W 01:00 -o ${jobname}_out.txt -e ${jobname}_err.txt Rscript sample_parameter_space_args.R $filename $run_id
    echo "$filename submitted"
    sleep 1
done

