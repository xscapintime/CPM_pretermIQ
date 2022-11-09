#!/bin/bash

## if you have a fixed subject list text file
#awk '{print NR  " " $s}' /data/project/BIPP/laila/AP_8yo_fmri_task/lists/list_liyang_new.txt > list_tmp

## if you have subjects with continuous ID  
#echo BIPP{044..051} | sed -e 's/  */\n/g' | awk '{print NR  " " $s}' > list_tmp

## if there's only few subject
echo BIPP043 | awk '{print NR  " " $s}' > list_tmp ## 10/10/22 

cat list_tmp | while read id
do

	n=`echo $id | cut -d " " -f 1`
	sub=`echo $id | cut -d " " -f 2`

	if [ $n == 1 ]
	then
		echo fmriprep.${n} for ${sub} submitted

		qsub -v subj=${sub} -N fmriprep.${n} STEP_4_fMRIprep_run_forLiyang_batch.sh
		sleep 0.5
	else
		prej=`expr $n - 1`
		echo fmriprep.${n} for ${sub} submitted, will start after fmriprep.${prej} is finised
		## if still not enough storage
		## job dependency on SGE system, the next job will start running after the one before it finished
		qsub -v subj=${sub} -N fmriprep.${n} -hold_jid fmriprep.${prej} STEP_4_fMRIprep_run_forLiyang_batch.sh
		sleep 0.5

	fi		
done

rm list_tmp
