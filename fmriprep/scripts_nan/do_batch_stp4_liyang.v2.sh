#!/bin/bash

#awk '{print NR  " " $s}' /data/project/BIPP/laila/AP_8yo_fmri_task/lists/list_liyang_new.txt > list_tmp
echo BIPP{038..051} | sed -e 's/  */\n/g' | awk '{print NR  " " $s}' > list_tmp ## 20/09/22 add 14 more

cat list_tmp | while read id
do

	n=`echo $id | cut -d " " -f 1`
	sub=`echo $id | cut -d " " -f 2`

	if [ $n == 1 ]
	then
		echo fmriprep.${n} for ${sub} submitted
		#echo $sub
		qsub -v subj=${sub} -N fmriprep.${n} STEP_4_fMRIprep_run_forLiyang_batch.v2.sh
		sleep 0.5
	else
		prej=`expr $n - 1`
		echo fmriprep.${n} for ${sub} submitted, will start after fmriprep.${prej} is finised
		#echo $sub
		qsub -v subj=${sub} -N fmriprep.${n} -hold_jid fmriprep.${prej} STEP_4_fMRIprep_run_forLiyang_batch.v2.sh ## if still not enough storage
		sleep 0.5

	fi		
done

rm list_tmp
