for list in /data/project/BIPP/laila/AP_8yo_fmri_task/lists/*_liyang_*.txt
do
	bn=`basename $list`
	n=`echo $bn | cut -d _ -f 3 | sed "s/.txt//g"`
	nl=`cat $list | wc -l`
	echo $bn $nl

	qsub -v list=${list} -N fmriprep_Parallel_${n} -t 1-${nl} STEP_4_fMRIprep_run_forLiyang_batch.sh
	# qsub STEP_4_fMRIprep_batch.sh -v list=${list} -N Liyang_fmriprep_Parallel_${n} -t 1-${nl} -hold_jid ${previous job name} ## if still not enough storage
done
