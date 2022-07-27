for list in /data/project/BIPP/laila/AP_8yo_fmri_task/lists/*_liyang_*.txt
do
	bn=`basename $list`
	n=`echo $bn | cut -d _ -f 3 | sed "s/.txt//g"`
	nl=`cat $list | wc -l`

	if [ $n == 1 ]
	then
		echo fmriprep_Parallel_${n} submitted
		echo $bn $nl
		qsub -v list=${list} -N fmriprep_Parallel_${n} -t 1-${nl} STEP_4_fMRIprep_run_forLiyang_batch.sh
		sleep 0.5
	else
		prej=`expr $n - 1`
		echo fmriprep_Parallel_${n} will be submitted after fmriprep_Parallel_${prej} is finised
		echo $bn $nl
		qsub -v list=${list} -N Liyang_fmriprep_Parallel_${n} -t 1-${nl} -hold_jid fmriprep_Parallel_${prej} STEP_4_fMRIprep_run_forLiyang_batch.sh ## if still not enough storage
		sleep 0.5

	fi		
done
