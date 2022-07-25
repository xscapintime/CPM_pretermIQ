for list in *liyang_*.txt
do
	n=`echo $list | cut -d _ -f 3 | sed "s/.txt//g"`
	nl=`cat $list | wc -l`
	qsub -v list=${list} -N Liyang_fmriprep_Parallel_${n} -t 1-${nl} STEP_4_fMRIprep_batch.sh
	# qsub STEP_4_fMRIprep_batch.sh -v list=${list} -N Liyang_fmriprep_Parallel_${n} -t 1-${nl} -hold_jid ${previous job name} ## if still not enough storage
done

