#!/bin/bash
# loop for list of subjects

# cat ../../lists/list_liyang_[1-3].txt > list_tmp

#module load cuda # not working for bash, load them in command line
#module load fsl  # cannot revoke another shell from scripts
#module load afni


#echo "BIPP026" > list_tmp
for id in `echo BIPP026` #`cat ../../lists/list_liyang_[1-3].txt`
do
	path=/data/project/BIPP/laila/AP_8yo_fmri_task/derivatives/fmriprep/sub-${id}/func

	#echo $path
	
	if  [[ ! -e ${path}/sub-${id}_cleanPreproc_rsOutput_tempFiltered_inT1wSpace.nii.gz  &&  -e ${path}/sub-${id}_task-Exrest_acq-TR1160_desc-confounds_regressors.tsv ]]; then

		echo creating confounds file for ${id}

		./fsl_regfilt_script_2_batch.sh $id

		#qsub -v subj=$id -N fsl_regfilt_${id} fsl_regfilt_script_2_batch.sh
	fi
done

# rm list_tmp #somehomw can not found rm....
