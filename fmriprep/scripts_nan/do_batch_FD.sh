#!/bin/bash


#for id in `cat ../../lists/list_liyang_[1-3].txt`
for id in `ls /data/project/BIPP/laila/AP_8yo_fmri_task/derivatives/fmriprep/ | grep -v html | grep sub- | grep -v EAP |sed "s/sub-//g"` # check all AP and BIPP
do
        path=/data/project/BIPP/laila/AP_8yo_fmri_task/derivatives/fmriprep/sub-${id}/func

        #echo $path

	if [ -f $path/sub-${id}_task-Exrest_acq-TR1160_desc-confounds_regressors.tsv ] && [ ! -s $path/sub-${id}_framwise_displacement.tsv ]; then

		echo $id get FD
		./get_FD_batch.sh $id

	fi

done

