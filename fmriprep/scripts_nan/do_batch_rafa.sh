#!/bin/bash

#for id in `cat ../../lists/list_liyang_[1-3].txt`
for id in `ls /data/project/BIPP/laila/AP_8yo_fmri_task/derivatives/fmriprep/ | grep -v html | grep sub- | grep -v EAP |sed "s/sub-//g"` # check all AP and BIPP

do
        path=/data/project/BIPP/laila/AP_8yo_fmri_task/derivatives/freesurfer/sub-${id}/

        #echo $path

        if  [ ! -f ${path}/parcellation/*_seq.nii.gz ]; then

		echo figshare anno sub-${id}
		#echo $path
		
                ./rafa_script_with_figshare_annot_batch.sh sub-${id} /data/project/BIPP/laila/AP_8yo_fmri_task/derivatives/freesurfer 

                #qsub -v subj=$id -N fsl_regfilt_${id} fsl_regfilt_script_2_batch.sh
        fi
done
