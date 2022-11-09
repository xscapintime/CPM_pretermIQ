#!/bin/bash



#for id in `cat ../../lists/list_liyang_[1-3].txt`
for id in `ls /data/project/BIPP/laila/AP_8yo_fmri_task/derivatives/fmriprep/ | grep html | grep sub- | grep -v EAP |sed "s/sub-//g" | sed "s/.html//g"` # check all AP and BIPP
do
	
	path=/data/project/BIPP/laila/AP_8yo_fmri_task/derivatives/fmriprep/sub-${id}/func
	pathpar=/data/project/BIPP/laila/AP_8yo_fmri_task/derivatives/freesurfer/sub-${id}/parcellation
	

	if [ -f ${path}/sub-${id}_cleanPreproc_rsOutput_tempFiltered_inT1wSpace.nii.gz ] && [ ! -s ${pathpar}/sub-${id}_glasserROI_timeseries_sept29.txt ]; then

		
		echo $id extracting TS
		./extract_ROI_timeseries_2_batch.sh $id

	fi
done
