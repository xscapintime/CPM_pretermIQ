#!/bin/bash

#module load fsl ## somehow module: command not found
#module load afni


## subject
subj=$1
echo $subj

#### Script for regressing out motion confounds ####
# Seven motion confounds chosen: global_signal,trans_x,trans_y,trans_z,rot_x,rot_y,rot_z

## input and output dir
path=/data/project/BIPP/laila/AP_8yo_fmri_task/derivatives/fmriprep/sub-${subj}/func

echo $path


# This created a file with columns of only the 7 motion variables we are interested in:
awk -F'\t' -vcols=global_signal,trans_x,trans_y,trans_z,rot_x,rot_y,rot_z '(NR==1){n=split(cols,cs,",");for(c=1;c<=n;c++){for(i=1;i<=NF;i++)if($(i)==cs[c])ci[c]=i}}{for(i=1;i<=n;i++)printf "%s" FS,$(ci[i]);printf "\n"}' ${path}/sub-${subj}_task-Exrest_acq-TR1160_desc-confounds_regressors.tsv > ${path}/only_confounds_with_header.tsv



#Remove header row from the only_confounds_with_header.tsv file:
awk '{if(NR>1)print}' ${path}/only_confounds_with_header.tsv > ${path}/only_confounds_no_header.tsv
echo "created file"

#now run fsl_regfilt:
echo "fslregfilt ${subj}"

#rm sub-${subj}_confoundsRegressed_task-rest-preproc_bold.nii.gz

fsl_regfilt -i ${path}/sub-${subj}_task-Exrest_acq-TR1160_space-T1w_desc-preproc_bold.nii.gz -d ${path}/only_confounds_no_header.tsv -f "1, 2, 3, 4,5 ,6 ,7" -o ${path}/sub-${subj}_confoundsRegressed_task-rest-preproc_bold.nii.gz

#now, run the bandpass filter

echo "bandpass ${subj}"

3dBandpass -band 0.01 0.1 -input ${path}/sub-${subj}_confoundsRegressed_task-rest-preproc_bold.nii.gz -prefix ${path}/sub-${subj}_cleanPreproc_rsOutput_tempFiltered_inT1wSpace.nii.gz





