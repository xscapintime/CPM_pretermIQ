#!/bin/bash


# loop for list of subjects
#for subj in `cat /data/project/BIPP/laila/AP_8yo_fmri_task/lists/list_BIPP_resting.txt`
#do

#echo "creating FD file for ${subj}"

#cd /data/project/BIPP/laila/AP_8yo_fmri_task/derivatives/fmriprep/sub-${subj}/func/


## subject
subj=$1

## input and output dir
path=/data/project/BIPP/laila/AP_8yo_fmri_task/derivatives/fmriprep/sub-${subj}/func



# This created a file with column name of FD:
awk -F'\t' -vcols=framewise_displacement '(NR==1){n=split(cols,cs,",");for(c=1;c<=n;c++){for(i=1;i<=NF;i++)if($(i)==cs[c])ci[c]=i}}{for(i=1;i<=n;i++)printf "%s" FS,$(ci[i]);printf "\n"}' ${path}/sub-${subj}_task-Exrest_acq-TR1160_desc-confounds_regressors.tsv > ${path}/sub-${subj}only_FD_with_header.tsv 


#Remove header row from the only_FD_with_header.tsv file and the first line (that's why NR>2):
awk '{if(NR>2)print}' ${path}/sub-${subj}only_FD_with_header.tsv > ${path}/sub-${subj}_framwise_displacement.tsv

#done



