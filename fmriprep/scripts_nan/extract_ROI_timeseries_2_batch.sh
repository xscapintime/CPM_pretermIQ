#!/bin/bash

#module load ants
#module load afni
#module load fsl

#for subj in `cat /data/project/BIPP/laila/AP_8yo_fmri_task/lists/list_BIPP_resting.txt` 
#do

##### ==== 1 ==== First, we are going to reduce the dimensionality of the 4-d image to 3-d by getting the mean of the timeseries
subj=$1
echo "getting mean of ${subj} func image"

#cd /data/project/BIPP/laila/AP_8yo_fmri_task/derivatives/fmriprep/sub-${subj}/func/ 

path=/data/project/BIPP/laila/AP_8yo_fmri_task/derivatives/fmriprep/sub-${subj}/func

fslmaths ${path}/sub-${subj}_cleanPreproc_rsOutput_tempFiltered_inT1wSpace.nii.gz -Tmean ${path}/sub-${subj}_cleanPreproc_rsOutput_tempFiltered_inT1wSpace_MEAN.nii.gz

##### ==== 2 ==== We are going to achieve two things in one step:
# 1. put the atlas in T1-weighted space 
# 2. and downsample the atlas (to give the resting state image and the atlas the same dimensions and resolution) 
# reference image: the mean (in 3d) of the 4d resting state fully pre-processed image in T1-w space, input image: parcellated atlas from rafa script

#cd /data/project/BIPP/laila/AP_8yo_fmri_task/derivatives/freesurfer/sub-${subj}/parcellation/

pathpar=/data/project/BIPP/laila/AP_8yo_fmri_task/derivatives/freesurfer/sub-${subj}/parcellation
echo "downsampling the atlas and putting in T1-weighted space for ${subj}"

antsApplyTransforms -d 3 -r ${path}/sub-${subj}_cleanPreproc_rsOutput_tempFiltered_inT1wSpace_MEAN.nii.gz -i ${pathpar}/HCPMMP1_seq.nii.gz -o ${pathpar}/HCPMMP1_seq_inT1w_and_downsampled.nii.gz -t [/data/project/BIPP/laila/AP_8yo_fmri_task/derivatives/fmriprep/sub-${subj}/anat/sub-${subj}_from-fsnative_to-T1w_mode-image_xfm.txt, 0] -n NearestNeighbor


###### ===== 3 ==== now, we can extrtact the time series using the Glasser atlas ROIs for each subject
echo "extracting ROI time series for ${subj}"

rm ${pathpar}/sub-${subj}_glasserROI_timeseries_sept29.txt

3dROIstats -mask ${pathpar}/HCPMMP1_seq_inT1w_and_downsampled.nii.gz ${path}/sub-${subj}_cleanPreproc_rsOutput_tempFiltered_inT1wSpace.nii.gz > ${pathpar}/sub-${subj}_glasserROI_timeseries_sept29.txt


#done




