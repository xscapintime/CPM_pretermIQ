#! /bin/csh
#$ -S /bin/csh
#$ -cwd
#$ -pe smp 1
#$ -j y
#$ -tc 528
#$ -l h_vmem=32G
#$ -V


set NSLOTS=1
setenv MKL_NUM_THREADS $NSLOTS
setenv OMP_NUM_THREADS $NSLOTS
setenv NUMEXPR_NUM_THREADS $NSLOTS

module load fsl
module load freesurfer
module load fmriprep/20.1.1 
module list

#set list = /data/project/BIPP/laila/AP_8yo_fmri_task/lists/list_liyang.txt

set subj=`sed "${SGE_TASK_ID}q;d" $list`

echo ${subj}

fmriprep /data/project/BIPP/laila/AP_8yo_fmri_task/Nifti/ /data/project/BIPP/liyang_fmriprep/AP_8yo_fmri_task/output/ participant --participant-label ${subj} --fs-license-file /home/k21188249/license/fmriprep/license.txt --output-spaces T1w MNIPediatricAsym:res-1:cohort-3 MyCustom --skip_bids_validation

echo ${subj} finished running fmriprep full

echo deleting ${subj} workflow files

# rm -r /data/project/BIPP/laila/AP_8yo_fmri_task/scripts/FINAL_scripts/logs/work/fmriprep_wf/single_subject_${subj}_wf
### rm -r $PWD/logs/work/fmriprep_wf/single_subject_${subj}_wf ###

echo done deleting ${subj} workflow files

