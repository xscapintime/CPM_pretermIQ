# rsync -arv 'k21188249@login1.nan.kcl.ac.uk:/data/project/BIPP/laila/AP_8yo_fmri_task/liyang_files_for_FC/*_framwise_displacement.tsv' .

## 
rsync -arv --exclude 'sub-EPA*/*' 'k21188249@login1.nan.kcl.ac.uk:/data/project/BIPP/laila/AP_8yo_fmri_task/derivatives/fmriprep/sub-*/func/*_framwise_displacement.tsv' .
