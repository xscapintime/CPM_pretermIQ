# rsync -arv 'k21188249@login1.nan.kcl.ac.uk:/data/project/BIPP/laila/AP_8yo_fmri_task/liyang_files_for_FC/*.txt' .

## 
rsync -arv --exclude 'sub-EAP*/*' 'k21188249@login1.nan.kcl.ac.uk:/data/project/BIPP/laila/AP_8yo_fmri_task/derivatives/freesurfer/sub-*/parcellation/*_timeseries_sept29.txt' .
