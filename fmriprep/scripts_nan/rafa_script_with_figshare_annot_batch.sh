#!/bin/bash
#####
# Transform the parcellation from standard space to native space
# Rafael Romero Garcia
# rr480@cam.ac.uk
# University of Cambridge 2017
#####

#First paratemer $sub is the name of the folder of the subject that want to be processed
#Second parameter is the folder where the freesurfer folders are located

sub=$1 
SUBJECTS_DIR=$2

fsaverage_subid=fsaverage
cd $SUBJECTS_DIR
fsaverage_path=/data/project/BIPP/laila/AP_8yo_fmri_task/derivatives/freesurfer/$fsaverage_subid/
#rm ${SUBJECTS_DIR}/$fsaverage_subid
#ln -s $fsaverage_path ${SUBJECTS_DIR}/$fsaverage_subid

if [ -f ${SUBJECTS_DIR}/${sub}/label/rh.aparc.annot ] && [ ! "${sub}" = 'fsaverage' ]; then
	echo ${SUBJECTS_DIR}/${sub}
	#for parcellation in PALS_B12_Lobes ; do 
	for parcellation in HCPMMP1 ; do 
		for hemi in lh rh ; do
			if [ ! -f ${SUBJECTS_DIR}/${sub}/label/${hemi}.${parcellation}.annot ]; then
			mri_surf2surf --srcsubject ${fsaverage_subid} \
				            --sval-annot ${SUBJECTS_DIR}/${fsaverage_subid}/label/hemi/${hemi}.${parcellation} \
				            --trgsubject ${sub} \
				            --trgsurfval ${SUBJECTS_DIR}/${sub}/label/${hemi}.${parcellation} \
				            --hemi ${hemi}

			fi
		done

		if [ ! -f ${SUBJECTS_DIR}/${sub}/parcellation/${parcellation}.nii.gz ]; then
			mkdir ${SUBJECTS_DIR}/${sub}/parcellation/
			mri_aparc2aseg --s ${sub} \
                        	--o ${SUBJECTS_DIR}/${sub}/parcellation/${parcellation}.nii.gz \
                        	--annot ${parcellation} \
                        	--rip-unknown \
                        	--hypo-as-wm
		fi

		if [ ! -f ${SUBJECTS_DIR}/${sub}/parcellation/${parcellation}_seq.nii.gz ]; then
			echo "addpath('/data/project/neurodev/laila/UCCHILD'); renumDesikan_sub('${SUBJECTS_DIR}/${sub}/parcellation/${parcellation}.nii.gz',1);exit" > temp.m
			matlab -nodisplay -r "temp"
		fi

		for hemi in lh rh ; do
			if [ ! -s ${SUBJECTS_DIR}/${sub}/stats/${hemi}.${parcellation}.log ]; then
		        mris_anatomical_stats -a ${SUBJECTS_DIR}/${sub}/label/${hemi}.${parcellation}.annot -b ${sub} ${hemi} > ${SUBJECTS_DIR}/${sub}/stats/${hemi}.${parcellation}.log
			fi
		done


	done
fi
#rm ${SUBJECTS_DIR}/$fsaverage_subid
