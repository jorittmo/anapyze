import warnings

warnings.filterwarnings("ignore")

import os, sys
from os.path import join, exists, isdir

anapyze_dir = r'/home/procesadoneuroimagen/WORK/repos/anapyze'
anapyze_rsc = join(anapyze_dir,'resources')
sys.path.insert(0,anapyze_dir)

from anapyze_procesing.spm import CAT12

# CONFIG
#---------------------------------------------------------------------------------------------------------------
dir_subjects = r'/media/procesadoneuroimagen/862A42592A424681'


# CONFIG
# ---------------------------------------------------------------------------------------------------------------
dir_subjects = r'/Volumes/txusser_data/PVallecas/raws'


# Standalone SPM Paths
spm_path = r'/home/procesadoneuroimagen/WORK/software/CAT12.9_R2017b_MCR_Linux'
mcr_path = r'/home/procesadoneuroimagen/WORK/software/CAT12.9_R2017b_MCR_Linux/v93'

# Change your templates here if necessary \toolbox\cat12\templates_MNI152NLin2009cAsym

tpm = r'/home/procesadoneuroimagen/WORK/software/spm12/tpm/TPM.nii'
template_volumes = r'/home/procesadoneuroimagen/WORK/software/spm12/toolbox/cat12/templates_MNI152NLin2009cAsym/Template_0_GS.nii'
#---------------------------------------------------------------------------------------------------------------


list_subjects = os.listdir(dir_subjects)
list_subjects.sort()

# This will create a cat_12seg_longit.m
images = []

for subj in list_subjects:
    print('Checking subject %s' % subj)

    subject = []
    check_cat_processing = []

    dir_subj = join(dir_subjects, subj)

    if isdir(dir_subj):

        list_visits = os.listdir(dir_subj)
        list_visits.sort()

        if len(list_visits)>2:

            for visit in list_visits:

                t1_file = join(dir_subj,visit,'T1','%s_T1_%s.nii' % (subj,visit))

                if exists(t1_file):
                    pass
                else:
                    t1_file = join(dir_subj,visit,'T1','%s_T1_%s.nii.gz' % (subj,visit))

                if exists(t1_file):
                    final_gm = join(dir_subj,visit,'T1','mri','mwmwp1r%s_T1_%s.nii' % (subj,visit))
                    check_cat_processing.append(final_gm)
                    subject.append(t1_file)

            # If any new visit not segmented, the patients will be reprocessed.
            # Probably there is a more efficient way to do this, still figuring out

            if all(exists(j) for j in check_cat_processing):
                pass

            else:
                images.append(subject)
                print("New images found for subject!!")

spm_proc = SPM(spm_path=spm_path,mcr_path=mcr_path)
cat_12_proc = spm_proc.cat12seg_longit(images, tpm, template_volumes, number_of_cores = 12, surface_processing=0, run=False)
