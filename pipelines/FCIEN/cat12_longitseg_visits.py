import warnings
warnings.filterwarnings("ignore")

import os,sys
import shutil
import gzip
import nibabel as nib
from os.path import join, exists, isdir, basename
import numpy as np
import time

anapyze_dir = r'/Users/jsilva/repositories/anapyze'
anapyze_rsc = join(anapyze_dir,'resources')
sys.path.insert(0,anapyze_dir)
from spm import CAT12

# CONFIG
#---------------------------------------------------------------------------------------------------------------
dir_subjects = r'/Volumes/txusser_data/PVallecas/raws'

# Standalone SPM Paths
spm_path = r'/Users/jsilva/software/CAT12.9/'
mcr_path = r'/Users/jsilva/software/CAT12.9/R2023b'

# Change your templates here if necessary \toolbox\cat12\templates_MNI152NLin2009cAsym
tpm = r'/Users/jsilva/software/spm/tpm/TPM.nii'
template_volumes = r'/Users/jsilva/software/spm12/toolbox/cat12/templates_MNI152NLin2009cAsym/Template_0_GS.nii'
#---------------------------------------------------------------------------------------------------------------

list_subjects = os.listdir(dir_subjects)
list_subjects.sort()

#This will create a cat_12seg_longit.m
images = []

for subj in list_subjects:

    subject = []
    check_cat_processing = []

    dir_subj = join(dir_subjects,subj)

    if isdir(dir_subj):

        list_visits = os.listdir(dir_subj)
        list_visits.sort()

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

spm_proc = CAT12(spm_path,mcr_path)
cat_12_proc = spm_proc.cat12seg_longit(images, tpm, template_volumes, number_of_cores = 12, surface_processing=0, run=False)