import warnings, os
import pandas as pd
from os.path import join,exists
import nibabel as nib

import shutil

warnings.filterwarnings("ignore")

import os, sys
from os.path import join, exists, isdir

anapyze_dir = r'/mnt/WORK/repos/anapyze'
anapyze_rsc = join(anapyze_dir,'resources')
sys.path.insert(0,anapyze_dir)

from anapyze_processing.spm import CAT12

csv = '/mnt/nasneuro_share/analysis/jsilva/dataframe_ptau217.csv'
df = pd.read_csv(csv,sep=';')

df_baseline = df.loc[df['mri_bl']==1]


# CONFIG
#---------------------------------------------------------------------------------------------------------------
dir_subjects = r'/mnt/nasneuro_share/data/derivatives/CAT12/PVallecas'
baseline_dir = r'/mnt/nasneuro_share/data/derivatives/CAT12_baseline'

spm_path = '/mnt/WORK/software'
mcr_path = '/mnt/WORK/software'

# Change your templates here if necessary \toolbox\cat12\templates_MNI152NLin2009cAsym

tpm = r'/mnt/WORK/software/spm12/tpm/TPM.nii'
template_volumes = r'/mnt/WORK/software/spm12/toolbox/cat12/templates_MNI152NLin2009cAsym/Template_0_GS.nii'
#---------------------------------------------------------------------------------------------------------------

imgs_to_segment = []

for index_, row_ in df_baseline.iterrows():

    subject = row_['subj']
    visit = row_['vis']

    fsubject = str(subject).zfill(4)
    fvisit = f'V0{visit}'

    nifti_img_orginal = join(dir_subjects,fsubject,fvisit,'T1',f'{fsubject}_T1_{fvisit}.nii')
    format = 'nii'

    if exists(nifti_img_orginal):
        pass
    else:
        nifti_img_orginal = join(dir_subjects,fsubject,fvisit,'T1',f'{fsubject}_T1_{fvisit}.nii.gz')
        format = 'nii.gz'

    if exists(nifti_img_orginal):

        dest_dir = join(baseline_dir,fsubject,fvisit,'T1')
        
        final_gm = join(dest_dir,'mri','mwp1%s_T1_V0%s.nii' % (fsubject,visit))
        print(final_gm)

        if not exists(final_gm):
        
            if not exists(dest_dir):
                os.makedirs(dest_dir)

            nifti_dest = join(dest_dir,f'{fsubject}_T1_{fvisit}.nii')

            if format == 'nii.gz':
                img = nib.load(nifti_img_orginal)
                nib.save(img, nifti_dest)
            else:
                shutil.copy(nifti_img_orginal,nifti_dest)

            print(nifti_dest)

            if exists(nifti_dest):
                print('OK!')

            imgs_to_segment.append(nifti_dest)
    
    else:
        print(f"Final GM already exists: {fsubject}")

spm_proc = CAT12(spm_path=spm_path,mcr_path=mcr_path)

cat_12_cross = spm_proc.cat12seg_imgs(imgs_to_segment, tpm, template_volumes, number_of_cores = 12, surface_processing=0, run=False)






