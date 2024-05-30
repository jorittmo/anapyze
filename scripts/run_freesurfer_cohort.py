# /usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Runs Freesurfer basic pipeline (recon_all) iteratively on a cohort directory
"""
import os
import shutil
from os.path import join

# Dir where my patients are
DIR_PATIENTS = "/mnt/lustre/scratch/home/usc/mp/paf/Tau_Project/nifti_files_pga1"
# Basename of T1 images
T1_NAME = "t1.nii"
# TODO:Implement the segmentation with T1+T2

dir_list = sorted(os.listdir(DIR_PATIENTS))

for i in dir_list:
    print(i)
    current_index = dir_list.index(i)

    print("Processing {current_index} patient of {len(dir_list)}: {i}")

    patient_dir = join(DIR_PATIENTS, i)
    output_dir = join(patient_dir, "FS_out")
    my_results_file = join(output_dir, "mri", "aparc+aseg.mgz")

    if os.path.exists(my_results_file):
        print("Patient {i} already has aparc+aseg.mgz")

    else:
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)

        # T1_hdr = join(patient_dir,"MRI.img")
        T1_nii = join(patient_dir, T1_NAME)

        my_script_name = os.path.join(patient_dir, i + "_script.sh")
        with open(my_script_name, "w") as my_script:
            my_script.write("#!/usr/bin/env bash" + "\n")
            my_script.write(f"export SUBJECTS_DIR={patient_dir}\n")
            my_script.write(f"recon-all -i {T1_nii} -s FS_out -all\n")
        my_script.close()

        os.system(f"sh {my_script_name}")
