#/usr/bin/env python2.7
# -*- coding: utf-8 -*-

import os, sys, shutil, datetime, time
from os.path import join

# Dir where my patients are
dir_patients = "/mnt/lustre/scratch/home/usc/mp/paf/Tau_Project/nifti_files_pga1"

# Cesga Params
cesga_max_time = "30:00:00"
cesga_cores = "1"
cesga_mem = "4096"
max_jobs = 60

# Here is my script

dir_list = sorted(os.listdir(dir_patients))

for i in dir_list:

	print(i)


	os.system("squeue -h -u uscmppaf -o %18i > mycue.txt")
	f = open("mycue.txt", "r")
	lines = f.readlines()
	njobs = len(lines)

	while njobs >= max_jobs:
		os.system("squeue -h -u uscmppaf -o %18i > mycue.txt")
		f = open("mycue.txt", "r")
		lines = f.readlines()
		njobs = len(lines)
		print(njobs)
		print("The cue is full, I am waiting")
		time.sleep(300)


	current_index = dir_list.index(i)

	print("Processing %s patient of %s: %s" % (current_index, len(dir_list), i)) 

	patient_dir = join(dir_patients,i)
	output_dir = join(patient_dir,"FS_out")
	my_results_file = join(output_dir,"mri","aparc+aseg.mgz")

	if os.path.exists(my_results_file):
		print("Patient %s already has aparc+aseg.mgz" % i)

	else:
		if os.path.exists(output_dir):
			shutil.rmtree(output_dir)
	
		#T1_hdr = join(patient_dir,"MRI.img")
		T1_nii = join(patient_dir,"MPRAGE_%s.nii" % i)
		
		my_script_name = os.path.join(patient_dir, i + "_script.sh")
		my_script = open(my_script_name,"w")
		my_script.write("#!/usr/bin/env bash" + "\n")
		my_script.write("export SUBJECTS_DIR=%s\n" % patient_dir)
		#my_script.write("mri_convert %s %s -it analyze -ot nii\n" % (T1_hdr,T1_nii))
		my_script.write("recon-all -i %s -s FS_out -all\n" % (T1_nii))
		#my_script.write("gtmseg --s Freesurfer_out --keep-cc --subsegwm \n")
		#my_script.write("mv %s %s\n" % (my_results_file,patient_dir))
		#my_script.write("tar -czvf %s/FS_out.tgz -C %s .\n" % (patient_dir,output_dir))
		#my_script.write("rm -r %s\n" % (output_dir))
		my_script.close()
    
		print("sbatch -t %s -c %s  --mem=%s --get-user-env --output=%s.out %s" % (
        		cesga_max_time, cesga_cores, cesga_mem, i, my_script_name))
		os.system("sbatch -t %s -c %s  --mem=%s --get-user-env --output=%s.out %s" % (
	       		cesga_max_time, cesga_cores, cesga_mem, i, my_script_name))
