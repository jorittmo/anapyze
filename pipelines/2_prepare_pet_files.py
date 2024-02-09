import os
from os.path import join

dir_patients = '/mnt/c/Users/jesus/Downloads/Invivo_2024/OUT'

list_dirs = os.listdir(dir_patients)

for i in list_dirs:

    dir_subj = join(dir_patients, i)

    pet_image = False

    files_ = os.listdir(dir_subj)

    for file_ in files_:

        # if file_[0:2] == 'IN':
        if file_[-3:] == 'nii':
            # if 'FDG' in file_:
            pet_image = join(dir_subj, file_)

    if pet_image:
        # pet_dir = join(dir_subj, 'fdg')

        # if not exists(pet_dir):
        #    os.makedirs(pet_dir)

        # if exists(pet_dir):
        #    shutil.rmtree(pet_dir)
        # os.makedirs(pet_dir)

        pet_in = join(dir_subj, 'pet.img')
        os.system('niitoanalyze %s %s' % (pet_image, pet_in))
