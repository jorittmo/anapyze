import os
from os.path import join, exists

dir_patients = '/mnt/d/IBIS_DATA/Reorder_New'

list_dirs = os.listdir(dir_patients)

for i in list_dirs:

    dir_subj = join(dir_patients, i)

    files_ = os.listdir(dir_subj)

    for file_ in files_:

        if file_[0:2] == 'IN' and file_[-3:] == 'nii' and 'FDG' in file_:

            pet_in = join(dir_subj, file_)
            date = file_[-12:-4]

            pet_dir = join(dir_subj, 'fdg_%s' % date)

            if not exists(pet_dir):
                os.makedirs(pet_dir)

            pet_ = join(pet_dir, 'pet.img')
            os.system('niitoanalyze %s %s' % (pet_in, pet_))
            print("Created %s" % pet_)
