from anapyze.io import adni

# Directory containing the ADNI data just downloaded

dir_data_ = '/home/jsilva/data/ADNI_download'

output_dir_ = '/home/jsilva/data/ADNI_reorder'

# This command will convert the data from DICOM to NIfTI and create a pseudo-BIDS file structure

adni.reorder_ADNI_data(dir_data_, output_dir_, dcm2niix = 'dcm2niix')