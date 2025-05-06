import warnings
import os, sys
import nibabel as nib
from nibabel import processing as nibproc
from os.path import join, exists, isdir
import numpy as np

anapyze_dir = "PLACE HERE YOUR REPO PATH"

sys.path.insert(0,anapyze_dir)

from spm import SPM

warnings.filterwarnings("ignore")

anapyze_rsc = join('/Users/jsilva/repositories/anapyze', 'resources')
dir_patients = r'/Volumes/txusser_data/IBIS_DATA/Reorder_New/'

# CONTINUE HERE.....