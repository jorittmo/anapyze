"""
I/O helpers for:
  - CAT12 (segmentation scripts)
  - ADNI (data reordering, biomarker extraction, etc.)
  - Generic dcm→nii conversion (via dcm2niix)
  - SPM (coregistration, normalization, smoothing, design matrix, etc.)
"""

from .cat12 import (
    generate_mfile_cat12_segmentation_crossec,
    generate_mfile_cat12_segmentation_longit,
    generate_mfile_cat12_new_tiv_model,
)
from .adni import (
    reorder_ADNI_data,
    filter_ADNI_mri_csv,
    is_ADNI_subject_amyloid_PET_positive,
    get_csf_biomarkers_ADNI,
    get_genetics_data_ADNI,
    get_cognition_data_ADNI,
    get_neuropsychological_battery_ADNI,
    get_wmh_ADNI,
)
from .io import dcm_nii_dcm2niix
from .spm import (
    generate_mfile_coregister,
    generate_mfile_old_normalize,
    generate_mfile_old_deformations,
    generate_mfile_new_normalize,
    generate_mfile_new_deformations,
    generate_mfile_smooth_imgs,
    generate_mfile_model,
    generate_mfile_estimate_model,
    generate_mfile_contrast,
)

__all__ = [
    # from cat12.py
    'generate_mfile_cat12_segmentation_crossec',
    'generate_mfile_cat12_segmentation_longit',
    'generate_mfile_cat12_new_tiv_model',
    # from adni.py
    'reorder_ADNI_data',
    'filter_ADNI_mri_csv',
    'is_ADNI_subject_amyloid_PET_positive',
    'get_csf_biomarkers_ADNI',
    'get_genetics_data_ADNI',
    'get_cognition_data_ADNI',
    'get_neuropsychological_battery_ADNI',
    'get_wmh_ADNI',
    # from io.py
    'dcm_nii_dcm2niix',
    # from spm.py
    'generate_mfile_coregister',
    'generate_mfile_old_normalize',
    'generate_mfile_old_deformations',
    'generate_mfile_new_normalize',
    'generate_mfile_new_deformations',
    'generate_mfile_smooth_imgs',
    'generate_mfile_model',
    'generate_mfile_estimate_model',
    'generate_mfile_contrast',
]