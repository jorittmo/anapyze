# Anapyze


[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![GitHub issues](https://img.shields.io/github/issues/Fundacion-CIEN/anapyze)](https://github.com/Fundacion-CIEN/anapyze/issues)
[![Last Commit](https://img.shields.io/github/last-commit/Fundacion-CIEN/anapyze)](https://github.com/Fundacion-CIEN/anapyze/commits/main)

**Anapyze** is a modular Python package providing end-to-end neuroimaging utilities—ranging from core image-processing routines to statistical analysis workflows and I/O wrappers for popular software (SPM, CAT12, FreeSurfer, etc.). 


## Table of Contents

1. [Highlights](#highlights)  
2. [Prerequisites](#prerequisites)  
3. [Installation](#installation)  
4. [Quick Start](#quick-start)  
5. [Core Features](#core-features)  
   - [Core Processing Utilities](#core-processing-utilities)  
   - [Statistical Analysis Routines](#statistical-analysis-routines)  
   - [I/O Helpers](#io-helpers)  
   - [Pipelines (FCIEN & IBIS)](#pipelines-fcien--ibis)  
6. [Usage Examples](#usage-examples)  
7. [Project Structure](#project-structure)  
8. [Development & Contribution](#development--contribution)  
9. [License](#license)  
10. [Contact](#contact)

---

## Highlights

- **“src” Layout**: All code lives under `src/anapyze`.
- **Platforms**: Tested on Linux and macOS; Windows support via WSL or native setup.  

---

## Prerequisites

- **Python**: ≥ 3.6  
- **MATLAB** (R2020a or higher recommended) if running SPM/CAT12 scripts  
- **External Software** (optional, depending on workflow):  
  - [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)   
  - [CAT12](http://www.neuro.uni-jena.de/cat/)   
  - [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/)   
  - [dcm2niix](https://github.com/rordenlab/dcm2niix) 

> **Note**: If you only use NIfTI-level functions (resampling, z-scoring, etc.), MATLAB/SPM is _not_ required.

---

## Installation

1. **Clone the Repository**  
   ```bash
   git clone https://github.com/Fundacion-CIEN/anapyze.git
   cd anapyze
   ```

2. **Install Python Dependencies**  
   ```bash
   pip install --upgrade pip
   pip install -r requirements.txt
   ```

3. **Install Anapyze in Editable Mode**  
   ```bash
   pip install -e .
   ```
   Now you try can `import anapyze` from Python.


## Quick Start

```bash
# Launch Python REPL after editable install
python
```
```python
>>> import anapyze
>>> from anapyze.core import run_matlab_command, coregister_spm
>>> from anapyze.analysis import run_2sample_ttest_spm
>>> from anapyze.io import generate_mfile_coregister

# Example: Run a MATLAB .m script that you already have.
>>> out = run_matlab_command("/path/to/your_spm_batch.m")
# Example: Create a coregistration .m and run it.
>>> coregister_spm(
...     "/data/mean_T1.nii",
...     "/data/func_subject01.nii",
...     "/data/coregister.m",
...     spm_path="/Users/jsilva/software/spm12"
... )
# This will automatically create /data/coregister.m and run it
# A co-registered rfunc_subject01.nii will appear in /data
```

---

## Core Features

### Core Processing Utilities

- **MATLAB Orchestration**:  
  `run_matlab_command(mfile_path, matlab_cmd=...)` – execute `.m` scripts non-interactively.  
- **PET Intensity Normalization**:  
  - `intensity_normalize_pet_histogram(img_path, ref_histogram, ...)`  
  - `intensity_normalize_pet_ref_region(img_path, ref_mask, ...)`  
- **Histogram Matching**:  
  - `histogram_matching(source_img, target_img, ...)`  
  - `logpow_histogram_matching(source_img, target_img, ...)`  
- **SPM-Based Wrappers** (via MATLAB):  
  - `coregister_spm(...)`  
  - `old_normalize_spm(...)` / `new_normalize_spm(...)`  
  - `old_deformations(...)` / `new_deformations(...)`  
  - `smooth_images_spm(img_list, fwhm, output_dir)`  
- **CAT12 & FreeSurfer Helpers**:  
  - `cat12_segmentation_crossec(...)` / `cat12_segmentation_longit(...)`  
  - `recon_all_freesurfer(...)` / `recon_all_freesurfer_whole_cohort(...)`  
  - `synthstrip_skull_striping_freesurfer(...)`  
- **Utility Functions** (in `utils.py`):  
  - `check_input_image_shape(img_path, expected_dims)`  
  - `change_image_dtype(img, new_dtype)`  
  - `resample_image_by_matrix_size(img, new_size)`  
  - `resample_image_by_voxel_sizes(img, new_voxelsize)`  
  - `remove_nan_negs(img_data)`  
  - `add_poisson_noise(img_data, lam)`  
  - `create_mean_std_imgs(img_list, output_dir)`  
  - `create_atlas_csv_from_normals_imgs(img_list, atlas_labels)`  
  - `transform_img_to_atlas_zscores(img, atlas_mask)`  
  - `estimate_fwhm_mizutani(spm_res, dim)`  
  - `spm_map_2_cohens_d(spm_t_map, group_sizes)`  
  - `get_fdr_thresholds_from_spmt(spm_t_map, alpha=0.05)`  
  - `get_tiv_from_cat12_xml_report(xml_report_path)`  
  - `get_weighted_average_iqrs_from_cat12_xml_report(xml_report_path)`

### Statistical Analysis Routines

- **Two-Sample Voxel-Wise t-Test (SPM)**  
  `run_2sample_ttest_spm(spm_path, save_dir, group1, group2, group1_ages, group2_ages, covar2_name=False, group1_covar2=False, group2_covar2=False, mask=None, covar1_name="age")`  
  - Automatically generates an SPM batch (`.m`), calls MATLAB, and saves statistical maps.  
- **Voxel-Wise Correlations**  
  `voxel_wise_corr_images_vs_scale(img_list, scale_scores, brain_mask, output_dir, fdr_alpha=0.05)`  
  - Calculates Pearson’s _r_ across subjects at every voxel versus a continuous measure (e.g., cognitive score).  
  - Returns _r_-map, _p_-map, and FDR-thresholded mask.

### I/O Helpers

- **CAT12 Segmentation Script Generators**  
  - `generate_mfile_cat12_segmentation_crossec(subject_list, img_paths, output_dir)`  
  - `generate_mfile_cat12_segmentation_longit(subject_list, img_paths, timepoints, output_dir)`  
  - `generate_mfile_cat12_new_tiv_model(subject_list, img_paths, output_dir)`  
- **ADNI Utilities**  
  - `reorder_ADNI_data(source_dir, dest_dir)`  
  - `filter_ADNI_mri_csv(csv_path, output_csv)`  
  - `is_ADNI_subject_amyloid_PET_positive(subject_id, pet_csv)`  
  - `get_csf_biomarkers_ADNI(subject_id, csf_csv)`  
  - `get_genetics_data_ADNI(subject_id, genetics_csv)`  
  - `get_cognition_data_ADNI(subject_id, cognition_csv)`  
  - `get_neuropsychological_battery_ADNI(subject_id, neuropsych_csv)`  
  - `get_wmh_ADNI(subject_id, wmh_csv)`  
- **DICOM→NIfTI Conversion**  
  `dcm_nii_dcm2niix(input_dir, output_dir, options=None)`  
- **SPM Batch Generators**  
  - `generate_mfile_coregister(reference, source, output_dir)`  
  - `generate_mfile_old_normalize(img_to_norm, template, output_dir)`  
  - `generate_mfile_old_deformations(deformation_field, output_dir)`  
  - `generate_mfile_new_normalize(img_to_norm, template, output_dir, opts)`  
  - `generate_mfile_new_deformations(deformation_field, output_dir)`  
  - `generate_mfile_smooth_imgs(img_list, fwhm, output_dir)`  
  - `generate_mfile_model(design_mat, contrasts, output_dir)`  
  - `generate_mfile_estimate_model(spm_mat, output_dir)`  
  - `generate_mfile_contrast(spm_mat, contrast_definitions, output_dir)`

### Pipelines (FCIEN)

- **`pipelines/FCIEN`**: Scripts to preprocess DTI and PET for the FCIEN Vallecas cohort (Work in progress)   

Both folders are registered as namespace packages and can be imported (e.g., `import pipelines.FCIEN.run_preprocess_dti`).

---

## Usage Examples

> The following snippets illustrate common workflows. Adapt file paths and parameters to your data.

### 1. Run a Two-Sample t-Test in SPM

```python
from anapyze.analysis import run_2sample_ttest_spm

# Subject lists and covariates
group1_imgs   = ["/data/subj1_IAV.nii", "/data/subj2_IAV.nii"]
group2_imgs   = ["/data/subjA_IAV.nii", "/data/subjB_IAV.nii"]
group1_ages   = [72.3, 68.9]
group2_ages   = [75.1, 70.4]

run_2sample_ttest_spm(
    spm_path="/Applications/MATLAB_R2023b.app/bin/spm",
    save_dir="/results/t_test",
    group1=group1_imgs,
    group2=group2_imgs,
    group1_ages=group1_ages,
    group2_ages=group2_ages,
    covar2_name=False,
    group1_covar2=False,
    group2_covar2=False,
    mask=None,               # e.g. "/templates/GM_mask.nii"
    covar1_name="age"
)
```

### 2. Resample an Image by Voxel Size

```python
from anapyze.core.utils import resample_image_by_voxel_sizes
import nibabel as nib

img          = nib.load("/data/subj01_func.nii")
new_vox_size = (2.0, 2.0, 2.0)  # mm
output_img   = resample_image_by_voxel_sizes(img, new_vox_size)
nib.save(output_img, "/data/subj01_func_resampled.nii")
```

---

## Project Structure

```
anapyze/
├── LICENSE
├── README.md
├── requirements.txt
├── setup.py
├── pipelines/
│   ├── FCIEN/
│   │   ├── __init__.py
│   │   ├── 1_Preprocess_DTI.py
│   │   └── 2_Preprocess_PET.py
│   └── IBIS/
│       ├── __init__.py
│       ├── preprocess_DTI.py
│       ├── preprocess_T1.py
│       └── preprocess_PET.py
└── src/
    └── anapyze/
        ├── __init__.py
        ├── core/
        │   ├── __init__.py
        │   ├── processor.py
        │   └── utils.py
        ├── analysis/
        │   ├── __init__.py
        │   ├── two_samples.py
        │   └── correlations.py
        └── io/
            ├── __init__.py
            ├── adni.py
            ├── cat12.py
            ├── io.py
            └── spm.py
```



## License

This project is distributed under the **MIT License**. See [LICENSE](LICENSE) for full terms.

---

## Contact

- **Maintainer**: Jesús Silva (jesus.bubuchis@gmail.com)  
- **GitHub**: [Fundacion-CIEN/anapyze](https://github.com/Fundacion-CIEN/anapyze)  
- **Issues & Feature Requests**: Use GitHub Issues to report bugs or request enhancements.
