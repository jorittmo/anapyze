import os
from os.path import exists, join
import shutil
import subprocess
from io import spm
from io import cat12
from core import processor
from core import utils
def run_2sample_ttest_spm(spm_path,save_dir,
        group1: list[str],
        group2: list[str],
        group1_ages: list[float],
        group2_ages: list[float],
        covar2_name = False,
        group1_covar2 = False,
        group2_covar2 = False,
        mask: str = False,
        contrast_name: str = "contrast",
        contrast: str = "[1 -1 0 0]",
        ):

    if exists(save_dir):
        shutil.rmtree(save_dir)

    os.makedirs(save_dir)

    print("Creating SPM model....")

    mfile_model = join(save_dir, "model.m")

    spm.generate_mfile_model(spm_path, mfile_model, save_dir, group1, group2,
                         covar1_name='Age', group1_covar1=group1_ages, group2_covar1=group2_ages,
                         covar2_name=covar2_name, group1_covar2=group1_covar2, group2_covar2=group2_covar2,
                         mask=mask)

    processor.run_matlab_command(mfile_model)

    print("Estimating model....")

    mfile_estimate = join(save_dir, "estimate.m")
    spm_mat = join(save_dir, "SPM.mat")

    spm.generate_mfile_estimate_model(mfile_estimate, spm_mat)

    processor.run_matlab_command(mfile_model)

    print("Calculating results....")

    mfile_results = join(save_dir, "results.m")
    spm_mat = join(save_dir, "SPM.mat")

    spm.generate_mfile_contrast(spm_path, mfile_results, spm_mat, contrast_name = contrast_name, contrast = contrast)
    processor.run_matlab_command(mfile_model)

    print("Converting results to Cohens d....")

    out_t_values = join(save_dir, "spmT_0001.nii")
    out_cohens = join(save_dir, "cohens_d.nii")

    utils.spm_map_2_cohens_d(out_t_values, out_cohens, len_1 = len(group1), len_2 = len(group2))

    print("Calculating thresholds for Cohens d (FDR corrected ....")
    t_thres, cohens_d_thres = utils.get_fdr_thresholds_from_spmt(out_t_values, n1 = len(group1), n2 = len(group2))
    print(t_thres,cohens_d_thres)

def run_2sample_ttest_cat12_new_tiv_model(spm_path, save_dir,
        group1: list[str],
        group2: list[str],
        group1_ages: list[float],
        group2_ages: list[float],
        group1_tivs: list[float],
        group2_tivs: list[float],
        mask: str = False,
        contrast_name: str = "contrast",
        contrast: str = "[1 -1 0 0]",
        ):

    if exists(save_dir):
        shutil.rmtree(save_dir)

    os.makedirs(save_dir)

    print("Creating SPM model....")

    mfile_model = join(save_dir, "model_cat12.m")

    cat12.generate_mfile_model(spm_path, save_dir, group1, group1_ages, group1_tivs,
                               group2,group2_ages, group2_tivs, mask)

    processor.run_matlab_command(mfile_model)

    print("Estimating model....")

    mfile_estimate = join(save_dir, "estimate.m")
    spm_mat = join(save_dir, "SPM.mat")

    spm.generate_mfile_estimate_model(mfile_estimate, spm_mat)

    processor.run_matlab_command(mfile_model)

    print("Calculating results....")

    mfile_results = join(save_dir, "results.m")
    spm_mat = join(save_dir, "SPM.mat")

    spm.generate_mfile_contrast(spm_path, mfile_results, spm_mat, contrast_name = contrast_name, contrast = contrast)
    processor.run_matlab_command(mfile_model)

    print("Converting results to Cohens d....")

    out_t_values = join(save_dir, "spmT_0001.nii")
    out_cohens = join(save_dir, "cohens_d.nii")

    utils.spm_map_2_cohens_d(out_t_values, out_cohens, len_1 = len(group1), len_2 = len(group2))

    print("Calculating thresholds for Cohens d (FDR corrected ....")
    t_thres, cohens_d_thres = utils.get_fdr_thresholds_from_spmt(out_t_values, n1 = len(group1), n2 = len(group2))
    print(t_thres,cohens_d_thres)