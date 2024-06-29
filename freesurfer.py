import os
from concurrent.futures import ThreadPoolExecutor
from os.path import exists, join


class FreeSurfer:
    """
    This class is a wrapper for some FreeSurfer functionality.
    """

    @staticmethod
    def recon_all_cohort_fs(cohort_dir, pats: dict, n_parallel: int = 2) -> None:
        """
        Runs recon-all in parallel for a cohort of subjects.
        :param cohort_dir: Directory containing all subjects' subdirs.
        :param pats: Dictionary of patients subdir names, T1 and T2 file names.
        :param n_parallel: Number of parallel processes to run.
        """

        def process_patient(item: tuple) -> None:

            pat, t1_name, t2_name = item
            pat_dir = join(cohort_dir, pat)
            t1_nii = join(cohort_dir, pat, t1_name)
            t2_nii = join(cohort_dir, pat, t2_name)

            if exists(t1_nii):
                os.system(f"recon-all -sd {pat_dir} -i {t1_nii} -s FS_out -all\n")
            if exists(t1_nii) and exists(t2_nii):
                # TODO : hippocampal subfields
                pass

        with ThreadPoolExecutor(max_workers = n_parallel) as executor:
            executor.map(process_patient, pats.items())

    # TODO : Similar function for SAMSEG

    @staticmethod
    def synthstrip_skull_striping(img_to_strip: str, out_: str = False, includes_csf = True):

        if not out_:
            out_path, out_name = os.path.split(img_to_strip)
            out_: str = join(out_path, 'skull_strip_' + out_name)

        csf_flag = ''
        if not includes_csf:
            csf_flag = '--no-csf'

        if exists(img_to_strip):
            os.system(f"mri_synthstrip -i {img_to_strip} -o {out_} {csf_flag}")

        else:
            raise FileExistsError("Input image does not exist")

    @staticmethod
    def check_freesurfer_env():

        if "FREESURFER_HOME" in os.environ:
            return os.environ["FREESURFER_HOME"]
        else:
            return False