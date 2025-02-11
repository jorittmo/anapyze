import os
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from itertools import repeat
from os.path import exists, join
from typing import Tuple
import xmltodict

from applePy import Analysis


class SPM:
    """
    This class will create .m matlab_scripts and parse them to standalone SPM.
    This would allow to more easily update to new SPM/CAT12 versions.
    """

    def __init__(
            self,
            matlab_cmd = '/Applications/MATLAB_R2023b.app/bin/matlab',
            spm_path = '/Users/jsilva/software/spm12'
            ):

        anapyze_proc_dir = os.path.dirname(os.path.abspath(__file__))
        self.template_dir = join(anapyze_proc_dir, 'templates')

        self.matlab_cmd = matlab_cmd
        self.spm_path = spm_path
        self.default_bounding_box = [-84, -120, -72, 84, 84, 96]  # CAT12 Default

    def coregister(
            self,
            reference_image: str,
            source_image: str,
            other_images: list[str] = False,
            wrap = '0 0 0',
            prefix = 'r'
            ) -> str:
        """
        Runs a co-register using spm standalone
        :param reference_image: Input reference image (typically MRI)
        :param source_image: Input source image (typically PET)
        :param other_images: List of other images to co-register (optional)
        :param wrap: Wrapping options. '0 1 0' wraps in y, '1 1 1' wraps the three axes
        :param prefix: Prefix for the output image
        :return: Path to the co-registered image
        """
        template = join(self.template_dir, "coregister_template.m")

        source_img_path, source_img_name = os.path.split(source_image)
        mfile_name: str = join(source_img_path, "coregister.m")

        with open(template) as f:
            text = f.read()
        f.close()

        text = text.replace('SPM_PATH', self.spm_path)
        text = text.replace('REFERENCE_IMAGE', reference_image)
        text = text.replace('SOURCE_IMAGE', source_image)

        other_images_block = "''"
        if other_images:
            other_images_block = '\n'
            for i in other_images:
                other_images_block += "'" + i + ",1'\n"

        text = text.replace('OTHER_IMAGES', other_images_block)
        text = text.replace('WRAP_XYZ', wrap)
        text = text.replace('PREFIX', prefix)

        with open(mfile_name, 'w') as f:
            f.write(text)
        f.close()

        command = f"{self.matlab_cmd} -nosplash -sd {source_img_path} -batch coregister"
        result = subprocess.run(command, shell = True, capture_output = True, text = True)

        if len(result.stderr) == 0:
            print(f"{source_image} coregistered!")
        else:
            raise Exception(result.stderr)

        output: str = os.path.join(source_img_path, prefix + source_img_name)
        return output

    def coregister_parallelization(
            self,
            reference_images: list[str],
            source_images: list[str],
            wrap: str = '0 0 0',
            prefix: str = 'r',
            n_parallel: int = 8
            ) -> None:

        with ThreadPoolExecutor(max_workers = n_parallel) as executor:
            executor.map(self.coregister, reference_images, source_images, repeat(False), repeat(wrap), repeat(prefix))

    def old_normalize(
            self,
            image_to_norm: str,
            template_image: str,
            images_to_write: list[str] = False,
            bb: list[int] = False,
            write_vox_size: float = 1.0,
            wrapping: str = '0 0 0',
            interpolation: str = 'splines',
            prefix: str = 'w'
            ) -> list[str, str]:

        template = join(self.template_dir, "old_normalize_template.m")

        if not bb:
            bb = self.default_bounding_box

        if interpolation == 'nearest':
            interp_val = 0
        elif interpolation == 'trilinear':
            interp_val = 1
        else:
            interp_val = 4  # 4TH DEGREE SPLINES

        source_img_path, source_img_name = os.path.split(image_to_norm)
        mfile_name = join(source_img_path, "old_normalize.m")

        with open(template) as f:
            text = f.read()
        f.close()

        text = text.replace('SPM_PATH', self.spm_path)
        text = text.replace('TEMPLATE', template_image)
        text = text.replace('SOURCE_IMAGE', image_to_norm)
        text = text.replace('WRAP_XYZ', wrapping)
        text = text.replace('PREFIX', prefix)
        text = text.replace('INTERPOLATION', str(interp_val))
        text = text.replace('VOXELSIZE', str(write_vox_size))

        resample_images_block = '\n'
        if images_to_write:
            images_to_write.append(image_to_norm)
        else:
            images_to_write = [image_to_norm]

        for i in images_to_write:
            resample_images_block += "'" + i + ",1'\n"

        bb_block = (str(bb[0]) + " " + str(bb[1]) + " " + str(bb[2]) + "\n"
                    + str(bb[3]) + " " + str(bb[4]) + " " + str(bb[5]))

        text = text.replace('RESAMPLE_IMAGES', resample_images_block)
        text = text.replace('BOUNDING_BOX', bb_block)

        with open(mfile_name, 'w') as f:
            f.write(text)
        f.close()

        command = f"{self.matlab_cmd} -nosplash -sd {source_img_path} -batch old_normalize"
        result = subprocess.run(command, shell = True, capture_output = True, text = True)

        if len(result.stderr) == 0:
            print(f"{image_to_norm} normalized!")
        else:
            raise Exception(result.stderr)

        output = os.path.join(source_img_path + prefix + source_img_name)

        transformation_matrix = image_to_norm[0:-4] + "_sn.mat"

        return output, transformation_matrix

    def old_normalize_parallelization(
            self,
            images_to_norm: list[str],
            template: str,
            bb: list[int] = False,
            write_vox_size: float = 1.0,
            wrapping: str = '0 0 0',
            interpolation: str = 'splines',
            prefix: str = 'w',
            n_parallel: int = 8
            ) -> None:

        with ThreadPoolExecutor(max_workers = n_parallel) as executor:
            executor.map(
                    self.old_normalize, images_to_norm, repeat(template), repeat(False), repeat(bb),
                    repeat(write_vox_size), repeat(wrapping), repeat(interpolation), repeat(prefix)
                    )

    def old_deformations(self, def_matrix: str, base_image: str, images_to_deform: list[str], interpolation: int) \
            -> list[str]:

        source_img_path, source_img_name = os.path.split(images_to_deform[0])
        mfile_name = join(source_img_path, "deformations.m")

        design_type_comp = "matlabbatch{1}.spm.util.defs.comp{1}.inv."
        design_type_out = "matlabbatch{1}.spm.util.defs.out{1}."

        new_spm = open(mfile_name, "w")

        new_spm.write(
                design_type_comp + "comp{1}.sn2def.matname = {'" + def_matrix + "'};\n"
                + design_type_comp + "comp{1}.sn2def.vox = [NaN NaN NaN];\n"
                + design_type_comp + "comp{1}.sn2def.bb = [NaN NaN NaN\n" + "NaN NaN NaN];\n"
                + design_type_comp + "space = {'" + base_image + "'};\n"
                + design_type_out + "pull.fnames = {\n"
                )

        for image in images_to_deform:
            new_spm.write("'" + image + "'" + "\n")
        new_spm.write("};" + "\n")

        new_spm.write(
                design_type_out + "pull.savedir.saveusr = {'" + source_img_path + "/'};\n"
                + design_type_out + "pull.interp = " + str(interpolation) + ";\n"
                + design_type_out + "pull.mask = 1;\n"
                + design_type_out + "pull.fwhm = [0 0 0];\n"
                )

        new_spm.close()

        self.run_mfile(mfile_name)

        deformed_images = []

        for j in images_to_deform:
            components = os.path.split(j)
            output = os.path.join(components[0], "w" + components[1])
            deformed_images.append(output)

        return deformed_images

    def new_normalize(
            self,
            image_to_norm: str,
            images_to_write: list[str] = False,
            bb: list[int] = False,
            write_vox_size: str = "[1 1 1]",
            interpolation: str = 'splines',
            prefix: str = 'w'
            ) -> Tuple[str, str]:

        template = join(self.template_dir, "new_normalize_template.m")

        if not bb:
            bb = self.default_bounding_box

        if interpolation == 'nearest':
            interp_val = 0
        elif interpolation == 'trilinear':
            interp_val = 1
        else:
            interp_val = 4  # 4TH DEGREE SPLINES

        source_img_path, source_img_name = os.path.split(image_to_norm)
        mfile_name = join(source_img_path, "new_normalize.m")

        with open(template) as f:
            text = f.read()
        f.close()

        text = text.replace('SPM_PATH', self.spm_path)
        text = text.replace('IMAGE_TO_NORM', image_to_norm)
        text = text.replace('PREFIX', prefix)
        text = text.replace('INTERPOLATION', str(interp_val))
        text = text.replace('VOXELSIZE', str(write_vox_size))

        resample_images_block = '\n'
        if images_to_write:
            images_to_write.append(image_to_norm)
        else:
            images_to_write = [image_to_norm]

        for i in images_to_write:
            resample_images_block += "'" + i + ",1'\n"

        bb_block = (str(bb[0]) + " " + str(bb[1]) + " " + str(bb[2]) + "\n"
                    + str(bb[3]) + " " + str(bb[4]) + " " + str(bb[5]))

        text = text.replace('RESAMPLE_IMAGES', resample_images_block)
        text = text.replace('BOUNDING_BOX', bb_block)

        with open(mfile_name, 'w') as f:
            f.write(text)
        f.close()

        command = f"{self.matlab_cmd} -nosplash -sd {source_img_path} -batch old_normalize"
        result = subprocess.run(command, shell = True, capture_output = True, text = True)

        if len(result.stderr) == 0:
            print(f"{image_to_norm} normalized!")
        else:
            raise Exception(result.stderr)

        output = os.path.join(source_img_path + prefix + source_img_name)

    def new_deformations(
            self,
            def_matrix: str,
            images_to_deform: list[str],
            interpolation: int,
            prefix = "w"
            ) \
            -> list[str]:

        source_img_path, source_img_name = os.path.split(images_to_deform[0])
        # Set the output file name
        mfile_name = join(source_img_path, "deformations.m")

        design_type_comp = "matlabbatch{1}.spm.util.defs.comp{1}."
        design_type_out = "matlabbatch{1}.spm.util.defs.out{1}."

        new_spm = open(mfile_name, "w")

        new_spm.write(
                design_type_comp + "def = {'" + def_matrix + "'};\n"
                + design_type_out + "pull.fnames = {" + "\n"
                )

        for image in images_to_deform:
            new_spm.write("'" + image + "'\n")
        new_spm.write("};\n")

        new_spm.write(
                design_type_out + "pull.savedir.savesrc = 1;\n"
                + design_type_out + "pull.interp =" + str(interpolation) + ";\n"
                + design_type_out + "pull.mask = 0;\n"
                + design_type_out + "pull.fwhm = [0 0 0];\n"
                + design_type_out + "pull.prefix ='" + prefix + "';\n"
                )

        new_spm.close()

        deformed_images = []

        for j in images_to_deform:
            components = os.path.split(j)
            output = os.path.join(components[0], "w" + components[1])
            deformed_images.append(output)

        return deformed_images

    def apply_normalization_to_atlas(self, def_matrix: str, norm_mri: str, fs_atlas: str) -> str:
        source_img_path, source_img_name = os.path.split(fs_atlas)
        # Set the output file name
        mfile_name = join(source_img_path, "deformations.m")

        design_type = "matlabbatch{1}.spm.util.defs."

        new_spm = open(mfile_name, "w")

        new_spm.write(
                design_type + "comp{1}.inv.comp{1}.def = {'" + def_matrix + "'};\n"
                )
        new_spm.write(design_type + "comp{1}.inv.space = {'" + norm_mri + "'};\n")
        new_spm.write(design_type + "out{1}.push.fnames = {'" + fs_atlas + "'};\n")
        new_spm.write(design_type + "out{1}.push.weight = {''};\n")
        new_spm.write(design_type + "out{1}.push.savedir.savesrc = 1;\n")
        new_spm.write(design_type + "out{1}.push.fov.file = {'" + norm_mri + "'};\n")
        new_spm.write(design_type + "out{1}.push.preserve = 2;\n")
        new_spm.write(design_type + "out{1}.push.fwhm = [0 0 0];\n")
        new_spm.write(design_type + "out{1}.push.prefix = 'w';\n")

        new_spm.close()

        self.run_mfile(mfile_name)

        components = os.path.split(fs_atlas)
        output = os.path.join(components[0], "w" + components[1])

        return output

    def smooth_imgs(self, dir_proc: str, images_to_smooth: list[str], smoothing: list[int]) -> None:

        source_img_path, source_img_name = os.path.split(images_to_smooth[0])
        # Set the output file name
        mfile_name = join(dir_proc, "smooth.m")

        design_type = "matlabbatch{1}.spm.spatial.smooth."
        smoothing_array = (
                "[" + str(smoothing[0]) + " " + str(smoothing[1]) + " " + str(smoothing[2]) + "]"
        )

        new_spm = open(mfile_name, "w")

        new_spm.write(design_type + "data = {\n")

        for i in images_to_smooth:
            new_spm.write("'" + i + ",1'\n")

        new_spm.write("};" + "\n")
        new_spm.write(design_type + "fwhm =" + smoothing_array + ";" + "\n")
        new_spm.write(design_type + "dtype = 0;" + "\n")
        new_spm.write(design_type + "im = 0;" + "\n")
        new_spm.write(design_type + "prefix ='" + "s" + "';" + "\n")

        new_spm.close()

        self.run_mfile(mfile_name)

    def run_2sample_ttest(
            self,
            save_dir: str,
            group1: list[str],
            group2: list[str],
            group1_ages: list[float],
            group2_ages: list[float],
            group1_tiv: list[float] = False,
            group2_tiv: list[float] = False,
            mask: str = False,
            contrast_name: str = "contrast",
            contrast: str = "[1 -1 0 0]",
            ):

        if exists(save_dir):
            shutil.rmtree(save_dir)

        os.makedirs(save_dir)

        print("Creating SPM model....")

        mfile_model = join(save_dir, "model.m")

        self.create_mfile_model(
                mfile_model, save_dir, group1, group2, group1_ages, group2_ages, group1_tiv, group2_tiv, mask
                )

        command = f"{self.matlab_cmd} -nosplash -sd {save_dir} -batch model"
        result = subprocess.run(command, shell = True, capture_output = True, text = True)

        print(result.stdout)
        print(result.stderr)

        print("Estimating model....")

        mfile_estimate = join(save_dir, "estimate.m")
        spm_mat = join(save_dir, "SPM.mat")

        self.create_mfile_estimate_model(mfile_estimate, spm_mat)

        command = f"{self.matlab_cmd} -nosplash -sd {save_dir} -batch estimate"
        result = subprocess.run(command, shell = True, capture_output = True, text = True)

        print(result.stdout)
        print(result.stderr)

        print("Calculating results....")

        mfile_results = join(save_dir, "results.m")
        spm_mat = join(save_dir, "SPM.mat")

        self.create_mfile_contrast(mfile_results, spm_mat, contrast_name = contrast_name, contrast = contrast)

        command = f"{self.matlab_cmd} -nosplash -sd {save_dir} -batch results"
        result = subprocess.run(command, shell = True, capture_output = True, text = True)

        print(result.stdout)
        print(result.stderr)

        print("Converting results to Cohens d....")

        out_t_values = join(save_dir, "spmT_0001.nii")
        out_cohens = join(save_dir, "cohens_d.nii")

        Analysis.spm_map_2_cohens_d(out_t_values, out_cohens, len(group1), len(group2))

        print("Calculating thresholds for Cohens d (FDR corrected ....")
        Analysis.get_cohens_d_thresholds_fdr(out_t_values, len(group1), len(group2))

    def create_mfile_model(self,
            mfile_name: str,
            save_dir: str,
            group1: list[str],
            group2: list[str],
            group1_ages: list[float],
            group2_ages: list[float],
            group1_tiv: list[float] = False,
            group2_tiv: list[float] = False,
            mask: str = False,
            ) -> None:


        new_spm = open(mfile_name, "w")

        new_spm.write(f"addpath ('{self.spm_path}');\n\n")

        design_type = "matlabbatch{1}.spm.stats.factorial_design."

        new_spm.write(
                design_type + "dir = {'" + save_dir + "/'};\n"
                + design_type + "des.t2.scans1 = {\n"
                )

        for image in group1:
            new_spm.write("'" + image + ",1'" + "\n")
        new_spm.write("};" + "\n")

        new_spm.write(design_type + "des.t2.scans2 = {" + "\n")

        for image in group2:
            new_spm.write("'" + image + ",1'" + "\n")
        new_spm.write("};" + "\n")

        new_spm.write(
                design_type
                + "des.t2.dept = 0;\n"
                + design_type + "des.t2.variance = 1;\n"
                + design_type + "des.t2.gmsca = 0;\n"
                + design_type + "des.t2.ancova = 0;\n"
                )

        new_spm.write(design_type + "cov(1).c = [")
        for age in group1_ages:
            new_spm.write(str(age) + "\n")
        for age in group2_ages:
            new_spm.write(str(age) + "\n")
        new_spm.write("];\n")

        new_spm.write(
                design_type + "cov(1).cname = 'Age';\n"
                + design_type + "cov(1).iCFI = 1;\n"
                + design_type + "cov(1).iCC = 5;\n"
                )

        if group1_tiv:
            new_spm.write(design_type + "cov(2).c = [")
            for tiv in group1_tiv:
                new_spm.write(str(tiv) + "\n")
            for tiv in group2_tiv:
                new_spm.write(str(tiv) + "\n")
            new_spm.write("];\n")

            new_spm.write(
                    design_type + "cov(2).cname = 'TIV';\n"
                    + design_type + "cov(2).iCFI = 1;\n"
                    + design_type + "cov(2).iCC = 1;\n"
                    )

        new_spm.write(
                design_type + "multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});\n"
                + design_type + "masking.tm.tm_none = 1;\n"
                + design_type + "masking.im = 0;\n"
                + design_type + "masking.em = {'" + mask + ",1'};\n"
                + design_type + "globalc.g_omit = 1;\n"
                + design_type + "globalm.gmsca.gmsca_no = 1;\n"
                + design_type + "globalm.glonorm = 1;\n\n"
                )

        new_spm.write("spm('defaults','fmri');\n")
        new_spm.write("spm_jobman('initcfg');\n")
        new_spm.write("spm_jobman('run',matlabbatch);\n")

        new_spm.close()

    def create_mfile_estimate_model(self, mfile_name: str, spm_mat: str) -> None:

        new_spm = open(mfile_name, "w")

        new_spm.write(f"addpath ('{self.spm_path}');\n\n")

        design_type = "matlabbatch{1}.spm.stats.fmri_est."

        new_spm.write(design_type + "spmmat = {'" + spm_mat + "'};\n")
        new_spm.write(design_type + "write_residuals = 0;\n")
        new_spm.write(design_type + "method.Classical = 1;\n\n")

        new_spm.write("spm('defaults','fmri');\n")
        new_spm.write("spm_jobman('initcfg');\n")
        new_spm.write("spm_jobman('run',matlabbatch);\n")

        new_spm.close()

    def create_mfile_contrast(
            self, mfile_name: str, spm_mat: str, contrast_name: str = "contrast", contrast: str = "[1 -1 0]"
            ) -> None:

        new_spm = open(mfile_name, "w")

        new_spm.write(f"addpath ('{self.spm_path}');\n\n")

        design_type = "matlabbatch{1}.spm.stats.con."

        new_spm.write(design_type + "spmmat = {'" + spm_mat + "'};\n")
        new_spm.write(design_type + "consess{1}.tcon.name = '" + contrast_name + "';\n")
        new_spm.write(design_type + "consess{1}.tcon.weights =" + contrast + ";\n")
        new_spm.write(design_type + "consess{1}.tcon.sessrep = 'none';\n")
        new_spm.write(design_type + "delete = 0;\n\n")

        new_spm.write("spm('defaults','fmri');\n")
        new_spm.write("spm_jobman('initcfg');\n")
        new_spm.write("spm_jobman('run',matlabbatch);\n")

        new_spm.close()


class CAT12:
    """
    This class will create .m matlab_scripts and parse them to matlab to run SPM.
    I have modified this class to make use of matlab instead os SPM standalone.
    This would allow to more easily update to new SPM versions.
    Also should be the class more directly usable in other OS such as MAC_OS or Linux
    However, this requires matlab to be in your $PATH.
    """

    def __init__(
            self,
            spm_path = "/home/jsilva/software/cat12_standalone",
            mcr_path = "/home/jsilva/software/MATLAB_MCR/v93",
            ):
        self.spm_path = spm_path
        self.mcr_path = mcr_path

        if not exists(self.spm_path):
            raise FileNotFoundError(f"{self.spm_path} is not found.")

        if not exists(self.spm_path):
            raise FileNotFoundError(f"{self.mcr_path} is not found.")

        self.spm_run = "%s/run_spm12.sh %s batch" % (self.spm_path, self.mcr_path)

    def run_mfile(self, mfile: str) -> None:
        os.system("%s %s" % (self.spm_run, mfile))

    def cat12seg_imgs(
            self,
            images_to_seg: list[str],
            template_tpm: str,
            template_volumes: str,
            number_of_cores: int = 0,
            output_vox_size: float = 1.5,
            bounding_box: str = "cat12",
            surface_processing: int = 0,
            run = False,
            ) -> str:
        """
        This function creates a mfile to later run with MATLAB.
        You can run it within spm.py stating run=True but multithreading is not available.
        mfile: Destination of the created mfile
        images_to_seg = A list of all the images you want to segment (Nifti (.nii) of Analyze (.img)).
        template_tpm = Template image to normalize the images to (something like ... PATH_TO/spm12/tpm/TPM.nii)
        template_volumes = Templete volumes image from CAT12 (something like ...
        PATH_TO/spm12/toolbox/cat12/template_volumes/Template_0_IXI555_MNI152_GS.nii)
        rest of the parameters = Consult CAT12 segment
        """

        design_type = "matlabbatch{1}.spm.tools.cat.estwrite."

        mfile_name = "cat12seg.m"
        new_spm = open(mfile_name, "w")

        new_spm.write(design_type + "data = {\n")
        for i in range(len(images_to_seg)):
            new_spm.write("'" + images_to_seg[i] + ",1'\n")
        new_spm.write("};" + "\n")

        new_spm.write(design_type + "data_wmh = {''};" + "\n")
        new_spm.write(design_type + "nproc = " + str(number_of_cores) + ";\n")
        new_spm.write(design_type + "useprior = '';" + "\n")
        new_spm.write(design_type + "opts.tpm = {'" + template_tpm + "'};\n")

        design_type = "matlabbatch{1}.spm.tools.cat.estwrite.opts."

        new_spm.write(design_type + "affreg = 'mni';" + "\n")
        new_spm.write(design_type + "biasacc = " + "0.5" + ";\n")
        new_spm.write(design_type + "accstr = " + "0.5" + ";\n")

        design_type = "matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation."

        new_spm.write(design_type + "restypes.optimal = " + "[1 0.3]" + ";\n")
        new_spm.write(design_type + "setCOM = " + "1" + ";\n")
        new_spm.write(design_type + "APP = " + "1070" + ";\n")
        new_spm.write(design_type + "affmod = " + "0" + ";\n")
        new_spm.write(design_type + "NCstr = " + "-Inf" + ";\n")
        new_spm.write(design_type + "spm_kamap = " + "0" + ";\n")
        new_spm.write(design_type + "LASstr = " + "0.5" + ";\n")
        new_spm.write(design_type + "LASmyostr = " + "0" + ";\n")
        new_spm.write(design_type + "gcutstr = " + "2" + ";\n")
        new_spm.write(design_type + "cleanupstr = " + "0.5" + ";\n")
        new_spm.write(design_type + "BVCstr = " + "0.5" + ";\n")
        new_spm.write(design_type + "WMHC = " + "2" + ";\n")
        new_spm.write(design_type + "SLC = " + "0" + ";\n")
        new_spm.write(design_type + "mrf = " + "1" + ";\n")

        design_type = "matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration."

        new_spm.write(
                design_type
                + "regmethod.shooting.shootingtpm = {'"
                + template_volumes
                + "'};\n"
                )
        new_spm.write(design_type + "regmethod.shooting.regstr = " + "0.5" + ";\n")
        new_spm.write(design_type + "vox = " + str(output_vox_size) + ";\n")

        if bounding_box == "cat12":
            new_spm.write(design_type + "bb = " + "12" + ";\n")
        elif bounding_box == "spm":
            new_spm.write(design_type + "bb = " + "16" + ";\n")
        else:
            raise ValueError("Bounding Box must be set to spm or cat12")

        design_type = "matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface."

        new_spm.write(design_type + "pbtres = " + "0.5" + ";\n")
        new_spm.write(design_type + "pbtmethod = " + "'pbt2x'" + ";\n")
        new_spm.write(design_type + "SRP = " + "22" + ";\n")
        new_spm.write(design_type + "reduce_mesh = " + "1" + ";\n")
        new_spm.write(design_type + "vdist = " + "2" + ";\n")
        new_spm.write(design_type + "scale_cortex = " + "0.7" + ";\n")
        new_spm.write(design_type + "add_parahipp = " + "0.1" + ";\n")
        new_spm.write(design_type + "close_parahipp = " + "1" + ";\n")

        design_type = "matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin."

        new_spm.write(design_type + "experimental = " + "0" + ";\n")
        new_spm.write(design_type + "new_release = " + "0" + ";\n")
        new_spm.write(design_type + "lazy = " + "0" + ";\n")
        new_spm.write(design_type + "ignoreErrors = " + "1" + ";\n")
        new_spm.write(design_type + "verb = " + "2" + ";\n")
        new_spm.write(design_type + "print = " + "2" + ";\n")

        design_type = "matlabbatch{1}.spm.tools.cat.estwrite.output."

        new_spm.write(design_type + "BIDS.BIDSno = " + "1" + ";\n")
        new_spm.write(design_type + "surface = " + str(surface_processing) + ";\n")
        new_spm.write(
                design_type + "surf_measures = " + str(surface_processing) + ";\n"
                )

        design_type = "matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases."

        new_spm.write(design_type + "neuromorphometrics = " + "0" + ";\n")
        new_spm.write(design_type + "lpba40 = 0;\n")
        new_spm.write(design_type + "cobra = 0;\n")
        new_spm.write(design_type + "hammers = 0;\n")
        new_spm.write(design_type + "thalamus = 0;\n")
        new_spm.write(design_type + "thalamic_nuclei = 0;\n")
        new_spm.write(design_type + "suit = 0;\n")
        new_spm.write(design_type + "ibsr = 0;\n")
        new_spm.write(design_type + "aal3 = 0;\n")
        new_spm.write(design_type + "mori = 0;\n")
        new_spm.write(design_type + "anatomy3 = 0;\n")
        new_spm.write(design_type + "julichbrain = 0;\n")
        new_spm.write(design_type + "Schaefer2018_100Parcels_17Networks_order = 0;\n")
        new_spm.write(design_type + "Schaefer2018_200Parcels_17Networks_order = 0;\n")
        new_spm.write(design_type + "Schaefer2018_400Parcels_17Networks_order = 0;\n")
        new_spm.write(design_type + "Schaefer2018_600Parcels_17Networks_order = 0;\n")
        new_spm.write(design_type + "ownatlas = " + "{''}" + ";\n")

        design_type = "matlabbatch{1}.spm.tools.cat.estwrite.output."

        new_spm.write(design_type + "GM.native = 0;\n")
        new_spm.write(design_type + "GM.warped = 0;\n")
        new_spm.write(design_type + "GM.mod = 1;\n")
        new_spm.write(design_type + "GM.dartel = 0;\n")

        new_spm.write(design_type + "WM.native = 0;\n")
        new_spm.write(design_type + "WM.warped = 0;\n")
        new_spm.write(design_type + "WM.mod = 1;\n")
        new_spm.write(design_type + "WM.dartel = 0;\n")

        new_spm.write(design_type + "CSF.native = 0;\n")
        new_spm.write(design_type + "CSF.warped = 0;\n")
        new_spm.write(design_type + "CSF.mod = 1;\n")
        new_spm.write(design_type + "CSF.dartel = 0;\n")

        new_spm.write(design_type + "ct.native = 0;\n")
        new_spm.write(design_type + "ct.warped = 0;\n")
        new_spm.write(design_type + "ct.dartel = 0;\n")

        new_spm.write(design_type + "pp.native = 0;\n")
        new_spm.write(design_type + "pp.warped = 0;\n")
        new_spm.write(design_type + "pp.dartel = 0;\n")

        new_spm.write(design_type + "WMH.native = 0;\n")
        new_spm.write(design_type + "WMH.warped = 0;\n")
        new_spm.write(design_type + "WMH.mod = 0;\n")
        new_spm.write(design_type + "WMH.dartel = 0;\n")

        new_spm.write(design_type + "SL.native = 0;\n")
        new_spm.write(design_type + "SL.warped = 0;\n")
        new_spm.write(design_type + "SL.mod = 0;\n")
        new_spm.write(design_type + "SL.dartel = 0;\n")

        new_spm.write(design_type + "TPMC.native = 0;\n")
        new_spm.write(design_type + "TPMC.warped = 0;\n")
        new_spm.write(design_type + "TPMC.mod = 0;\n")
        new_spm.write(design_type + "TPMC.dartel = 0;\n")

        new_spm.write(design_type + "atlas.native = 0;\n")
        new_spm.write(design_type + "label.native = 1;\n")
        new_spm.write(design_type + "label.warped = 0;\n")
        new_spm.write(design_type + "label.dartel = 0;\n")
        new_spm.write(design_type + "labelnative = 1;\n")

        new_spm.write(design_type + "bias.native = 0;\n")
        new_spm.write(design_type + "bias.warped = 1;\n")
        new_spm.write(design_type + "bias.dartel = 0;\n")

        new_spm.write(design_type + "las.native = 0;\n")
        new_spm.write(design_type + "las.warped = 1;\n")
        new_spm.write(design_type + "las.dartel = 0;\n")

        new_spm.write(design_type + "jacobianwarped = 0;\n")
        new_spm.write(design_type + "warps = [1 0];\n")
        new_spm.write(design_type + "rmat = 0;\n")

        new_spm.close()

        if run:
            self.run_mfile(mfile_name)

        return mfile_name

    def cat12seg_longit(
            self,
            images_to_seg: list[str],
            template_tpm: str,
            template_volumes: str,
            longmodel: int = 2,
            number_of_cores: int = 2,
            output_vox_size: float = 1.5,
            bounding_box: str = "cat12",
            surface_processing: int = 0,
            run = False,
            ):
        """
        This function creates a mfile to later run with MATLAB.
        You can run it within spm.py stating run=True but multithreading is not available.
        mfile: Destination of the created mfile
        images_to_seg = A list of all the images of all the subjects you want to segment (Nifti (.nii/.nii.gz) of
        Analyze (.img)).
        images_to_seg should be provided as [[S1_image1,S1_image2,...],[S2_image1,S2_image2,....]
        template_tpm = Template image to normalize the images to (something like ... PATH_TO/spm12/tpm/TPM.nii)
        template_volumes = Templete volumes image from CAT12 (something like ...
        PATH_TO/spm12/toolbox/cat12/template_volumes/Template_0_IXI555_MNI152_GS.nii)
        rest of the parameters = Consult CAT12 segment
        """

        design_type = "matlabbatch{1}.spm.tools.cat.long."

        mfile_name = "cat12seg_longit.m"
        new_spm = open(mfile_name, "w")

        new_spm.write(design_type + "datalong.subjects = {\n")
        for i in range(len(images_to_seg)):
            new_spm.write("{\n")
            for j in range(len(images_to_seg[i])):
                new_spm.write("'" + images_to_seg[i][j] + "'\n")
            new_spm.write("}\n")
        new_spm.write("};\n")

        new_spm.write(design_type + "longmodel = " + str(longmodel) + ";\n")
        new_spm.write(design_type + "enablepriors = 1;\n")
        new_spm.write(design_type + "prepavg = 2;\n")
        new_spm.write(design_type + "bstr = 0;\n")
        new_spm.write(design_type + "avgLASWMHC = 0;\n")
        new_spm.write(design_type + "nproc = " + str(number_of_cores) + ";\n")
        new_spm.write(design_type + "opts.tpm = {'" + template_tpm + "'};\n")
        new_spm.write(design_type + "opts.affreg = 'mni';\n")
        new_spm.write(design_type + "opts.biasacc = 0.5;\n")

        design_type = "matlabbatch{1}.spm.tools.cat.long.extopts."

        new_spm.write(design_type + "restypes.optimal = [1 0.3];\n")
        new_spm.write(design_type + "setCOM = 1;\n")
        new_spm.write(design_type + "APP = 1070;\n")
        new_spm.write(design_type + "affmod = 0;\n")
        new_spm.write(design_type + "LASstr = 0.5;\n")
        new_spm.write(design_type + "LASmyostr = 0;\n")
        new_spm.write(design_type + "gcutstr = 2;\n")
        new_spm.write(design_type + "WMHC = 2;\n")
        new_spm.write(design_type + "registration.shooting.shootingtpm = {'" + template_volumes + "'};\n")
        new_spm.write(design_type + "registration.shooting.regstr = 0.5;\n")
        new_spm.write(design_type + "vox = " + str(output_vox_size) + ";\n")

        if bounding_box == "cat12":
            new_spm.write(design_type + "bb = 12;\n")
        elif bounding_box == "spm":
            new_spm.write(design_type + "bb = 16;\n")
        else:
            raise ValueError("Bounding Box must be set to spm or cat12")

        new_spm.write(design_type + "SRP = 22;\n")
        new_spm.write(design_type + "ignoreErrors = 1;\n")

        design_type = "matlabbatch{1}.spm.tools.cat.long."

        new_spm.write(design_type + "output.BIDS.BIDSno = 1;\n")
        new_spm.write(design_type + "output.surface = " + str(surface_processing) + ";\n")
        new_spm.write(design_type + "ROImenu.noROI = struct([]);\n")
        new_spm.write(design_type + "longTPM = 1;\n")
        new_spm.write(design_type + "modulate = 1;\n")
        new_spm.write(design_type + "dartel = 0;\n")
        new_spm.write(design_type + "printlong = 2;\n")
        new_spm.write(design_type + "delete_temp = 1;\n")

        new_spm.write("\n")
        new_spm.write("spm('defaults','fmri');\n")
        new_spm.write("spm_jobman('initcfg');\n")
        new_spm.write("spm_jobman('run',matlabbatch);\n")

        new_spm.close()

        if run:
            self.run_mfile(mfile_name)

        return mfile_name

    def create_cat12_new_model(
            self,
            save_dir: str,
            group1: list[str],
            group1_ages: list[float],
            group1_tivs: list[float],
            group2: list[str],
            group2_ages: list[float],
            group2_tivs: list[float],
            mask: str,
            ):
        if exists(save_dir):
            shutil.rmtree(save_dir)

        os.makedirs(save_dir)

        mfile_name = join(save_dir, "cat_12_vbm.m")
        design_type = "matlabbatch{1}.spm.tools.cat.factorial_design."

        new_spm = open(mfile_name, "w")

        new_spm.write("addpath('%s');\n" % self.spm_path)

        new_spm.write(
                design_type + "dir = {'" + save_dir + "/'};\n"
                + design_type + "des.t2.scans1 = {\n"
                )

        for image in group1:
            new_spm.write("'" + image + ",1'\n")
        new_spm.write("};\n")

        new_spm.write(design_type + "des.t2.scans2 = {\n")

        for image in group2:
            new_spm.write("'" + image + ",1'\n")
        new_spm.write("};\n")

        new_spm.write(
                design_type + "des.t2.dept = 0;\n"
                + design_type + "des.t2.variance = 1;\n"
                + design_type + "des.t2.gmsca = 0;\n"
                + design_type + "des.t2.ancova = 0;\n"
                )

        new_spm.write(design_type + "cov.c = [")
        for age in group1_ages:
            new_spm.write(str(age) + "\n")
        for age in group2_ages:
            new_spm.write(str(age) + "\n")
        new_spm.write("];\n")

        new_spm.write(
                design_type + "cov.cname = 'Age';\n"
                + design_type + "cov.iCFI = 1;\n"
                + design_type + "cov.iCC = 5;\n"
                )

        new_spm.write(
                design_type + "multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});\n"
                )

        new_spm.write(
                design_type + "masking.tm.tm_none = 1;\n"
                + design_type + "masking.im = 1;\n"
                )

        new_spm.write(design_type + "masking.em = {'" + mask + "'}\n;")

        new_spm.write(design_type + "globals.g_ancova.global_uval = [")
        for tiv in group1_tivs:
            new_spm.write(str(tiv) + "\n")
        for tiv in group2_tivs:
            new_spm.write(str(tiv) + "\n")
        new_spm.write("];\n")

        new_spm.write(
                design_type + "check_SPM.check_SPM_zscore.do_check_zscore.use_unsmoothed_data = 1;\n"
                + design_type + "check_SPM_zscore.do_check_zscore.adjust_data = 1;\n"
                + design_type + "check_SPM.check_SPM_ortho = 1;\n"
                )

        spm_mat = join(save_dir, "SPM.mat")
        design_type = "matlabbatch{2}.spm.stats.fmri_est."

        new_spm.write(design_type + "spmmat = {'" + spm_mat + "'};\n")
        new_spm.write(design_type + "write_residuals = 0;")
        new_spm.write(design_type + "method.Classical = 1;")

        design_type = "matlabbatch{3}.spm.stats.con."
        contrast_name = "Atrophy"
        contrast = "[1 -1 0 0]"

        new_spm.write(design_type + "spmmat = {'" + spm_mat + "'};\n")
        new_spm.write(design_type + "consess{1}.tcon.name = '" + contrast_name + "';\n")
        new_spm.write(design_type + "consess{1}.tcon.weights =" + contrast + ";\n")
        new_spm.write(design_type + "consess{1}.tcon.sessrep = 'none';\n")
        new_spm.write(design_type + "delete = 0;\n")

        new_spm.close()

    @staticmethod
    def get_tiv_from_xml_report(cat_xml_filepath: str) -> float:
        """parse the information from a list-like object of "cat_*.xml" filepaths to
        a list of dictionaries for more easy data handling"""

        with open(cat_xml_filepath) as file:
            cat_xml_dict = xmltodict.parse(file.read())

            try:
                vol_tiv = cat_xml_dict["S"]["subjectmeasures"]["vol_TIV"]
                return vol_tiv

            except KeyError:
                raise KeyError("Could not extract TIV")

    @staticmethod
    def get_weighted_average_iqrs_from_xml_report(cat_xml_filepath: str) -> float:
        """Extract Weighted Average IQRS from dictionaries produced by _parse_xml_files_to_dict"""

        with open(cat_xml_filepath) as file:
            cat_xml_dict = xmltodict.parse(file.read())

        pattern = "(Image Quality Rating \(IQR\):) (\d\d.\d\d)% \(([A-Z](-|\+)?)\)"

        try:
            catlog = cat_xml_dict["S"]["catlog"]["item"]
            for e in catlog:
                if e.startswith("Image Quality Rating (IQR)"):
                    match = re.search(pattern, e)
                    weighted_average_iqr = match.group(2)
                    weighted_average_iqr = float(weighted_average_iqr)
                    return weighted_average_iqr

        except:
            raise KeyError("Could not extract IQRS")