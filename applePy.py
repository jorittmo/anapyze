import os
import subprocess
from os.path import dirname, exists, join
from textwrap import dedent

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd
import pingouin as pg
import pwlf
from scipy.fftpack import fftn, fftshift
from scipy.ndimage import zoom
from scipy.stats import pearsonr, spearmanr, t

import spm


class Format:
    """
    Tools to convert data from and to different formats (Nifti, Analyze, DICOM)
    """

    @staticmethod
    def nii_hdr_convert(img_name, compress_out=True):
        """
        This function converts an image file from Nifti to Analyze vice versa.
        Depending on your input, it will output the other format.
        It also works with .nii.gz

        :param img_name: the name of the nii/nii.gz/hdr/img to be converted
        :type img_name: str
        :param compress_out: Will output nii.gz if the input is Analyze
        :type compress_out: bool
        :return: the name of the result file that was created
        :rtype: str
        """

        # Load the NIFTI image using nibabel
        img, data = Utils.load_nifti(img_name)

        # Determine the output filename based on the input file format
        if img_name[-4:] == ".hdr" or img_name[-4:] == ".img":
            # Convert from Analyze to NIFTI
            if compress_out:
                out_name = img_name[:-4] + ".nii.gz"
            else:
                out_name = img_name[:-4] + ".nii"

            Utils.save_nifti(img, data, out_name)
            print(
                f"File {img_name} was converted to Nifti format and saved as {out_name}"
            )
            return out_name

        elif img_name[-4:] == ".nii" or img_name[-7:] == ".nii.gz":
            # Convert from NIFTI to Analyze
            file_basename = img_name[:-4] if img_name[-4:] == ".nii" else img_name[:-7]
            out_name = file_basename + ".hdr"

            Utils.save_nifti(img, data, out_name)

            print(
                f"File {img_name} was converted to Analyze format and saved as {out_name}"
            )

            return out_name

        else:
            raise TypeError(
                dedent(
                    """Unsupported file format.
                    Please provide a NIFTI (.nii, .nii.gz)
                    or Analyze (.hdr, .img) file."""
                )
            )

    @staticmethod
    def dcm_nii_dcm2niix(input_folder, output_folder, output_filename):
        """
        Converts a dicom folder into a nifti using dcm2niix
        dcm2niix should be installed and added to the PATH environment variable
        :param input_folder: Input DICOM folder
        :param output_folder: Where the Nifti images will go
        :param output_filename: Naming of the output imaegs
        :return:
        """

        # Create the output folder if it doesn't exist
        os.makedirs(output_folder, exist_ok=True)

        # Run dcm2niix command to convert DICOM to NIfTI
        command = f"dcm2niix -o {output_folder} -f {output_filename} {input_folder}"
        subprocess.run(command, shell=True)


class Analysis:
    """
    Basic utilities for image analysis (statistics, etc...)
    """

    @staticmethod
    def create_mean_std_imgs(images, output_mean, output_std):
        """This function creates the mean and standard deviation images for a list of images,
        and saves them in nifti files.

        :param images: a list of strings (nifti files paths)
        :type images: list
        :param output_mean: the name of the nifti file that will contain the mean image
        :type output_mean: str
        :param output_std: the name of the nifti file that will contain the standard deviation image
        :type output_std: str
        :return: None
        :rtype: None
        """

        sample_nii = images[0]
        sample_img, sample_data = Utils.load_nifti(sample_nii)

        mean_data = sample_data * 0
        std_data = sample_data * 0

        data = sample_data[:, :, :, np.newaxis]

        for i in range(1, len(images)):
            img, img_data = Utils.load_nifti(images[i])
            img_data = img_data[:, :, :, np.newaxis]
            data = np.append(data, img_data, axis=3)

        for i in range(sample_data.shape[0]):
            for j in range(sample_data.shape[1]):
                for k in range(sample_data.shape[2]):
                    values = data[i, j, k, :]
                    mean_val = np.mean(values)
                    std_val = np.std(values)
                    mean_data[i, j, k] = mean_val
                    std_data[i, j, k] = std_val

        # print("Finished %s of %s slices." % (i, sample_data.shape[0]))

        Utils.save_nifti(sample_img, mean_data, output_mean)
        Utils.save_nifti(sample_img, std_data, output_std)

    @staticmethod
    def transform_img_to_voxel_zscores(img_, mean, std, out):
        """
        Transforms an image to voxel-by-voxel z-scores using normative data
        :param img_: Input image
        :param mean: Mean image from a normative group
        :param std: Standard deviation image from a normative group
        :param out: Z-scores image (Nifti, Analyze)
        :return:
        """
        if exists(img_) and exists(mean) and exists(std):
            img, data = Utils.load_nifti(img_)
            mean_img, mean_data = Utils.load_nifti(mean)
            std_img, std_data = Utils.load_nifti(std)

            z_scores_data = (data - mean_data) / std_data
            Utils.save_nifti(img, z_scores_data, out)

        else:
            raise FileNotFoundError("File not found at")

    @staticmethod
    def transform_img_to_atlas_zscores(img_, out, atlas_csv, atlas_hdr):
        """This function transforms an image to atlas z-scores, based on an atlas csv and hdr file,
        and saves the output image in a nifti file.

        :param img_: the name of the nifti file that contains the image to be transformed
        :type img_: str
        :param out: the name of the nifti file that will contain the output image
        :type out: str
        :param atlas_csv: CSV containing the atlas information (ROI_NUM, MEAN, STD)
        :type atlas_csv: str
        :param atlas_hdr: the name of the hdr file that contains the atlas image
        :type atlas_hdr: str
        :return: None
        :rtype: None
        """

        atlas_df = pd.read_csv(atlas_csv)
        atlas_img, atlas_data = Utils.load_nifti(atlas_hdr)

        if exists(img_):
            img, data = Utils.load_nifti(img_)

            pat_atlas = np.zeros(atlas_data.shape)

            for indx, row in atlas_df.iterrows():
                roi_num = row["ROI_NUM"]
                roi_mean = row["ROI_MEAN"]
                roi_std = row["ROI_STD"]

                index_ = np.where(atlas_data == roi_num)
                value = np.mean(data[index_])
                z_score = (value - roi_mean) / roi_std
                pat_atlas[index_] = z_score

            Utils.save_nifti(atlas_img, pat_atlas, out)

            print(f"Done transforming {img} to atlas z-scores!")

        else:
            raise FileNotFoundError("File not found at: " + img_)

    @staticmethod
    def create_atlas_csv_from_normals_imgs(normals, output_csv, atlas_csv, atlas_hdr):
        """This function creates a csv file with the mean and
        standard deviation of the ROI values for a list of
        images, based on an atlas csv and hdr file.

        :param normals: a list of strings (nifti files paths)
        :type normals: list
        :param output_csv: the name of the csv file that will contain the output data
        :type output_csv: str
        :param atlas_csv:the CSV file that contains the atlas information (ROI_NUM, ROI_NAME)
        :type atlas_csv: str
        :param atlas_hdr: the name of the hdr file that contains the atlas image
        :type atlas_hdr: str
        :return: None
        :rtype: None
        """

        atlas_df = pd.read_csv(atlas_csv, sep=";")
        atlas_img, atlas_data = Utils.load_nifti(atlas_hdr)

        for indx_, row_ in atlas_df.iterrows():
            roi_num = row_["ROI_NUM"]
            roi_name = row_["ROI_NAME"]
            index_ = np.where(atlas_data == roi_num)

            roi_values = []

            for img_ in normals:
                img, img_data = Utils.load_nifti(img_)
                value = np.mean(img_data[index_])
                roi_values.append(value)

            roi_mean = np.mean(roi_values)
            roi_std = np.std(roi_values)

            atlas_df.loc[indx_, "ROI_MEAN"] = roi_mean
            atlas_df.loc[indx_, "ROI_STD"] = roi_std

            print(f"{roi_name}: {roi_mean} +/- {roi_std}")

        atlas_df.to_csv(output_csv)
        print("Done!")

    @staticmethod
    def calculate_z_scores_array(image, atlas):
        """Calculates the z-scores for each ROI in an atlas.
        :param image: the path to the input image
        :type image: str
        :param atlas: the path to the atlas image
        :type atlas: str
        :return: None
        :rtype: None
        """

        img, img_data = Utils.load_nifti(image)
        atlas, atlas_data = Utils.load_nifti(atlas)

        values = []

        rois = np.unique(atlas_data)
        for i in rois:
            indx = np.where(atlas_data == i)
            mean = np.mean(img_data[indx])
            values.append(mean)

        return values

    @staticmethod
    def voxel_wise_partial_pearson_images_scale(
        images, scale, covariate, mask, output_rs
    ):
        """Calculates the partial pearson correlation coefficient between images and scalars
        (i.e. a neurpsychological scale)
        for each voxel in the mask image.
        :param images: a list of paths for the  images
        :type images: list
        :param scale: a list of scale values
        :type scale: list
        :param covariate: a list of covariate values
        :type covariate: list
        :param mask: the path to the mask image
        :type mask: str
        :param output_rs: the path to the output file for the correlation coefficients
        :type output_rs: str
        """
        mask_img, mask_data = Utils.load_nifti(mask)

        sample_img, sample_data = Utils.load_nifti(images[0])

        corr_r = sample_data * 0

        data = sample_data[..., np.newaxis]

        for image in images[1:]:
            img, img_data = Utils.load_nifti(image)
            img_data = img_data[..., np.newaxis]
            data = np.append(data, img_data, axis=3)

        print("Data shape: ", data.shape)
        print("Mask shape: ", mask_data.shape)
        print("Scale values: ", len(scale))

        for i in range(sample_data.shape[0]):
            for j in range(sample_data.shape[1]):
                for k in range(sample_data.shape[2]):
                    if mask_data[i, j, k] == 0:
                        pass

                    else:
                        data_array = data[i, j, k, :]

                        df = pd.DataFrame()
                        df["var1"] = pd.Series(data_array)
                        df["var2"] = pd.Series(scale)
                        df["cov"] = pd.Series(covariate)

                        res = pg.partial_corr(data=df, x="var1", y="var2", covar="cov")

                        corr_r[i, j, k] = res.r

        Utils.save_nifti(sample_img, corr_r, output_rs)

        print("Correlation Analysis Finished")

    @staticmethod
    def voxel_wise_pearson_images_scale(
        images, scale, mask, output_rs, output_ps, corr="pearson"
    ):
        """Calculates the pearson correlation coefficient and p-value between images and scalars (i.e. a neurpsychological scale)
        for each voxel in the mask image.
        :param images: a list of paths for the images
        :type images: list
        :param scale: a list of scale values
        :type scale: list
        :param mask: the path to the mask image
        :type mask: str
        :param output_rs: the path to the output file for the correlation coefficients
        :type output_rs: str
        :param output_ps: the path to the output file for the p-values
        :type output_ps: str
        :param corr: the type of correlation coefficient to calculate
        :type corr: str
        """

        mask_data = Utils.load_nifti(mask, only_data=True)
        sample_img, sample_data = Utils.load_nifti(images[0])

        corr_r = sample_data * 0
        corr_p = sample_data * 0

        data = sample_data[..., np.newaxis]

        for image in images[1:]:
            img_data = Utils.load_nifti(image, only_data=True)
            img_data = img_data[..., np.newaxis]
            data = np.append(data, img_data, axis=3)

        print("Data shape: ", data.shape)
        print("Mask shape: ", mask_data.shape)
        print("Scale values: ", len(scale))

        for i in range(sample_data.shape[0]):
            for j in range(sample_data.shape[1]):
                for k in range(sample_data.shape[2]):
                    if mask_data[i, j, k] == 0:
                        pass

                    else:
                        data_array = data[i, j, k, :]

                        if corr == "pearson":
                            r_, p_ = pearsonr(data_array, scale)

                        elif corr == "spearman":
                            r_, p_ = spearmanr(data_array, scale)

                        else:
                            raise ValueError(
                                "corr variable must be pearson or spearman"
                            )

                        corr_r[i, j, k] = r_
                        corr_p[i, j, k] = p_

        Utils.save_nifti(sample_img, corr_r, output_rs)
        Utils.save_nifti(sample_img, corr_p, output_ps)

        print("Correlation Analysis Finished")

    @staticmethod
    def image_to_image_corr_atlas_based_spearsman(image_1, image_2, atlas):
        """Calculates the spearman correlation coefficient and p-value between two images using
        the atlas ROI-values extracted for each region of interest (ROI) defined by an atlas.

        :param image_1: the path to the first image
        :type image_1: str
        :param image_2: the path to the second image
        :type image_2: str
        :param atlas: the path to the atlas image that defines the ROIs
        :type atlas: str
        :return: the spearman correlation coefficient and the p-value for each ROI
        :rtype: float, float
        """

        # Load the two images + atlas using nibabel library
        img_1_data = Utils.load_nifti(image_1, only_data=True)
        img_2_data = Utils.load_nifti(image_2, only_data=True)
        atlas_data = Utils.load_nifti(atlas, only_data=True)

        # Get the unique values of the atlas, which correspond to the ROIs
        rois = np.unique(atlas_data)

        # Initialize lists to store the mean values of the images for each ROI
        img_1_array = []
        img_2_array = []

        # Loop over the ROIs
        for i in rois:
            if i != 0:
                # Get the indices of the voxels that belong to the current ROI
                indx = np.where(atlas_data == i)

                # Calculate the mean value of the first image for the current ROI
                mean_1 = np.mean(img_1_data[indx])
                # Calculate the mean value of the second image for the current ROI
                mean_2 = np.mean(img_2_data[indx])

                # Append the mean values to the lists
                img_1_array.append(mean_1)
                img_2_array.append(mean_2)

        # Calculate the spearman correlation coefficient and p-value between the mean values of the images for each ROI
        rho, p = spearmanr(img_1_array, img_2_array)
        return rho, p

    @staticmethod
    def ncc(img1_path, img2_path):
        """Calculates the normalized cross correlation (NCC) between two images.

        :param img1_path: the path to the first image
        :type img1_path: str
        :param img2_path: the path to the second image
        :type img2_path: str
        :return: the NCC value between the two images
        :rtype: float
        """

        # Load Nifti images
        img1_data = Utils.load_nifti(img1_path, only_data=True)
        img2_data = Utils.load_nifti(img2_path, only_data=True)

        # Get image data as numpy arrays

        # Subtract the mean from each image
        img1_data = img1_data - np.mean(img1_data)
        img2_data = img2_data - np.mean(img2_data)

        # Calculate the numerator of the NCC equation
        numerator = np.sum(img1_data * img2_data)

        # Calculate the denominator of the NCC equation
        denominator = np.sqrt(np.sum(img1_data**2) * np.sum(img2_data**2))

        # Calculate and return the NCC
        return numerator / denominator

    @staticmethod
    def spm_map_2_cohens_d(img, out, len_1, len_2):
        """
        Converts an image to cohens_d
        Expected input is a spmT_0001.nii file
        """
        img, data = Utils.load_nifti(img)
        d_coeff = np.sqrt(1 / len_1 + 1 / len_2)
        data = data * d_coeff

        Utils.save_nifti(img, data, out)

    @staticmethod
    def get_tvalue_thresholds_FDR(img_, n1, n2):
        # Load the NIFTI image
        data = Utils.load_nifti(img_, only_data=True)

        # Flatten the data to get a 1D array of t-values
        t_values = data.flatten()
        t_values = abs(np.unique(t_values[t_values != 0]))
        t_values = np.sort(t_values)

        df = n1 + n2 - 2
        p_values = t.sf(abs(t_values), df=df)
        indx = np.where(p_values < 0.05)
        thresholded = t_values[indx]
        thres = np.percentile(thresholded, 5)

        d_coeff = np.sqrt(1 / n1 + 1 / n2)
        cohens_thres = thres * d_coeff

        print(cohens_thres)
        new_file_name = join(dirname(img_), "cohensd_thres.txt")
        with open(new_file_name, "w") as new_file:
            new_file.write(str(cohens_thres))
        new_file.close()


class Preprocessing:
    @staticmethod
    def normalize_histogram(input_image, template, mask, output):
        """Normalizes an image using the mode of an intensity histogram.
        More info at: https://pubmed.ncbi.nlm.nih.gov/32771619/

        :param input_image: the path to the input image
        :type input_image: str
        :param template: the path to the template image
        :type template: str
        :param mask: the path to the mask image
        :type mask: str
        :param output: the path to the output image
        :type output: str
        :return: the normalization value used to scale the input image
        :rtype: float
        """

        fdg, fdg_data = Utils.load_nifti(input_image)
        template, template_data = Utils.load_nifti(template)
        mask, mask_data = Utils.load_nifti(mask)

        if len(mask.shape) == 4:
            mask_data = mask_data[:, :, :, 0]

        indx = np.where(mask_data == 1)
        mean_template = np.mean(template_data[indx])
        mean_fdg = np.mean(fdg_data[indx])

        fdg_data = fdg_data * (mean_template / mean_fdg)

        division = template_data[indx] / fdg_data[indx]
        values, bins = np.histogram(division, 200, range=[0.5, 2])
        amax = np.amax(values)
        indx = np.where(values == amax)
        norm_value = bins[indx][0]
        fdg_data = fdg_data * norm_value

        Utils.save_nifti(fdg, fdg_data, output)

        return norm_value

    @staticmethod
    def normalize_using_ref_region(input_image, output_image, ref_region):
        """Normalizes an image using a reference region.

        :param input_image: the path to the input image
        :type input_image: str
        :param output_image: the path to the output image
        :type output_image: str
        :param ref_region: the path to the reference region image
        :type ref_region: str
        :return: None
        :rtype: None
        """
        pons_img = Utils.load_nifti(ref_region, only_data=True)

        if len(pons_img.shape) == 4:
            pons_img = pons_img[:, :, :, 0]
        pons_vox = np.where(pons_img == 1)

        input_img, img_data = Utils.load_nifti(input_image)
        if len(input_img.shape) == 4:
            img_data = img_data[:, :, :, 0]

        pons_value = np.mean(img_data[pons_vox])
        normalized_img = img_data / pons_value

        Utils.save_nifti(input_img, normalized_img, output_image)

    @staticmethod
    def histogram_matching(reference_nii, input_nii, output_nii):
        """Matches the histogram of an input image to a reference image.

        :param reference_nii: the path to the reference image
        :type reference_nii: str
        :param input_nii: the path to the input image
        :type input_nii: str
        :param output_nii: the path to the output image
        :type output_nii: str
        :return: None
        :rtype: None
        """

        nt_data = Utils.load_nifti(reference_nii, only_data=True)
        patient, pt_data = Utils.load_nifti(input_nii)

        # Stores the image data shape that will be used later
        oldshape = pt_data.shape

        # Converts the data arrays to single dimension and normalizes by the maximum
        nt_data_array = nt_data.ravel()
        pt_data_array = pt_data.ravel()

        # get the set of unique pixel values and their corresponding indices and counts
        s_values, bin_idx, s_counts = np.unique(
            pt_data_array, return_inverse=True, return_counts=True
        )
        t_values, t_counts = np.unique(nt_data_array, return_counts=True)

        # take the cumsum of the counts and normalize by the number of pixels to
        # get the empirical cumulative distribution functions for the source and
        # template images (maps pixel value --> quantile)
        s_quantiles = np.cumsum(s_counts).astype(np.float64)
        s_quantiles /= s_quantiles[-1]
        t_quantiles = np.cumsum(t_counts).astype(np.float64)
        t_quantiles /= t_quantiles[-1]

        # interpolate linearly to find the pixel values in the template image
        # that correspond most closely to the quantiles in the source image
        interp_t_values = np.interp(s_quantiles, t_quantiles, t_values)

        # Reshapes the corresponding values to the indexes and reshapes the array to input
        final_image_data = interp_t_values[bin_idx].reshape(oldshape)
        # final_image_data[indx] = 0

        Utils.save_nifti(patient, final_image_data, output_nii)

    @staticmethod
    def logpow_histogram_matching(
        reference_nii, input_nii, output_nii, alpha=1, beta=3
    ):
        """Matches the histogram of an input image to a reference image using a log-power transformation.
        More info: https://doi.org/10.1117/1.JEI.23.6.063017

        :param reference_nii: the path to the reference image
        :type reference_nii: str
        :param input_nii: the path to the input image
        :type input_nii: str
        :param output_nii: the path to the output image
        :type output_nii: str
        :param alpha: the additive constant for the log transformation, defaults to 1
        :type alpha: int, optional
        :param beta: the power exponent for the log transformation, defaults to 3
        :type beta: int, optional
        :return: the path to the output image
        :rtype: str
        """
        nt_data = Utils.load_nifti(reference_nii, only_data=True)
        patient, pt_data = Utils.load_nifti(input_nii)

        # Stores the image data shape that will be used later
        oldshape = pt_data.shape

        # Converts the data arrays to single dimension and normalizes by the maximum
        nt_data_array = nt_data.ravel()
        pt_data_array = pt_data.ravel()

        # get the set of unique pixel values and their corresponding indices and counts
        s_values, bin_idx, s_counts = np.unique(
            pt_data_array, return_inverse=True, return_counts=True
        )
        t_values, t_counts = np.unique(nt_data_array, return_counts=True)

        s_counts = np.power(np.log10(s_counts + alpha), beta)
        t_counts = np.power(np.log10(t_counts + alpha), beta)

        # take the cumsum of the counts and normalize by the number of pixels to
        # get the empirical cumulative distribution functions for the source and
        # template images (maps pixel value --> quantile)
        s_quantiles = np.cumsum(s_counts).astype(np.float64)
        s_quantiles /= s_quantiles[-1]
        t_quantiles = np.cumsum(t_counts).astype(np.float64)
        t_quantiles /= t_quantiles[-1]

        # interpolate linearly to find the pixel values in the template image
        # that correspond most closely to the quantiles in the source image
        interp_t_values = np.interp(s_quantiles, t_quantiles, t_values)

        # Reshapes the corresponding values to the indexes and reshapes the array to input
        final_image_data = interp_t_values[bin_idx].reshape(oldshape)
        # final_image_data[indx] = 0

        Utils.save_nifti(patient, final_image_data, output_nii)

    @staticmethod
    def estimate_fwhm_mizutani(
        nifti, bin_size=5, n_segs=2, orientation="axial", plot=False
    ):
        """
        This function estimates image xy resolution based in
        a previous method implemented by MIZUTANI ET AL. for microscopy imaging:

        https://onlinelibrary.wiley.com/doi/10.1111/jmi.12315

        The idea of expanding this method to PET imaging was proposed by Carbonell et al:
        https://n.neurology.org/content/94/15_Supplement/4301

        The FWHMs are estimated from multiple regression of the logarithm of the square norm of
        the image Fourier transform against the square distance from the origin in the Fourier
        domain(essentially from what is known as a Wilson plot).
        The lower frequencies in the Wilson plot are linearly fitted
        and the PSF is estimated from the slope.

        :param orientation: plane to walk traverse the image (axial, sagittal, coronal).
        :param nifti: input nifti image

        :param bin_size: In the original paper from Mizutani,
        the logarithm of the average squared norm in each
        5 × 5 pixel bin of the Fourier transform was used.
        You can vary this to check what fits your data better.
        Values between 5 and 10 seem to work well for PET images.

        :param plot: If True, the function will plot the Wilson
        plots and fitted lines for quality check.

        :return: The average FWHM across all slices in the image.
        """

        img, volume = Utils.load_nifti(nifti)
        image_xy_size = np.abs(img.affine[0, 0])
        fwhm_values = []
        sigma_values = []

        if orientation == "axial":
            orient = 2
        elif orientation == "coronal":
            orient = 1
        elif orientation == "sagittal":
            orient = 0
        else:
            raise ValueError("orientation must be set to axial, coronal or sagital")

        for i in range(volume.shape[orient]):
            if orientation == "axial":
                slice_2d = volume[:, :, i]
            elif orientation == "coronal":
                slice_2d = volume[:, i, :]
            elif orientation == "sagittal":
                slice_2d = volume[i, :, :]

            # Compute the 2D Fourier transform of the input slice
            fourier_slice = fftshift(fftn(slice_2d))

            # Compute the square norm of the Fourier transform
            square_norm = np.abs(fourier_slice) ** 2

            # Calculate logarithm of the square norm
            log_square_norm = np.log(square_norm + np.finfo(float).eps)

            # Calculate k^2 values (distances from the origin in Fourier domain)
            k_values = np.array(
                np.meshgrid(
                    np.fft.fftshift(np.fft.fftfreq(slice_2d.shape[0])),
                    np.fft.fftshift(np.fft.fftfreq(slice_2d.shape[1])),
                    indexing="ij",
                )
            )
            k_squared = np.sum(k_values**2, axis=0)

            # Calculate the number of bins along each axis
            num_bins_x = slice_2d.shape[0] // bin_size
            num_bins_y = slice_2d.shape[1] // bin_size

            # Initialize binned arrays for the logarithm of the average squared norm and k^2 values
            log_square_norm_binned = np.zeros((num_bins_x, num_bins_y))
            k_squared_binned = np.zeros((num_bins_x, num_bins_y))

            for x in range(num_bins_x):
                for y in range(num_bins_y):
                    # Compute the average of log_square_norm and k_squared within the 5x5 bin
                    log_square_norm_binned[x, y] = np.mean(
                        log_square_norm[
                            x * bin_size : (x + 1) * bin_size,
                            y * bin_size : (y + 1) * bin_size,
                        ]
                    )
                    k_squared_binned[x, y] = np.mean(
                        k_squared[
                            x * bin_size : (x + 1) * bin_size,
                            y * bin_size : (y + 1) * bin_size,
                        ]
                    )

            # Flatten the binned arrays and remove zero elements
            k_squared_flat = k_squared_binned.flatten()
            log_square_norm_flat = log_square_norm_binned.flatten()

            non_zero_mask = k_squared_flat > 0
            k_squared_flat = k_squared_flat[non_zero_mask]
            log_square_norm_flat = log_square_norm_flat[non_zero_mask]

            # Instantiate the piecewise linear regression model
            model = pwlf.PiecewiseLinFit(k_squared_flat, log_square_norm_flat)

            # Determine the optimal number of line segments
            n_segments = n_segs
            breakpoints = model.fit(n_segments)

            # Predict the fitted line
            x_pred = np.linspace(k_squared_flat.min(), k_squared_flat.max(), 100)
            y_pred = model.predict(x_pred)

            # Extract the slope and intercept of the fitted line
            slope, intercept = model.slopes[0], model.intercepts[0]

            # Calculate the width σ of the standard PSF
            sigma = np.sqrt(-slope / (4 * np.pi**2))
            sigma_values.append(sigma)

            # Calculate the FWHM of the PSF
            fwhm = 2 * np.sqrt(2 * np.log(2)) * sigma * image_xy_size
            fwhm_values.append(fwhm)

            if i == volume.shape[2] // 3 and plot:
                # Plot the Wilson plot with the fitted line overlaid
                plt.scatter(k_squared_flat, log_square_norm_flat, s=1)
                plt.plot(x_pred, y_pred, color="black", linewidth=1)
                plt.xlabel("|k|^2")
                plt.ylabel("ln|F(k)|^2")
                plt.title("Wilson plot with fitted line")
                plt.ylim(
                    [np.min(log_square_norm_flat) - 5, np.max(log_square_norm_flat) + 5]
                )
                plt.show()

        return np.nanmean(fwhm_values), np.nanmean(sigma_values)


class Utils:
    """Small Utilities for imaging work"""

    @staticmethod
    def load_nifti(img_path, only_data=False):
        """Gets image path, returns image object and data"""
        img = nib.load(img_path)
        data = img.get_fdata()

        if only_data:
            return data
        return img, data

    @staticmethod
    def save_nifti(template_img, data, out_path):
        """Gets template image object for header, data and outpath
        Saves the image"""
        img = nib.Nifti1Image(data, template_img.affine, template_img.header)
        nib.save(img, out_path)

    @staticmethod
    def _check_input_image_shape(img_data):
        if img_data.ndim != 3:
            img_data = img_data[:, :, :, 0]

        return img_data

    @staticmethod
    def change_image_dtype(input_filepath, output_filepath, new_dtype=np.float32):
        """
        Options for new_dtype
        np.uint8: 8-bit unsigned integer
        np.uint16: 16-bit unsigned integer
        np.uint32: 32-bit unsigned integer
        np.int8: 8-bit signed integer
        np.int16: 16-bit signed integer
        np.int32: 32-bit signed integer
        np.float16: 16-bit floating-point
        np.float32: 32-bit floating-point (single precision)
        np.float64: 64-bit floating-point (double precision)
        """

        # Load NIfTI or Analyze image
        img, img_data = Utils.load_nifti(input_filepath)

        # Change the data type of the image data
        new_img_data = img_data.astype(new_dtype)

        # Create a new NIfTI or Analyze image with the changed data type
        if isinstance(img, nib.nifti1.Nifti1Image):
            new_img = nib.Nifti1Image(new_img_data, img.affine, header=img.header)
        elif isinstance(img, nib.analyze.AnalyzeImage):
            header = img.header.copy()
            header.set_data_dtype(new_dtype)
            new_img = nib.analyze.AnalyzeImage(new_img_data, img.affine, header)

        # Save the image with the changed data type to the specified output filepath
        nib.save(new_img, output_filepath)

    @staticmethod
    def resample_image_by_matrix_size(
        input_filepath, output_filepath, target_shape, interpolation="linear"
    ):
        """
        # Example usage
        input_image_file = '/path/to/input/image/file.nii.gz'
        output_resampled_image_file = '/path/to/output/resampled/image/file.nii.gz'
        target_shape = (128, 128, 64)
        interpolation_method = 'quadratic'  # Specify the desired interpolation method
        resample_image(input_image_file, output_resampled_image_file, target_shape, interpolation_method)

        Interpolation options: nearest, linear, cubic, quadratic
        """

        # Load NIfTI or Analyze image
        img, img_data = Utils.load_nifti(input_filepath)
        if len(img_data.shape) == 3:
            target_shape = target_shape[0:3]

        # Calculate zoom factors for each dimension
        zoom_factors = np.array(target_shape) / np.array(img_data.shape)

        # Determine the order of interpolation based on the given method
        interpolation_methods = {
            "nearest": 0,
            "linear": 1,
            "cubic": 3,
            "quadratic": 4,
            # Add more interpolation methods and their corresponding order values here
        }
        order = interpolation_methods.get(
            interpolation, 1
        )  # Default to linear if method not recognized

        # Resample the image data using scipy's zoom function
        resampled_img_data = zoom(img_data, zoom_factors, order=order)

        # Update the affine matrix to reflect the new voxel dimensions
        voxel_sizes = np.array(img.header.get_zooms()) * (
            np.array(img_data.shape) / target_shape
        )
        new_affine = img.affine.copy()
        np.fill_diagonal(new_affine, voxel_sizes)

        # Create a new NIfTI or Analyze image with the resampled data and updated affine
        if isinstance(img, nib.nifti1.Nifti1Image):
            resampled_img = nib.Nifti1Image(
                resampled_img_data, new_affine, header=img.header
            )
        elif isinstance(img, nib.analyze.AnalyzeImage):
            header = img.header.copy()
            header.set_data_shape(target_shape)
            resampled_img = nib.analyze.AnalyzeImage(
                resampled_img_data, new_affine, header
            )

        # Save the resampled image to the specified output filepath
        nib.save(resampled_img, output_filepath)

    @staticmethod
    def resample_image_by_voxel_sizes(
        input_filepath, output_filepath, target_voxel_sizes, interpolation="linear"
    ):
        """
        # Example usage
        input_image_file = '/path/to/input/image/file.nii.gz'
        output_resampled_image_file = '/path/to/output/resampled/image/file.nii.gz'
        target_voxel_sizes = (2.0, 2.0, 2.0)  # Specify the desired voxel sizes
        interpolation_method = 'linear'  # Specify the desired interpolation method

        resample_image_with_voxel_sizes(input_image_file, output_resampled_image_file, target_voxel_sizes, interpolation_method)
        """

        # Load NIfTI or Analyze image
        img = nib.load(input_filepath)
        img_data = img.get_fdata()

        # Calculate zoom factors based on target voxel sizes
        current_voxel_sizes = np.array(img.header.get_zooms())
        zoom_factors = current_voxel_sizes / target_voxel_sizes

        # Determine the order of interpolation based on the given method
        interpolation_methods = {
            "nearest": 0,
            "linear": 1,
            "cubic": 3,
            "quadratic": 4,
            # Add more interpolation methods and their corresponding order values here
        }
        order = interpolation_methods.get(
            interpolation, 1
        )  # Default to linear if method not recognized

        # Resample the image data using scipy's zoom function
        resampled_img_data = zoom(img_data, zoom_factors, order=order)

        # Update the affine matrix to reflect the new voxel sizes
        new_affine = img.affine.copy()
        np.fill_diagonal(new_affine, target_voxel_sizes)

        # Create a new NIfTI or Analyze image with the resampled data and updated affine
        if isinstance(img, nib.nifti1.Nifti1Image):
            resampled_img = nib.Nifti1Image(
                resampled_img_data, new_affine, header=img.header
            )
        elif isinstance(img, nib.analyze.AnalyzeImage):
            header = img.header.copy()
            resampled_img = nib.analyze.AnalyzeImage(
                resampled_img_data, new_affine, header
            )

        # Save the resampled image to the specified output filepath
        nib.save(resampled_img, output_filepath)

    @staticmethod
    def remove_nan_negs(input_filepath, output_filepath):
        """
        Removes nan and neg values from image
        :param input_filepath: Input image filepath as Nifti or Analyze
        :param output_filepath: Output image filepath as Nifti or Analyze
        :return:
        """

        # Load NIfTI or Analyze image
        img = nib.load(input_filepath)
        img_data = img.get_fdata()

        # Replace NaN and negative values with zeros
        img_data[np.isnan(img_data) | (img_data < 0)] = 0

        # Create a new NIfTI or Analyze image with cleaned data
        if isinstance(img, nib.nifti1.Nifti1Image):
            cleaned_img = nib.Nifti1Image(img_data, img.affine, header=img.header)
        elif isinstance(img, nib.analyze.AnalyzeImage):
            header = img.header.copy()
            cleaned_img = nib.analyze.AnalyzeImage(img_data, img.affine, header)

        # Save the cleaned image to the specified output filepath
        nib.save(cleaned_img, output_filepath)

    @staticmethod
    def reorient_and_clean(img_path, out_path):
        img = nib.load(img_path)

        # Reorient to closest canonical orientation
        img_canonical = nib.as_closest_canonical(img)

        # Get image data as numpy array
        data = img_canonical.get_fdata()

        # Set NaN values to 0
        data[np.isnan(data)] = 0

        # Set negative values to 0
        data[data < 0] = 0

        # Create a new NIfTI image with the modified data and save it
        new_img = nib.Nifti1Image(data, img_canonical.affine, img_canonical.header)
        nib.save(new_img, out_path)

    @staticmethod
    def add_poisson_noise(
        input_filepath, output_filepath, intensity_scaling_factor=1.0
    ):
        """Adds poison noise to an image"""
        # Load NIfTI or Analyze image
        img = nib.load(input_filepath)
        img_data = img.get_fdata()

        # Apply Poisson noise to the image
        noisy_img_data = np.random.poisson(img_data * intensity_scaling_factor)

        # Create a new NIfTI or Analyze image with the noisy data
        if isinstance(img, nib.nifti1.Nifti1Image):
            noisy_img = nib.Nifti1Image(noisy_img_data, img.affine, img.header)
            nib.save(noisy_img, output_filepath)

        elif isinstance(img, nib.analyze.AnalyzeImage):
            noisy_img = nib.analyze.AnalyzeImage(noisy_img_data, img.affine, img.header)
            nib.save(noisy_img, output_filepath)

        # Save the noisy image to the specified output filepath

    @staticmethod
    def apply_constant_to_img(image, c, operation, output=False):
        """This function applies a constant to an image using a
        specified operation, and saves the result in a new image.

        :param image: the name of the image file to be processed
        :type image: str
        :param c: the constant to be applied to the image
        :type c: float
        :param operation: Operation to be performed on the image. It can be 'mult', 'div' or 'sum'
        :type operation: str
        :param output: If False, the output file will have the same name as the input file
        :type output: str or bool
        :return: None
        """
        if output:
            output_file = output

        else:
            output_file = image

        img_, data = Utils.load_nifti(image)

        if operation == "mult":
            data = data * c

        elif operation == "div":
            data = data / c

        elif operation == "sum":
            data = data + c

        Utils.save_nifti(img_, data, output_file)
