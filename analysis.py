import os
import numpy as np
import nibabel as nib
from os.path import exists
import pandas as pd
from scipy import signal
from qc_utils import qc_utils
import pingouin as pg


def create_mean_std_imgs(images, output_mean, output_std):
    """This function creates the mean and standard deviation images for a list of images,
    and saves them in nifti files.

    :param images: a list of strings containing the names of the nifti files that contain the images to be processed
    :type images: list
    :param output_mean: the name of the nifti file that will contain the mean image
    :type output_mean: str
    :param output_std: the name of the nifti file that will contain the standard deviation image
    :type output_std: str
    :return: None
    :rtype: None
    """

    sample_nii = images[0]
    sample_img = nib.load(sample_nii)
    sample_data = sample_img.get_fdata()

    mean_data = sample_data * 0
    std_data = sample_data * 0

    data = sample_data[:, :, :, np.newaxis]

    for i in range(1, len(images)):
        img = nib.load(images[i])
        img_data = img.get_fdata()
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

    mean_img = nib.AnalyzeImage(mean_data, sample_img.affine, sample_img.header)
    nib.save(mean_img, output_mean)
    std_img = nib.AnalyzeImage(std_data, sample_img.affine, sample_img.header)
    nib.save(std_img, output_std)


def transform_img_to_voxel_zscores(img_, mean, std, out):

    if exists(img_) and exists(mean) and exists(std):

        img = nib.load(img_)
        data = img.get_fdata()

        mean_img = nib.load(mean)
        mean_data = mean_img.get_fdata()

        std_img = nib.load(mean)
        std_data = std_img.get_fdata()

        z_scores_data = (data - mean_data) / std_data
        img_out = nib.Nifti1Image(z_scores_data, img.affine, img.header)
        nib.save(img_out, out)

    else:
        raise FileNotFoundError("File not found at")


def transform_img_to_atlas_zscores(img_, out, atlas_csv, atlas_hdr):
    """This function transforms an image to atlas z-scores, based on an atlas csv and hdr file,
    and saves the output image in a nifti file.

    :param img_: the name of the nifti file that contains the image to be transformed
    :type img_: str
    :param out: the name of the nifti file that will contain the output image
    :type out: str
    :param atlas_csv: the name of the csv file that contains the atlas information, such as ROI_NUM, ROI_MEAN and ROI_STD
    :type atlas_csv: str
    :param atlas_hdr: the name of the hdr file that contains the atlas image
    :type atlas_hdr: str
    :return: None
    :rtype: None
    """

    atlas_csv = atlas_csv
    atlas_df = pd.read_csv(atlas_csv)
    atlas_img = nib.load(atlas_hdr)
    atlas_data = atlas_img.get_fdata()

    if exists(img_):

        img = nib.load(img_)
        data = img.get_fdata()

        pat_atlas = np.zeros(atlas_data.shape)

        for indx, row in atlas_df.iterrows():
            roi_num = row['ROI_NUM']
            roi_mean = row['ROI_MEAN']
            roi_std = row['ROI_STD']

            index_ = np.where(atlas_data == roi_num)
            value = np.mean(data[index_])
            z_score = (value - roi_mean) / roi_std
            pat_atlas[index_] = z_score

        img = nib.Nifti1Image(pat_atlas, atlas_img.affine, atlas_img.header)
        nib.save(img, out)

        print('Done transforming %s to atlas z-scores!' % img_)

    else:
        raise FileNotFoundError("File not found at: " + img_)


def create_atlas_csv_from_normals_imgs(normals, output_csv, atlas_csv, atlas_hdr):
    """This function creates a csv file with the mean and standard deviation of the ROI values for a list of
    images, based on an atlas csv and hdr file.

    :param normals: a list of strings containing the names of the nifti files that contain the images to be processed
    :type normals: list
    :param output_csv: the name of the csv file that will contain the output data
    :type output_csv: str
    :param atlas_csv: the name of the csv file that contains the atlas information, such as ROI_NUM and ROI_NAME
    :type atlas_csv: str
    :param atlas_hdr: the name of the hdr file that contains the atlas image
    :type atlas_hdr: str
    :return: None
    :rtype: None
    """

    atlas_df = pd.read_csv(atlas_csv,sep=';')
    atlas_img = nib.load(atlas_hdr)
    atlas_data = atlas_img.get_fdata()[:, :, :, 0]

    for indx_, row_ in atlas_df.iterrows():

        roi_num = row_['ROI_NUM']
        roi_name = row_['ROI_NAME']
        index_ = np.where(atlas_data == roi_num)

        roi_values = []

        for img_ in normals:
            img_data = nib.load(img_).get_fdata()
            value = np.mean(img_data[index_])
            roi_values.append(value)

        roi_mean = np.mean(roi_values)
        roi_std = np.std(roi_values)

        atlas_df.loc[indx_, 'ROI_MEAN'] = roi_mean
        atlas_df.loc[indx_, 'ROI_STD'] = roi_std

        print('%s: %s +- %s' % (roi_name, roi_mean, roi_std))

    atlas_df.to_csv(output_csv)
    print('Done!')


def calculate_z_scores_array(image, atlas):
    """Calculates the z-scores for each ROI in an atlas.
    :param image: the path to the input image
    :type image: str
    :param atlas: the path to the atlas image
    :type atlas: str
    :return: None
    :rtype: None
    """

    img_data = nib.load(image).get_fdata()
    atlas_data = nib.load(atlas).get_fdata()

    values = []

    rois = np.unique(atlas_data)
    for i in rois:
        indx = np.where(atlas_data==i)
        mean = np.mean(img_data[indx])
        values.append(mean)

    print(values)


def voxel_wise_partial_pearson_images_scale(images, scale, covariate, mask, output_rs):
    """Calculates the partial pearson correlation coefficient between images and scalars (i.e. a neurpsychological scale)
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

    mask_data = nib.load(mask).get_fdata()

    sample_img = nib.load(images[0])
    sample_data = sample_img.get_fdata()

    corr_r = sample_data * 0

    data = sample_data[..., np.newaxis]

    for image in images[1:]:
        img_data = nib.load(image).get_fdata()
        img_data = img_data[..., np.newaxis]
        data = np.append(data, img_data, axis=3)

    print('Data shape: ', data.shape)
    print('Mask shape: ', mask_data.shape)
    print('Scale values: ', len(scale))

    for i in range(sample_data.shape[0]):
        for j in range(sample_data.shape[1]):
            for k in range(sample_data.shape[2]):

                if mask_data[i, j, k] == 0:
                    pass

                else:
                    data_array = data[i, j, k, :]

                    df = pd.DataFrame()
                    df['var1'] = pd.Series(data_array)
                    df['var2'] = pd.Series(scale)
                    df['cov'] = pd.Series(covariate)

                    res = pg.partial_corr(data=df, x='var1', y='var2', covar='cov')

                    corr_r[i, j, k] = res.r

    r_img = nib.Nifti1Image(corr_r, sample_img.affine, sample_img.header)
    nib.save(r_img, output_rs)

    print('Correlation Analysis Finished')


def voxel_wise_pearson_images_scale(images, scale, mask, output_rs, output_ps, corr='pearson'):
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

    from scipy.stats import pearsonr, spearmanr

    mask_data = nib.load(mask).get_fdata()

    sample_img = nib.load(images[0])
    sample_data = sample_img.get_fdata()

    corr_r = sample_data * 0
    corr_p = sample_data * 0

    data = sample_data[..., np.newaxis]

    for image in images[1:]:
        img_data = nib.load(image).get_fdata()
        img_data = img_data[..., np.newaxis]
        data = np.append(data, img_data, axis=3)

    print('Data shape: ', data.shape)
    print('Mask shape: ', mask_data.shape)
    print('Scale values: ', len(scale))

    for i in range(sample_data.shape[0]):
        for j in range(sample_data.shape[1]):
            for k in range(sample_data.shape[2]):

                if mask_data[i, j, k] == 0:
                    pass

                else:
                    data_array = data[i, j, k, :]

                    if corr == 'pearson':

                        r_, p_ = pearsonr(data_array, scale)

                    elif corr == 'spearman':

                        r_, p_ = spearmanr(data_array, scale)

                    else:

                        raise ValueError('corr variable must be pearson or spearman')

                    corr_r[i, j, k] = r_
                    corr_p[i, j, k] = p_

    r_img = nib.Nifti1Image(corr_r, sample_img.affine, sample_img.header)
    nib.save(r_img, output_rs)

    ps_img = nib.Nifti1Image(corr_p, sample_img.affine, sample_img.header)
    nib.save(ps_img, output_ps)

    print('Correlation Analysis Finished')


def image_to_image_corr_atlas_based_spearsman(image_1, image_2, atlas):

    """Calculates the spearman correlation coefficient and p-value between two images using
    the atlas ROI-values extracted for each region of interest (ROI) defined by an atlas.

    :param image_1: the path to the first image
    :type image_1: str
    :param image_2: the path to the second image
    :type image_2: str
    :param atlas: the path to the atlas image that defines the ROIs
    :type atlas: str
    :return: the spearman correlation coefficient and the p-value between the two images for each ROI
    :rtype: float, float
    """

    from scipy.stats import spearmanr

    # Load the two images + atlas using nibabel library
    img_1 = nib.load(image_1)
    img_1_data = img_1.get_fdata()

    img_2 = nib.load(image_2)
    img_2_data = img_2.get_fdata()

    atlas_img = nib.load(atlas)
    atlas_data = atlas_img.get_fdata()[:, :, :, 0]

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
    img1 = nib.load(img1_path)
    img2 = nib.load(img2_path)

    # Get image data as numpy arrays
    img1_data = img1.get_fdata()
    img2_data = img2.get_fdata()

    # Subtract the mean from each image
    img1_data = img1_data - np.mean(img1_data)
    img2_data = img2_data - np.mean(img2_data)

    # Calculate the numerator of the NCC equation
    numerator = np.sum(img1_data * img2_data)

    # Calculate the denominator of the NCC equation
    denominator = np.sqrt(np.sum(img1_data ** 2) * np.sum(img2_data ** 2))

    # Calculate and return the NCC
    return numerator / denominator





