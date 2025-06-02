import numpy as np
import nibabel as nib
from scipy.stats import pearsonr, spearmanr

def voxel_wise_corr_images_vs_scale(images, scale, mask, corr = "pearson"):
    """Calculates the pearson correlation coefficient and p-value between images and scalars
    (i.e. a neuropsychological scale) for each voxel in the mask image.
    :param images: a list of paths for the images
    :param scale: a list of scale values
    :param mask: the path to the mask image
    :param output_rs: the path to the output file for the correlation coefficients
    :param output_ps: the path to the output file for the p-values
    :param corr: the type of correlation coefficient to calculate (pearson, spearman)
    """

    mask_data =  mask.get_fdata()
    
    sample_img = images[0]
    sample_data = sample_img.get_fdata()

    corr_r = sample_data * 0
    corr_p = sample_data * 0

    data = sample_data[..., np.newaxis]

    for image in images[1:]:
        
        img_data = image.get_fdata()[..., np.newaxis]
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

    corr_r_img = nib.Nifti1Image(corr_r, sample_img.affine, header = sample_img.header)
    corr_p_img = nib.Nifti1Image(corr_p, sample_img.affine, header = sample_img.header)
    
    return corr_r_img, corr_p_img
    
def image_to_image_corr_atlas_based_spearman(image_1, image_2, atlas):
    """Calculates the spearman correlation coefficient and p-value between two images using
    the atlas ROI-values extracted for each region of interest (ROI) defined by an atlas.

    :param image_1: the path to the first image
    :param image_2: the path to the second image
    :param atlas: the path to the atlas image that defines the ROIs
    :return rho: correlation coefficient
    :return p: p-value
    """

    # Load the two images + atlas using nibabel library
    img_1_data = image_1.get_fdata()
    img_2_data = image_2.get_fdata()
    atlas_data = atlas.get_fdata()

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

def normalized_cross_correlation_2images(img1, img2, no_zeros = True):
    """
    Calculates the normalized cross correlation (NCC) between two images.
    :param img1_path: the path to the first image
    :param img2_path: the path to the second image
    :return: the NCC value between the two images
    """

    # Load Nifti images
    img1_data = img1.get_fdata()
    img2_data = img2.get_fdata()

    if no_zeros:
    # Get non-zero indices
        non_zero_mask = (img1_data != 0) & (img2_data != 0)
        img1_data = img1_data[non_zero_mask]
        img2_data = img2_data[non_zero_mask]

    # Subtract the mean from each image
    img1_data = img1_data - np.mean(img1_data)
    img2_data = img2_data - np.mean(img2_data)

    # Calculate the numerator of the NCC equation
    numerator = np.sum(img1_data * img2_data)

    # Calculate the denominator of the NCC equation  
    denominator = np.sqrt(np.sum(img1_data ** 2) * np.sum(img2_data ** 2))

    # Calculate and return the NCC
    return numerator / denominator
