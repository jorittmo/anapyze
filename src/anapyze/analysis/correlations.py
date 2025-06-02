import numpy as np
import nibabel as nib
from scipy.stats import pearsonr, spearmanr

def voxel_wise_corr_images_vs_scale(images, scale, mask, corr = "pearson"):
    """
    Computes voxel-wise correlation between a list of 3D NIfTI images and a set of scalar values.

    Parameters
    ----------
    images : List[nibabel.Nifti1Image]
        List of loaded NIfTI image objects (all must have the same shape and affine). Each element
        represents a subject’s image for which you want to correlate voxel intensities with the
        provided scale values.
    scale : Sequence[float]
        A sequence of numeric values (e.g., neuropsychological test scores) of the same length as
        `images`. Each entry corresponds to the scalar value for the subject of the same index
        in `images`.
    mask : nibabel.Nifti1Image
        A binary NIfTI image whose nonzero voxels define the mask region. Correlation is computed
        only for voxels where `mask.get_fdata() != 0`. Must have the same spatial dimensions as
        the images in `images`.
    corr : {"pearson", "spearman"}, default "pearson"
        Specifies which correlation coefficient to compute at each voxel:
        - `"pearson"`: Pearson’s r and two-tailed p-value via `scipy.stats.pearsonr`.
        - `"spearman"`: Spearman’s rho and two-tailed p-value via `scipy.stats.spearmanr`.

    Returns
    -------
    corr_r_img : nibabel.Nifti1Image
        A NIfTI image where each voxel contains the correlation coefficient (r or rho) between
        the intensity values across subjects and the corresponding `scale` values. Voxels outside
        the mask are set to zero.
    corr_p_img : nibabel.Nifti1Image
        A NIfTI image where each voxel contains the two-tailed p-value associated with the
        computed correlation. Voxels outside the mask are set to zero.
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
    """
    Calculates the Spearman correlation coefficient and p-value between two 3D NIfTI images
    by averaging voxel intensities within each ROI defined by an atlas.

    Parameters
    ----------
    image_1 : nibabel.Nifti1Image
        First input NIfTI image. Must share the same spatial dimensions and affine with
        `image_2` and `atlas`.
    image_2 : nibabel.Nifti1Image
        Second input NIfTI image. Must align with `image_1` and `atlas`.
    atlas : nibabel.Nifti1Image
        Atlas NIfTI image defining ROIs. Nonzero integer labels indicate distinct regions;
        voxels labeled 0 are treated as background and ignored.

    Returns
    -------
    rho : float
        Spearman correlation coefficient between the mean ROI values of `image_1` and `image_2`.
    p : float
        Two‐tailed p‐value associated with the computed Spearman correlation.
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
    Computes the Normalized Cross-Correlation (NCC) between two 3D NIfTI images.

    Parameters
    ----------
    img1 : nibabel.Nifti1Image
        First input NIfTI image. Must have the same dimensions as `img2`.
    img2 : nibabel.Nifti1Image
        Second input NIfTI image. Must have the same dimensions as `img1`.
    no_zeros : bool, default True
        If True, only voxels where both images are nonzero will be included in the calculation.
        If False, all voxels (including zeros) are used.

    Returns
    -------
    float
        The normalized cross-correlation coefficient between `img1` and `img2`. Value range is
        between -1 and 1.
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
