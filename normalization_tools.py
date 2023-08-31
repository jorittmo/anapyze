import os
import numpy as np
import nibabel as nib
from scipy.fftpack import fftn, fftshift
from scipy import signal
import matplotlib.pyplot as plt
import pwlf
import spm


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

    fdg = nib.load(input_image)
    template = nib.load(template)
    mask = nib.load(mask)

    fdg_data = fdg.get_fdata()
    template_data = template.get_fdata()
    if len(mask.shape) == 4:
        mask_data = mask.get_fdata()[:, :, :, 0]
    else:
        mask_data = mask.get_fdata()[:, :, :]

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

    img = nib.Nifti1Image(fdg_data, fdg.affine, fdg.header)
    nib.save(img, output)

    return norm_value


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

    pons_img = nib.load(ref_region).get_fdata()
    if len(pons_img.shape) == 4:
        pons_img = pons_img[:, :, :, 0]
    pons_vox = np.where(pons_img == 1)

    input_img = nib.load(input_image)
    img_data = input_img.get_fdata()
    if len(input_img.shape) == 4:
        img_data = img_data[:, :, :, 0]

    pons_value = np.mean(img_data[pons_vox])
    normalized_img = img_data / pons_value

    new_img = nib.Nifti1Image(normalized_img, input_img.affine, input_img.header)
    nib.save(new_img, output_image)


def exact_histogram_matching(reference_nii, input_nii, output_nii, number_kernels=3, nbins=1024):
    """This function performs an exact histogram matching between a reference image and an input image,
    and saves the output image in a nifti file.
    More information: 10.1109/TIP.2005.864170

    :param reference_nii: the name of the nifti file that contains the reference image
    :type reference_nii: str
    :param input_nii: the name of the nifti file that contains the input image
    :type input_nii: str
    :param output_nii: the name of the nifti file that will contain the output image
    :type output_nii: str
    :param number_kernels: the number of kernels to use for the exact histogram matching, defaults to 3
    :type number_kernels: int
    :param nbins: the number of bins to use for the histogram matching, defaults to 1024
    :type nbins: int
    :return: the name of the nifti file that contains the output image
    :rtype: str
    """

    template = nib.load(reference_nii)
    nt_data = template.get_data()[:, :, :]
    scaled_nt_data = np.round((nbins - 1) * (nt_data / np.max(nt_data)))
    scaled_nt_data = np.asarray(scaled_nt_data, np.uint16)

    patient = nib.load(input_nii)
    pt_data = patient.get_data()[:, :, :]
    scaled_pt_data = np.round((nbins - 1) * (pt_data / np.max(pt_data)))
    scaled_pt_data = np.asarray(scaled_pt_data, np.uint16)

    reference_histogram = ExactHistogramMatcher.get_histogram(scaled_nt_data, 10)
    matched_img = ExactHistogramMatcher.match_image_to_histogram(pt_data, reference_histogram, number_kernels)

    img = nib.Nifti1Image(matched_img, patient.affine, patient.header)
    nib.save(img, output_nii)

    return output_nii


def bi_histogram_matching(reference_nii, input_nii, output_nii, nbins=256):
    """This function performs a bi-directional histogram matching between a reference image and an input image,
    and saves the output image in a nifti file.

    :param reference_nii: the name of the nifti file that contains the reference image
    :type reference_nii: str
    :param input_nii: the name of the nifti file that contains the input image
    :type input_nii: str
    :param output_nii: the name of the nifti file that will contain the output image
    :type output_nii: str
    :param nbins: the number of bins to use for the histogram matching, defaults to 256
    :type nbins: int
    :return: the path to the output normalized Nifti image
    :rtype: str
    """

    # Load and prepare the template image
    template = nib.load(reference_nii)
    nt_data = template.get_data()[:, :, :]
    indx = np.where(nt_data < 0)
    # Set array values in the interval to val
    nt_data[indx] = 0
    nt_max = np.amax(nt_data)
    nt_data = ((nbins - 1) / nt_max) * nt_data

    # Load and prepare the patient image
    patient = nib.load(input_nii)
    pt_data = patient.get_data()[:, :, :]
    indx = np.where(pt_data < 0)
    # Set array values in the interval to val
    pt_data[indx] = 0
    pt_max = np.amax(pt_data)
    pt_data = (255 / pt_max) * pt_data
    oldshape = pt_data.shape

    # We adjust first the low part of the histogram
    nt_data_low = nt_data
    indx = np.where(nt_data > (nbins / 2))
    nt_data_low[indx] = 0
    pt_data_low = pt_data
    indx = np.where(pt_data > (nbins / 2))
    pt_data_low[indx] = 0
    # We save data_low just for checking for the moment
    img = nib.Nifti1Image(pt_data_low, patient.affine, patient.header)
    nib.save(img, input_nii[0:-4] + '_low.nii')

    # Converts the data arrays to single dimension and histograms the first part of the data
    nt_data_array = nt_data_low.ravel()
    pt_data_array = pt_data_low.ravel()
    s_values, bin_idx, s_counts = np.unique(pt_data_array, return_inverse=True, return_counts=True)
    t_values, t_counts = np.unique(nt_data_array, return_counts=True)
    s_quantiles = np.cumsum(s_counts).astype(np.float64)
    s_quantiles /= s_quantiles[-1]
    t_quantiles = np.cumsum(t_counts).astype(np.float64)
    t_quantiles /= t_quantiles[-1]
    interp_t_values = np.interp(s_quantiles, t_quantiles, t_values)
    final_image_data_low = interp_t_values[bin_idx].reshape(oldshape)

    # We prepare now to adjust the second part of the histogram
    nt_data = template.get_data()[:, :, :]
    indx = np.where(nt_data < 0)
    nt_data[indx] = 0
    nt_max = np.amax(nt_data)
    nt_data = ((nbins - 1) / nt_max) * nt_data

    pt_data = patient.get_data()[:, :, :]
    indx = np.where(pt_data < 0)
    pt_data[indx] = 0
    pt_max = np.amax(pt_data)
    pt_data = ((nbins - 1) / pt_max) * pt_data
    oldshape = pt_data.shape

    nt_data_high = nt_data
    indx = np.where(nt_data <= (nbins / 2))
    nt_data_high[indx] = 0
    indx = np.where(nt_data > (nbins / 2))
    nt_data_high[indx] = nt_data_high[indx] - (nbins / 2)
    pt_data_high = pt_data
    indx = np.where(pt_data <= (nbins / 2))
    pt_data_high[indx] = 0
    indx_pth = np.where(pt_data_high > (nbins / 2))
    pt_data_high[indx_pth] = pt_data_high[indx_pth] - (nbins / 2)

    # Converts the data arrays to single dimension and histograms the second part of the data
    nt_data_array = nt_data_high.ravel()
    pt_data_array = pt_data_high.ravel()
    s_values, bin_idx, s_counts = np.unique(pt_data_array, return_inverse=True, return_counts=True)
    t_values, t_counts = np.unique(nt_data_array, return_counts=True)
    s_quantiles = np.cumsum(s_counts).astype(np.float64)
    s_quantiles /= s_quantiles[-1]
    t_quantiles = np.cumsum(t_counts).astype(np.float64)
    t_quantiles /= t_quantiles[-1]
    interp_t_values = np.interp(s_quantiles, t_quantiles, t_values)
    final_image_data_high = interp_t_values[bin_idx].reshape(oldshape)

    final_image_data_high[indx_pth] = final_image_data_high[indx_pth] + 128
    final_image_data = final_image_data_high + final_image_data_low

    # Saves the output data
    img = nib.Nifti1Image(final_image_data, patient.affine, patient.header)
    nib.save(img, output_nii)

    return output_nii


def logpow_histogram_matching(reference_nii, input_nii, output_nii, alpha=1, beta=3):
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

    # Load the template image
    template = nib.load(reference_nii)
    nt_data = template.get_data()[:, :, :]

    # Load the patient image
    patient = nib.load(input_nii)
    pt_data = patient.get_data()[:, :, :]

    # Stores the image data shape that will be used later
    oldshape = pt_data.shape

    # Converts the data arrays to single dimension and normalizes by the maximum
    nt_data_array = nt_data.ravel()
    pt_data_array = pt_data.ravel()

    # get the set of unique pixel values and their corresponding indices and counts
    s_values, bin_idx, s_counts = np.unique(pt_data_array, return_inverse=True, return_counts=True)
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

    # Saves the output data
    img = nib.Nifti1Image(final_image_data, patient.affine, patient.header)
    nib.save(img, output_nii)

    return output_nii


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

    # Load the template image
    template = nib.load(reference_nii)
    nt_data = template.get_data()[:, :, :]

    # Load the patient image
    patient = nib.load(input_nii)
    pt_data = patient.get_data()[:, :, :]

    # Stores the image data shape that will be used later
    oldshape = pt_data.shape

    # Converts the data arrays to single dimension and normalizes by the maximum
    nt_data_array = nt_data.ravel()
    pt_data_array = pt_data.ravel()

    # get the set of unique pixel values and their corresponding indices and counts
    s_values, bin_idx, s_counts = np.unique(pt_data_array, return_inverse=True, return_counts=True)
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

    # Saves the output data
    img = nib.Nifti1Image(final_image_data, patient.affine, patient.header)
    nib.save(img, output_nii)

    return


def optimize_smoothing_group_ncc(spm_path, img_group_paths, template_path, sigma_range):

    """Optimizes the smoothing of a group of images using a Gaussian filter and the spm library.

    :param spm_path: the path to the spm executable
    :type spm_path: str
    :param img_group_paths: a list of paths to the images to be smoothed
    :type img_group_paths: list
    :param template_path: the path to the image that is used as a template which the imgs at img_group_paths should match
    :type template_path: str
    :param sigma_range: a list of sigma values for the Gaussian filter
    :type sigma_range: list
    :return: the best sigma value for the Gaussian filter and the best normalized cross correlation value
    :rtype: float, float
    """

    best_ncc = -np.inf
    best_sigma = None

    template = nib.load(template_path)
    template_data = np.nan_to_num(template.get_fdata())
    indexes = np.where(template_data > 0.1 * np.max(template_data))

    template_data = template_data - np.mean(template_data[indexes])

    for sigma in sigma_range:
        print("Testing Gaussian Filter with FHWM %s" % sigma)
        ncc_sum = 0
        ncc_count = 0

        spm_proc = spm.spm(spm_path)
        spm_proc.smooth_imgs(img_group_paths, [sigma, sigma, sigma])

        smooth_img_group = []

        for img_path in img_group_paths:
            components = os.path.split(img_path)
            simg_ = os.path.join(components[0], 's' + components[1])
            smooth_img_group.append(simg_)

        for img_path in smooth_img_group:
            # Load Nifti images
            img1 = nib.load(img_path)
            # Get image data as numpy arrays
            img1_data = np.nan_to_num(img1.get_fdata())
            # Subtract the mean from each image
            img1_data = img1_data - np.mean(img1_data[indexes])

            # Calculate the numerator of the NCC equation
            numerator = np.sum(img1_data[indexes] * template_data[indexes])
            # Calculate the denominator of the NCC equation
            denominator = np.sqrt(np.sum(img1_data[indexes] ** 2) * np.sum(template_data[indexes] ** 2))
            # Calculate the NCC
            ncc_value = numerator / denominator
            ncc_sum += ncc_value
            ncc_count += 1

        avg_ncc = ncc_sum / ncc_count
        # Check if this NCC value is the highest so far
        if avg_ncc > best_ncc:
            best_ncc = avg_ncc
            best_sigma = sigma
    return best_sigma, best_ncc


def estimate_fwhm_mizutani(nifti, bin_size=5, n_segs=2, orientation='axial', plot=False):

    """
    This function estimates image xy resolution based in a previous method implemented by MIZUTANI ET AL.
    for microscopy imaging:

    https://onlinelibrary.wiley.com/doi/10.1111/jmi.12315

    The idea of expanding this method to PET imaging was proposed by Carbonell et al:
    https://n.neurology.org/content/94/15_Supplement/4301

    The FWHMs are estimated from multiple regression of the logarithm of the square norm of
    the image Fourier transform against the square distance from the origin in the Fourier domain
    (essentially from what is known as a Wilson plot). The lower frequencies in the Wilson plot are
    linearly fitted and the PSF is estimated from the slope.

    :param orientation: plane to walk traverse the image (axial, sagittal, coronal).
    :param nifti: input nifti image

    :param bin_size: In the original paper from Mizutani, the logarithm of the average squared norm in each
    5 × 5 pixel bin of the Fourier transform was used. You can vary this to check what fits your data better.
    Values between 5 and 10 seem to work well for PET images.

    :param plot: If True, the function will plot the Wilson plots and fitted lines for quality check.

    :return: The average FWHM across all slices in the image.
    """



    img = nib.load(nifti)
    image_xy_size = np.abs(img.affine[0, 0])
    volume = img.get_fdata()
    fwhm_values = []
    sigma_values = []

    if orientation == 'axial':
        orient = 2
    elif orientation == 'coronal':
        orient = 1
    elif orientation == 'sagittal':
        orient = 0
    else:
        raise ValueError('orientation must be set to axial, coronal or sagital')

    for i in range(volume.shape[orient]):

        if orientation == 'axial':
            slice_2d = volume[:, :, i]
        elif orientation == 'coronal':
            slice_2d = volume[:, i, :]
        elif orientation == 'sagittal':
            slice_2d = volume[i, :, :]

        # Compute the 2D Fourier transform of the input slice
        fourier_slice = fftshift(fftn(slice_2d))

        # Compute the square norm of the Fourier transform
        square_norm = np.abs(fourier_slice) ** 2

        # Calculate logarithm of the square norm
        log_square_norm = np.log(square_norm + np.finfo(float).eps)

        # Calculate k^2 values (distances from the origin in Fourier domain)
        k_values = np.array(np.meshgrid(
            np.fft.fftshift(np.fft.fftfreq(slice_2d.shape[0])),
            np.fft.fftshift(np.fft.fftfreq(slice_2d.shape[1])),
            indexing='ij'
        ))
        k_squared = np.sum(k_values ** 2, axis=0)

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
                    log_square_norm[x * bin_size:(x + 1) * bin_size, y * bin_size:(y + 1) * bin_size])
                k_squared_binned[x, y] = np.mean(
                    k_squared[x * bin_size:(x + 1) * bin_size, y * bin_size:(y + 1) * bin_size])

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
        sigma = np.sqrt(-slope / (4 * np.pi ** 2))
        sigma_values.append(sigma)

        # Calculate the FWHM of the PSF
        fwhm = 2 * np.sqrt(2 * np.log(2)) * sigma * image_xy_size
        fwhm_values.append(fwhm)

        if i == volume.shape[2] // 3 and plot:
            # Plot the Wilson plot with the fitted line overlaid
            plt.scatter(k_squared_flat, log_square_norm_flat, s=1)
            plt.plot(x_pred, y_pred, color='black', linewidth=1)
            plt.xlabel('|k|^2')
            plt.ylabel('ln|F(k)|^2')
            plt.title('Wilson plot with fitted line')
            plt.ylim([np.min(log_square_norm_flat) - 5, np.max(log_square_norm_flat) + 5])
            plt.show()

    return np.nanmean(fwhm_values), np.nanmean(sigma_values)


class ExactHistogramMatcher:
    _kernel1 = 1.0 / 5.0 * np.array([[0, 1, 0],
                                     [1, 1, 1],
                                     [0, 1, 0]])

    _kernel2 = 1.0 / 9.0 * np.array([[1, 1, 1],
                                     [1, 1, 1],
                                     [1, 1, 1]])

    _kernel3 = 1.0 / 13.0 * np.array([[0, 0, 1, 0, 0],
                                      [0, 1, 1, 1, 0],
                                      [1, 1, 1, 1, 1],
                                      [0, 1, 1, 1, 0],
                                      [0, 0, 1, 0, 0]])

    _kernel4 = 1.0 / 21.0 * np.array([[0, 1, 1, 1, 0],
                                      [1, 1, 1, 1, 1],
                                      [1, 1, 1, 1, 1],
                                      [1, 1, 1, 1, 1],
                                      [0, 1, 1, 1, 0]])

    _kernel5 = 1.0 / 25.0 * np.array([[1, 1, 1, 1, 1],
                                      [1, 1, 1, 1, 1],
                                      [1, 1, 1, 1, 1],
                                      [1, 1, 1, 1, 1],
                                      [1, 1, 1, 1, 1]])
    _kernel_mapping = {1: [_kernel1],
                       2: [_kernel1, _kernel2],
                       3: [_kernel1, _kernel2, _kernel3],
                       4: [_kernel1, _kernel2, _kernel3, _kernel4],
                       5: [_kernel1, _kernel2, _kernel3, _kernel4, _kernel5]}

    @staticmethod
    def get_histogram(image, image_bit_depth=8):
        """
        :param image: image as numpy array
        :param image_bit_depth: bit depth of the image. Most images have 8 bit.
        :return:
        """
        max_grey_value = pow(2, image_bit_depth)

        if len(image.shape) == 3:
            dimensions = image.shape[2]
            hist = np.empty((max_grey_value, dimensions))

            for dimension in range(0, dimensions):
                for gray_value in range(0, max_grey_value):
                    image_2d = image[:, :, dimension]
                    hist[gray_value, dimension] = len(image_2d[image_2d == gray_value])
        else:
            hist = np.empty((max_grey_value,))

            for gray_value in range(0, max_grey_value):
                hist[gray_value] = len(image[image == gray_value])

        return hist

    @staticmethod
    def _get_averaged_images(img, kernels):
        return np.array([signal.convolve2d(img, kernel, 'same') for kernel in kernels])

    @staticmethod
    def _get_average_values_for_every_pixel(img, number_kernels):
        """
        :param img: the image to be used in order to calculate averaged images
        :param number_kernels: number of kernels to be used in order to calculate the averaged images
        :return: averaged images with the shape:
                 (image height * image width, number averaged images)
                 Every row represents one pixel and its averaged values.
                 I. e. x[0] represents the first pixel and contains an array with k
                 averaged pixels where k are the number of used kernels.
        """
        kernels = ExactHistogramMatcher._kernel_mapping[number_kernels]
        averaged_images = ExactHistogramMatcher._get_averaged_images(img, kernels)
        img_size = averaged_images[0].shape[0] * averaged_images[0].shape[1]

        # shape of averaged_images: (number averaged images, height, width).
        # Reshape in a way, that one row contains all averaged values of pixel in position (x, y)
        reshaped_averaged_images = averaged_images.reshape((number_kernels, img_size))
        transposed_averaged_images = reshaped_averaged_images.transpose()
        return transposed_averaged_images

    @staticmethod
    def sort_rows_lexicographically(matrix):
        # Because lexsort in numpy sorts after the last row,
        # then after the second last row etc., we have to rotate
        # the matrix in order to sort all rows after the first column,
        # and then after the second column etc.

        rotated_matrix = np.rot90(matrix)

        # TODO lexsort is very memory hungry! If the image is too big, this can result in SIG 9!
        sorted_indices = np.lexsort(rotated_matrix)
        return matrix[sorted_indices]

    @staticmethod
    def _match_to_histogram(image, reference_histogram, number_kernels):
        """
        :param image: image as numpy array.
        :param reference_histogram: reference histogram as numpy array
        :param number_kernels: The more kernels you use in order to calculate average images,
                               the more likely it is, the resulting image will have the exact
                               histogram like the reference histogram
        :return: The image with the exact reference histogram.
        """
        img_size = image.shape[0] * image.shape[1]

        merged_images = np.empty((img_size, number_kernels + 2))

        # The first column are the original pixel values.
        merged_images[:, 0] = image.reshape((img_size,))

        # The last column of this array represents the flattened image indices.
        # These indices are necessary to keep track of the pixel positions
        # after they haven been sorted lexicographically according their values.
        indices_of_flattened_image = np.arange(img_size).transpose()
        merged_images[:, -1] = indices_of_flattened_image

        # Calculate average images and add them to merged_images
        averaged_images = ExactHistogramMatcher._get_average_values_for_every_pixel(image, number_kernels)
        for dimension in range(0, number_kernels):
            merged_images[:, dimension + 1] = averaged_images[:, dimension]

        # Sort the array according the original pixels values and then after
        # the average values of the respective pixel
        sorted_merged_images = ExactHistogramMatcher.sort_rows_lexicographically(merged_images)

        # Assign gray values according the distribution of the reference histogram
        index_start = 0
        for gray_value in range(0, len(reference_histogram)):
            index_end = int(index_start + reference_histogram[gray_value])
            sorted_merged_images[index_start:index_end, 0] = gray_value
            index_start = index_end

        # Sort back ordered by the flattened image index. The last column represents the index
        sorted_merged_images = sorted_merged_images[sorted_merged_images[:, -1].argsort()]
        new_target_img = sorted_merged_images[:, 0].reshape(image.shape)

        return new_target_img

    @staticmethod
    def match_image_to_histogram(image, reference_histogram, number_kernels=3):
        """
        :param image: image as numpy array.
        :param reference_histogram: reference histogram as numpy array
        :param number_kernels: The more kernels you use in order to calculate average images,
                               the more likely it is, the resulting image will have the exact
                               histogram like the reference histogram
        :return: The image with the exact reference histogram.
                 CAUTION: Don't save the image in a lossy format like JPEG,
                 because the compression algorithm will alter the histogram!
                 Use lossless formats like PNG.
        """
        if len(image.shape) == 3:
            # Image with more than one dimension. I. e. an RGB image.
            output = np.empty(image.shape)
            dimensions = image.shape[2]

            for dimension in range(0, dimensions):
                output[:, :, dimension] = ExactHistogramMatcher._match_to_histogram(image[:, :, dimension],
                                                                                    reference_histogram[:, dimension],
                                                                                    number_kernels)
        else:
            # Gray value image
            output = ExactHistogramMatcher._match_to_histogram(image,
                                                               reference_histogram,
                                                               number_kernels)

        return output
