import os
import subprocess
import nibabel as nib
import numpy as np

import pydicom
from pydicom2nifti.convert_dicom import dicom_series_to_nifti
from scipy.ndimage import zoom

# This is an effort to convert all the fruitcake library to python

def convert_format_nii_hdr(img_name, compress_out=True):
    """This function converts an image file from nifti to analyze and viceversa.
    Depending on your input, it will output the other format. It also works with .nii.gz

    :param img_name: the name of the nii/nii.gz/hdr file to be converted
    :type img_name: str
    :param compress_out: Will out nii.gz if the input is analyze
    :type img_name: bool
    :return: the name of the result file that was created
    :rtype: str
    """

    # Load the NIFTI image using nibabel
    img = nib.load(img_name)
    data = img.get_fdata()

    if img[-4:] == ".hdr" or img[-4:] == '.img':

        out_img = nib.Nifti1Image(data, img.affine, img.header)

        if compress_out:
            out_name = img[-4:] + '.nii.gz'
        else:
            out_name = img[-4:] + '.nii'

        nib.save(out_img,out_name)
        print(f"File {img_name} was converted to Nifti format and saved as {out_name}")

        return out_name

    elif img[-4:] == ".nii" or img[-7:] == '.nii.gz':
    
        analyze_img = nib.analyze.AnalyzeImage(data, img.affine, img.header)

        if img[-4:] == ".nii":
            out_name = img[-4:] + '.hdr'
        else:
            out_name = img[-7:] + '.hdr'

        nib.save(analyze_img, out_name)
        print(f"File {img_name} was converted to Analyze format and saved as {out_name}")

        return out_name
    
    else:

        raise TypeError("Wrong format")

def convert_dicom_to_nifti_dcm2niix(input_folder, output_folder, output_filename):
    # Create the output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    # Run dcm2niix command to convert DICOM to NIfTI
    command = f"dcm2niix -o {output_folder} -f {output_filename} {input_folder}"
    subprocess.run(command, shell=True)

def convert_dicom_to_nifti_python(input_folder, output_folder, output_filename):
    # Convert DICOM series to NIfTI using pydicom2nifti
    dicom_series_to_nifti(input_folder, output_folder)
    
    # Rename the output NIfTI file
    nifti_files = [f for f in os.listdir(output_folder) if f.endswith('.nii.gz')]
    if len(nifti_files) == 1:
        old_path = os.path.join(output_folder, nifti_files[0])
        new_path = os.path.join(output_folder, output_filename + '.nii.gz')
        os.rename(old_path, new_path)

def add_poisson_noise(input_filepath, output_filepath, intensity_scaling_factor=1.0):
    # Load NIfTI or Analyze image
    img = nib.load(input_filepath)
    img_data = img.get_fdata()

    # Calculate the mean of the image data
    img_mean = np.mean(img_data)

    # Apply Poisson noise to the image
    noisy_img_data = np.random.poisson(img_data * intensity_scaling_factor)

    # Create a new NIfTI or Analyze image with the noisy data
    if isinstance(img, nib.nifti1.Nifti1Image):
        noisy_img = nib.Nifti1Image(noisy_img_data, img.affine, img.header)
    
    elif isinstance(img, nib.analyze.AnalyzeImage):
        
        noisy_img = nib.analyze.AnalyzeImage(noisy_img_data, img.affine, img.header)

    # Save the noisy image to the specified output filepath
    nib.save(noisy_img, output_filepath)

def apply_constant_to_img(image, c, operation, output=False):
    """This function applies a constant to an image using a specified operation, and saves the result in a new image.

    :param image: the name of the image file to be processed
    :type image: str
    :param c: the constant to be applied to the image
    :type c: float
    :param operation: the operation to be performed on the image. It can be 'mult', 'div' or 'sum'
    :type operation: str
    :param output: the name of the output image file. If False, the output file will have the same name as the input file
    :type output: str or bool
    :return: None
    """
    if output:
        output_file = output

    else:
        output_file = image

    img_ = nib.load(image)
    data = img_.get_fdata()

    if operation == 'mult':

        data = data*c

    elif operation == 'div':

        data = data/c

    elif operation == 'sum':

        data = data+c

    new_img = nib.AnalyzeImage(data, img_.affine, img_.header)
    nib.save(new_img, output_file)

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
    img = nib.load(input_filepath)
    img_data = img.get_fdata()

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

def resample_image_by_matrix_size(input_filepath, output_filepath, target_shape, interpolation='linear'):
    
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
    img = nib.load(input_filepath)
    img_data = img.get_fdata()

    # Calculate zoom factors for each dimension
    zoom_factors = np.array(target_shape) / np.array(img_data.shape)

    # Determine the order of interpolation based on the given method
    interpolation_methods = {
        'nearest': 0,
        'linear': 1,
        'cubic': 3,
        'quadratic': 4
        # Add more interpolation methods and their corresponding order values here
    }
    order = interpolation_methods.get(interpolation, 1)  # Default to linear if method not recognized

    # Resample the image data using scipy's zoom function
    resampled_img_data = zoom(img_data, zoom_factors, order=order)

    # Update the affine matrix to reflect the new voxel dimensions
    voxel_sizes = np.array(img.header.get_zooms()) * (np.array(img_data.shape) / target_shape)
    new_affine = img.affine.copy()
    np.fill_diagonal(new_affine, voxel_sizes)

    # Create a new NIfTI or Analyze image with the resampled data and updated affine
    if isinstance(img, nib.nifti1.Nifti1Image):
        resampled_img = nib.Nifti1Image(resampled_img_data, new_affine, header=img.header)
    elif isinstance(img, nib.analyze.AnalyzeImage):
        header = img.header.copy()
        header.set_data_shape(target_shape)
        resampled_img = nib.analyze.AnalyzeImage(resampled_img_data, new_affine, header)

    # Save the resampled image to the specified output filepath
    nib.save(resampled_img, output_filepath)

def resample_image_by_voxel_sizes(input_filepath, output_filepath, target_voxel_sizes, interpolation='linear'):
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
        'nearest': 0,
        'linear': 1,
        'cubic': 3,
        'quadratic': 4
        # Add more interpolation methods and their corresponding order values here
    }
    order = interpolation_methods.get(interpolation, 1)  # Default to linear if method not recognized

    # Resample the image data using scipy's zoom function
    resampled_img_data = zoom(img_data, zoom_factors, order=order)

    # Update the affine matrix to reflect the new voxel sizes
    new_affine = img.affine.copy()
    np.fill_diagonal(new_affine, target_voxel_sizes)

    # Create a new NIfTI or Analyze image with the resampled data and updated affine
    if isinstance(img, nib.nifti1.Nifti1Image):
        resampled_img = nib.Nifti1Image(resampled_img_data, new_affine, header=img.header)
    elif isinstance(img, nib.analyze.AnalyzeImage):
        header = img.header.copy()
        resampled_img = nib.analyze.AnalyzeImage(resampled_img_data, new_affine, header)

    # Save the resampled image to the specified output filepath
    nib.save(resampled_img, output_filepath)

def remove_nan_negs(input_filepath, output_filepath):
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








