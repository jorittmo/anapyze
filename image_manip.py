import nibabel as nib

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

def convert_dicom_to_nifti(input_folder, output_folder, output_filename):
    # Create the output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    # Run dcm2niix command to convert DICOM to NIfTI
    command = f"dcm2niix -o {output_folder} -f {output_filename} {input_folder}"
    subprocess.run(command, shell=True)

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


