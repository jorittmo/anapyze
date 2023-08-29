import nibabel as nib


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


def convert_format(img_name, compress_out=True):
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