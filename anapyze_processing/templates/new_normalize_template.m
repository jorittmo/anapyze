addpath('SPM_PATH')

matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = {'IMAGE_TO_NORM,1'};
matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = {RESAMPLE_IMAGES};
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {'SPM_PATH/tpm/TPM.nii'};
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0 0.1 0.01 0.04];
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [BOUNDING_BOX];
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = [VOXELSIZE VOXELSIZE VOXELSIZE];
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = INTERPOLATION;
matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.prefix = 'PREFIX';

spm('defaults','fmri');
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);