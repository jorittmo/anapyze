addpath('SPM_PATH')

matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {'REFERENCE_IMAGE,1'};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {'SOURCE_IMAGE,1'};
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {OTHER_IMAGES};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [WRAP_XYZ];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'PREFIX';

spm('defaults','fmri');
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);