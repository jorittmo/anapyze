addpath('SPM_PATH')

matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.source = {'SOURCE_IMAGE,1'};
matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.wtsrc = '';
matlabbatch{1}.spm.tools.oldnorm.estwrite.subj.resample = {RESAMPLE_IMAGES};
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.template = {'TEMPLATE,1'};
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.weight = '';
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smosrc = 8;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.smoref = 0;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.regtype = 'mni';
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.cutoff = 25;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.nits = 16;
matlabbatch{1}.spm.tools.oldnorm.estwrite.eoptions.reg = 1;
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.preserve = 0;
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.bb = [BOUNDING_BOX];
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.vox = [VOXELSIZE VOXELSIZE VOXELSIZE];
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.interp = INTERPOLATION;
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.wrap = [WRAP_XYZ];
matlabbatch{1}.spm.tools.oldnorm.estwrite.roptions.prefix = 'PREFIX';

spm('defaults','fmri');
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);