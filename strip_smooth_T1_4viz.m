function strip_smooth_T1_4viz(RibbonIm, MaskIm, T1Im, OutputFolder)
% strip_smooth_T1_4viz: based on a cortical ribbon (from FreeSurfer), a
% gray matter tissue class and T1w images, write a stripped T1 (where
% the mask combines GM mask and ribbon) and
% smooth the outer part of the cortical surface such that it looks smooth
% in Cartool while preserving details in the inner part of the stripped
% brain
%
% Usage:
%-------------------------------------------------------------------------
% 
% strip_smooth_T1_4viz(RibbonIm, MaskIm, T1Im, OutputFolder)
%
% Inputs:
%-------------------------------------------------------------------------
%
% RibbonIm    : path to cortical ribbon image (.nii format)
% MaskIm      : path to gray matter mask image (.nii format)
% T1Im        : path to T1-weighted image (.nii format)
% OutputFolder: path to output folder for s2mm_stripped_T1.nii image
%
% Outputs:
%-------------------------------------------------------------------------
%
% s2mm_stripped_T1.nii file inside OutputFolder
%
%-------------------------------------------------------------------------
% Cartool: https://sites.google.com/site/cartoolcommunity
%-------------------------------------------------------------------------
% Renaud Marquis @ FBMlab, April 2018, updated September 2019
%-------------------------------------------------------------------------

% RibbonIm = fullfile(OutputFolder,'ribbon.nii');
% MaskIm = fullfile(OutputFolder,'T1w_class-GM_ribbon_hippo_amyg.nii');
% T1Im = fullfile(OutputFolder,'T1.nii');
Vtemp1name = 'temp_stripped_T1.nii';
prefixVtemp2 = 's2mm_';
Vtemp1 = fullfile(OutputFolder,Vtemp1name);
Vtemp2 = spm_file(Vtemp1,'prefix',prefixVtemp2);

matlabbatch{1}.spm.util.imcalc.input = {
                                        [RibbonIm,',1']
                                        [MaskIm,',1']
                                        [T1Im,',1']
                                        };
matlabbatch{1}.spm.util.imcalc.output = Vtemp1name;
matlabbatch{1}.spm.util.imcalc.outdir = {OutputFolder};
matlabbatch{1}.spm.util.imcalc.expression = '(((i1>0)+i2)>0).*i3';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
matlabbatch{2}.spm.spatial.smooth.data(1) = cfg_dep('Image Calculator: ImCalc Computed Image: output', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{2}.spm.spatial.smooth.fwhm = [2 2 2];
matlabbatch{2}.spm.spatial.smooth.dtype = 0;
matlabbatch{2}.spm.spatial.smooth.im = 0;
matlabbatch{2}.spm.spatial.smooth.prefix = prefixVtemp2;

spm_jobman('run',matlabbatch);

V = spm_vol(Vtemp1);
D = spm_read_vols(V);

SE = strel('sphere',1);
E = imerode((D>0),SE);
% V2 = V; V2.fname = spm_file(V2.fname,'prefix','eroded_');
% spm_write_vol(V2,E);

ED = E.*D;
V3 = spm_vol(Vtemp2);
S = spm_read_vols(V3);
NZ = ED > 0;
S2=S;S2(NZ)=ED(NZ);

V4 = V;
V4.fname = spm_file(V4.fname,'basename',[prefixVtemp2,'stripped_T1']);
spm_write_vol(V4,S2);

delete(Vtemp1,Vtemp2);%,V2.fname);

