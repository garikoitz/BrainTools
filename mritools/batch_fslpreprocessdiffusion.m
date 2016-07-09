function batch_fslpreprocessdiffusion(subName, DWIdir, retdir)
% GLU MINI project adapted from: BDE lab preprocessing for diffusion data
%
%
% Inputs
% basedir     - path to the base directory with the raw data
% t1dir       - path to the ac-pc aligned t1. Everything will be
%               coregistered to this
% doMakeNifti - Logical. Whether or not to make niftis from par/rec files
% doPreProc   - logical. Whether or not to do preprocessing
%
% Example:
%
% DWIdir = '/bcbl/home/public/Gari/MINI/ANALYSIS/DWI'
% retdir = '/bcbl/home/public/Gari/MINI/ANALYSIS/ret'
% 
% Edited by GLU on June 2016
% update instructions
% Now it only is for MINI fsl preprocessing


%% Set up directories and find files
basedir = fullfile(DWIdir,subName);
rawdir = fullfile(DWIdir,subName, 'raw');
if ~exist(rawdir, 'dir'), mkdir(rawdir), end
t1dir = fullfile(retdir, subName, 'anat');

% Note glu: conversion from dicom in the ipython notebook file with qsubs
d30 = dir(fullfile(rawdir,'*d35b1000*.nii'));
d60 = dir(fullfile(rawdir,'*d65b2500*.nii'));
b0 = dir(fullfile(rawdir,'*d6b0*.nii')); % grab post-anterior encoded file 


% temp: note that this pulls only 1 of each file type, some subjects have
% e.g. repeated measures for 64 or 32 dir data in a session
% Note glu: just left one type per every subject to avoid problems

% After making it work with 30 and 60 separately, now I am going to
% concatenate everything so that to mrTrix it arrives only one file. As the
% b is different the model has to consider it, only available in mrTrix3
% Go from smaller to bigger, consider it for the pe matrix
% dMRI60Files{1}=fullfile(rawdir,d60(1).name);
% dMRI30Files{1}=fullfile(rawdir,d30(1).name);
% Add the b0 with the reversed phase encode
% for ii = 1:length(b0)
%     dMRI60Files{1+ii}=fullfile(rawdir,b0(ii).name);
%     dMRI30Files{1+ii}=fullfile(rawdir,b0(ii).name);
% end 
% % Bvals and Bvecs files
% for ii = 1:length(dMRI64Files)
%     bvals60{ii} = [prefix(prefix(dMRI60Files{ii})) '_bvals'];
%     bvecs60{ii} = [prefix(prefix(dMRI60Files{ii})) '_bvecs'];
%     bvals30{ii} = [prefix(prefix(dMRI30Files{ii})) '_bvals'];
%     bvecs30{ii} = [prefix(prefix(dMRI30Files{ii})) '_bvecs'];
% end

dMRIFiles{1}=fullfile(rawdir,b0(1).name);
dMRIFiles{2}=fullfile(rawdir,d30(1).name);
dMRIFiles{3}=fullfile(rawdir,d60(1).name);
% Bvals and Bvecs files
for ii = 1:length(dMRIFiles)
    bvals{ii} = [prefix(prefix(dMRIFiles{ii})) '_bvals'];
    bvecs{ii} = [prefix(prefix(dMRIFiles{ii})) '_bvecs'];
end



% Phase encode matrix. This denotes, for each volume, which direction is
% the phase encode
% Edit GLU: See here the info used for MINI data in Siemens TRIO
% http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/TOPUP/Faq
% We have 5xb0 + 30 dirs & 5xb0 + 60 dirs with one A>>P (this is y=-1)
% and 6xb0 with P>>A (this is y=1)
% The dwell time is exactly the same of the example = 0.095 
% (Calculated from echo spacing 0.75 and EPI factor 128)
pe_mat = [0 1 0; 0 -1 0; 0 -1 0];
dwellTime = 0.095;

% Directory to save everything
% outdir60 = fullfile(basedir,'dmri60');
% outdir30 = fullfile(basedir,'dmri30');
outdir = fullfile(basedir,'dmri');

% break

%% Pre process: This is mostly done with command line calls to FSL
fsl_preprocess(dMRIFiles, bvecs, bvals, pe_mat, outdir, dwellTime);




