function bde_preprocessdiffusion(basedir, t1dir, doMakeNifti, doPreProc)
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
% basedir = '/Users/gari/Documents/BCBL_PROJECTS/MINI/ANALYSIS/DWI/S011'
% t1dir = '/Users/gari/Documents/BCBL_PROJECTS/MINI/ANALYSIS/ret/S011/anat'
%  doMakeNifti = 0;
%  doPreproc = 1;
%
% bde_preprocessdiffusion(basedir, t1dir)
% TO DO: streamline multi subs, par to nifti function, sge, dim error
% 
% Edited by GLU on June 2016

%% Argument checking
if ~exist('doMakeNifti','var') || isempty(doMakeNifti)
    % convert parrec files in raw directory to nifti (or =0, skip if these files already exist)
    doMakeNifti = 0;
end
if ~exist('doPreProc','var') || isempty(doPreProc)
    doPreProc = 0; % skip pre-processing in FSL if these files already exist
end

%% Set up directories and find files
rawdir = fullfile(basedir,'raw');
if ~exist(rawdir, 'dir'), mkdir(rawdir), end


% Note glu: conversion from dicom in the ipython notebook file with qsubs
d30 = dir(fullfile(rawdir,'*d35b1000*.nii'));
d60 = dir(fullfile(rawdir,'*d65b2500*.nii'));
b0 = dir(fullfile(rawdir,'*d6b0*.nii')); % grab post-anterior encoded file 


% temp: note that this pulls only 1 of each file type, some subjects have
% e.g. repeated measures for 64 or 32 dir data in a session
% Note glu: just left one type per every subject to avoid problems
dMRI60Files{1}=fullfile(rawdir,d60(1).name);
dMRI30Files{1}=fullfile(rawdir,d30(1).name);

% Add the b0 with the reversed phase encode
for ii = 1:length(b0)
    dMRI60Files{1+ii}=fullfile(rawdir,b0(ii).name);
    dMRI30Files{1+ii}=fullfile(rawdir,b0(ii).name);
end

% Bvals and Bvecs files
for ii = 1:length(dMRI64Files)
    bvals60{ii} = [prefix(prefix(dMRI60Files{ii})) '_bvals'];
    bvecs60{ii} = [prefix(prefix(dMRI60Files{ii})) '_bvecs'];
    bvals30{ii} = [prefix(prefix(dMRI30Files{ii})) '_bvals'];
    bvecs30{ii} = [prefix(prefix(dMRI30Files{ii})) '_bvecs'];
end

% Phase encode matrix. This denotes, for each volume, which direction is
% the phase encode
% Edit GLU: See here the info used for MINI data in Siemens TRIO
% http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/TOPUP/Faq
% We have 5xb0 + 30 dirs & 5xb0 + 60 dirs with one A>>P
% and 6xb0 with P>>A
% The dwell time is exactly the same of the example = 0.095 
% (Calculated from echo spacing 0.75 and EPI factor 128)
pe_mat = [0 -1 0; 0 1 0];

% Directory to save everything
outdir60 = fullfile(basedir,'dmri60');
outdir30 = fullfile(basedir,'dmri30');

% break

%% Pre process: This is mostly done with command line calls to FSL
if doPreProc
    fsl_preprocess(dMRI60Files, bvecs60, bvals60, pe_mat, outdir60);
    fsl_preprocess(dMRI30Files, bvecs30, bvals30, pe_mat, outdir30);
end

%% Run dtiInit to fit tensor model
% Turn off motion and eddy current correction, that was taken care of by FSL

% Set up t1 path and the params that are common
t1 = fullfile(t1dir,'t1_acpc.nii.gz'); % Path to the t1-weighted image
params = dtiInitParams; % Set up parameters for controlling dtiInit
params.eddyCorrect=-1; % This turns off eddy current and motion correction
params.rotateBvecsWithCanXform=1; % Siemens data requires this to be 1
params.phaseEncodeDir=2; % AP phase encode, 1 is for R>>L (2 = A/P 'col')
params.clobber=1; % Overwrite anything previously done
params.fitMethod='rt'; % 'ls, or 'rt' for robust tensor fitting (longer)

% First for 60 dir data 
dtEddy = fullfile(outdir60,'eddy','data.nii.gz'); % Path to the data
params.bvalsFile = fullfile(outdir60,'eddy','bvals'); % Path to bvals
params.bvecsFile = fullfile(outdir60,'eddy','bvecs'); % Path to the bvecs
dt6FileName{1} = dtiInit(dtEddy,t1,params); % Run dtiInit to preprocess data

% Then for 30 dir data
dtEddy = fullfile(outdir30,'eddy','data.nii.gz'); % Path to the data
params.bvalsFile = fullfile(outdir30,'eddy','bvals'); % Path to bvals
params.bvecsFile = fullfile(outdir30,'eddy','bvecs'); % Path to the bvecs
dt6FileName{2} = dtiInit(dtEddy,t1,params); % Run dtiInit to preprocess data

%% Run AFQ

% Cell array with paths to the dt6 directories
% % dt6dirs = horzcat(fileparts(dt6FileName{1}), fileparts(dt6FileName{2}));


dt6dirs = horzcat({fileparts(dt6FileName{1}{1})}, {fileparts(dt6FileName{2}{1})});
% dt6dirs = horzcat({'/Users/gari/Documents/BCBL_PROJECTS/MINI/ANALYSIS/DWI/S011/dmri60/dti60trilin'},{'/Users/gari/Documents/BCBL_PROJECTS/MINI/ANALYSIS/DWI/S011/dmri30/dti30trilin'})



% afq = AFQ_Create('sub_dirs',dt6dirs,'sub_group',[0 0],'clip2rois', 0);
% To run AFQ in test mode so it will go quickly
% afq = AFQ_Create('sub_dirs',dt6dirs,'sub_group',[0 0],'run_mode','test');

% To run AFQ using mrtrix for tractography
afq = AFQ_Create('sub_dirs',dt6dirs,...
                 'sub_group',[0 0],...
                 'computeCSD',1);
afq = AFQ_run([],[],afq);


% Find the VOF per every subject
wholebrainfgPath= '/Users/gari/Documents/BCBL_PROJECTS/MINI/ANALYSIS/DWI/S011/dmri60/dti60trilin/fibers'; 
fgMori = dtiReadFibers('MoriGroups.mat')
L_arcuate= fgMori(19);
R_arcuate= fgMori(20);
% Create lables from freesurfer. But this is has been ac-pc-ed, I think I should
% ac-pc the aparc+aseg as well.
  fsIn   = '/path/to/aparc+aseg.mgz';
  outDir = '/save/directory/rois';
  type   = 'mat';
  refT1  = '/path/to/t1Anatomical.nii.gz';
  fs_roisFromAllLabels(fsIn,outDir,type,refT1);fsROIdir= ;
outdir= ;
thresh= ;
v_crit= ;
dt= ;
savefiles= ;
arcThresh= ;
parcThresh= ;
% [L_VOF, R_VOF, L_pArc, R_pArc, L_pArc_vot, R_pArc_vot] = ...
%                                 AFQ_FindVOF(wholebrainfgPath,...
%                                             L_arcuate,...
%                                             R_arcuate,...
%                                             fsROIdir,...
%                                             outdir,...
%                                             thresh,...
%                                             v_crit, ...
%                                             dt, ...
%                                             savefiles, ...
%                                             arcThresh, ...
%                                             parcThresh)


% TO DO: integrate parallel version:
% afq = AFQ_run_sge_LH(afq, 2, 3); %tmp

save(fullfile(basedir, 'afqOut'), 'afq')

