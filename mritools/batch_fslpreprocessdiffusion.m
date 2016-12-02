function batch_fslpreprocessdiffusion(subName, AnalysisDir, doPreProc, doDtiInit, doAfqCreate, doAfqRun)
% GLU MINI project adapted from: BDE lab preprocessing for diffusion data
%
%
% Inputs
% basedir     - path to the base directory with the raw data
% t1dir       - path to the ac-pc aligned t1. Everything will be
%               coregistered to this
% doPreProc   - logical. Whether or not to do preprocessing. Maybe it failed in
% a previous step and we want to start right after it. 
%
% Example:
% 
% AnalysisDir = '/bcbl/home/public/Gari/MINI/ANALYSIS'
% % AnalysisDir = '/Users/gari/Documents/BCBL_PROJECTS/MINI/ANALYSIS'
% subName = 'S002'
% doPreProc = 0
% doDtiInit = 0
% doAfqCreate = 1
% doAfqRun = 1
%
% Edited by GLU on June 2016
% update instructions



%% Set up directories and find files
DWIdir = fullfile(AnalysisDir, 'DWI');
retdir = fullfile(AnalysisDir, 'ret');
FSdir = fullfile(AnalysisDir, 'freesurferacpc');
basedir = fullfile(DWIdir,subName);
rawdir = fullfile(DWIdir,subName, 'raw');
if ~exist(rawdir, 'dir'), mkdir(rawdir), end
t1dir = fullfile(retdir, subName, 'anat');
dmridir = fullfile(basedir,'dmri');
aparcAsegDir = fullfile(FSdir, subName, 'mri');

% Note glu: conversion from dicom in the ipython notebook file with qsubs
d30 = dir(fullfile(rawdir,'*d35b1000*.nii'));
d60 = dir(fullfile(rawdir,'*d65b2500*.nii'));
b0 = dir(fullfile(rawdir,'*d6b0*.nii')); % grab post-anterior encoded file 


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
outdir = fullfile(basedir,'dmri');


%% Pre process: This is mostly done with command line calls to FSL
if doPreProc
    fsl_preprocess(dMRIFiles, bvecs, bvals, pe_mat, outdir, dwellTime);
end
%% Run dtiInit to fit tensor model
% Turn off motion and eddy current correction, that was taken care of by FSL

% GLU: I shouldn't fit the tensor model since I have the multishell data

% OJO: decidir sobre el acpc o usar el de FS!

% Set up t1 path and the params that are common
% t1 = fullfile(t1dir,'t1_std_acpc.nii.gz'); % Path to the acpc t1-weighted image
t1 = fullfile(t1dir,'t1_std_acpc.nii.gz');
copyfile(t1, dmridir); % Otherwise it won't find it downstream, since in the dt6.mat files.t1 only the name is stores

aparcAseg = fullfile(aparcAsegDir, 'aparc+aseg.mgz');
copyfile(aparcAseg, dmridir);

params = dtiInitParams; % Set up parameters for controlling dtiInit
params.eddyCorrect=-1; % This turns off eddy current and motion correction
params.rotateBvecsWithCanXform=1; % Siemens data requires this to be 1
params.phaseEncodeDir=2; % AP phase encode, 1 is for R>>L (2 = A/P 'col')
params.clobber=1; % Overwrite anything previously done
params.fitMethod='ls'; % 'ls, or 'rt' for robust tensor fitting (longer)

dtEddy = fullfile(dmridir,'eddy','data.nii.gz'); % Path to the data
params.bvalsFile = fullfile(dmridir,'eddy','bvals'); % Path to bvals
params.bvecsFile = fullfile(dmridir,'eddy','bvecs'); % Path to the bvecs

if doDtiInit
    dt6FileName = dtiInit(dtEddy,t1,params); % Run dtiInit to preprocess data
else
    dtiDir = dir(fullfile(dmridir,'dti*trilin*'));
    if dtiDir.isdir && exist(fullfile(dmridir, dtiDir.name,'dt6.mat'),'file' )
        dt6FileName = {fullfile(dmridir, dtiDir.name,'dt6.mat')};
    else
        error('Cannot find dt6.mat file')
    end
    
end
%% Run AFQ

% Cell array with paths to the dt6 directories
dt6dirs = horzcat({fileparts(dt6FileName{1})});
%dt6dirs = horzcat({'/bcbl/home/public/Gari/MINI/ANALYSIS/DWI/S002/dmri/dti90trilin'});
%dt6dirs = horzcat({'/Users/gari/Documents/BCBL_PROJECTS/MINI/ANALYSIS/DWI/S011/dmri/dti90trilin'});

% afq = AFQ_Create('sub_dirs',dt6dirs,'sub_group',[0 0],'clip2rois', 0);
% To run AFQ in test mode so it will go quickly
% afq = AFQ_Create('sub_dirs',dt6dirs,'sub_group',[0 0],'run_mode','test');

% To run AFQ using mrtrix for tractography
if doAfqCreate
    afq = AFQ_Create('sub_dirs',dt6dirs,...
                     'sub_group', [0], ...
                     'sub_names', [subName],...
                     'computeCSD',1);
    save(fullfile(basedir, 'afqOut'), 'afq')
end
if doAfqRun
    load(fullfile(basedir, 'afqOut.mat'))
    afq = AFQ_run([],[],afq);
    save(fullfile(basedir, 'afqOut'), 'afq')
end



