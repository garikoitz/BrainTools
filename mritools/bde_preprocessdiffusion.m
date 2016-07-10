function bde_preprocessdiffusion(subName, DWIdir, retdir)
% GLU MINI project adapted from: BDE lab preprocessing for diffusion data
% This scripts asumes that the there nifti anatomical files and that all
% the preprocessing (topup, eddy) has been done already.
%
% Inputs
% subName     - Just the name of the subject with ' '
% DWIdir      - path to the base ANALYSIS/DWI directory, before the subjects
% retdir      - path to the ac-pc aligned t1. Everything will be
%               coregistered to this. In my case they are in the ret
%               directory since I had to do it there first. 
%
% 
% Edited by GLU on June 2016
% Eliminated makeNifti and doPreproc since it is being taking care outside
% - makeNifti: in the iPython notebook with qsubs calls to mcverter
% - dofslPreproc: notebook qsub calls to
% mritools/batch_fslpreprocessdiffusion.m, which has been adapted from
% here. 


% I moved everything that was single subject to batch_fslpreprocessdiffusion,
% now here I just have to do the AFQ create but for a large set of subjects. 




%% Set up directories and find files
basedir = fullfile(DWIdir, subName)
rawdir = fullfile(basedir,'raw');
dmridir = fullfile(basedir,'dmri');
t1dir = fullfile(retdir, subName, 'anat');
copiedt1 = fullfile(dmridir, 't1_std_acpc.nii.gz');



%% Run AFQ

% Cell array with paths to the dt6 directories
% % dt6dirs = horzcat(fileparts(dt6FileName{1}), fileparts(dt6FileName{2}));


% dt6dirs = horzcat({fileparts(dt6FileName{1}{1})}, {fileparts(dt6FileName{2}{1})});
% dt6dirs = horzcat({'/Users/gari/Documents/BCBL_PROJECTS/MINI/ANALYSIS/DWI/S011/dmri60/dti60trilin'},{'/Users/gari/Documents/BCBL_PROJECTS/MINI/ANALYSIS/DWI/S011/dmri30/dti30trilin'})
dt6dirs = horzcat({fileparts(dt6FileName{1})});
%dt6dirs = horzcat({'/bcbl/home/public/Gari/MINI/ANALYSIS/DWI/S002/dmri/dti90trilin'});
%dt6dirs = horzcat({'/Users/gari/Documents/BCBL_PROJECTS/MINI/ANALYSIS/DWI/S011/dmri/dti90trilin'});

% afq = AFQ_Create('sub_dirs',dt6dirs,'sub_group',[0 0],'clip2rois', 0);
% To run AFQ in test mode so it will go quickly
% afq = AFQ_Create('sub_dirs',dt6dirs,'sub_group',[0 0],'run_mode','test');

% To run AFQ using mrtrix for tractography
afq = AFQ_Create('sub_dirs',dt6dirs,...
                 'sub_group', [0], ...
                 'sub_names', ['S011'],...
                 'computeCSD',1);
afq = AFQ_run([],[],afq);
save(fullfile(basedir, 'afqOut'), 'afq')


%% Find the VOF per every subject

% Create lables from freesurfer. But this is has been ac-pc-ed, I think I should
% ac-pc the aparc+aseg as well. Use the same code used for ribbon
% TODO: fix it for every subject and for the cluster
fs_SUBJECTS_DIR = '/bcbl/home/public/Gari/MINI/ANALYSIS/freesurfer';
subName = 'S011';
path2anat = '/bcbl/home/public/Gari/MINI/ANALYSIS/ret/S011/anat';
subjID = fullfile(fs_SUBJECTS_DIR, subName, 'mri', 'aparc+aseg.mgz');
subjID2009 = fullfile(fs_SUBJECTS_DIR, subName, 'mri', 'aparc.a2009s+aseg.mgz');
outfile     = fullfile(path2anat, 't1_aparcaseg.nii.gz');
outfile2009     = fullfile(path2anat, 't1_aparcaseg2009.nii.gz');
fillWithCSF = true;
alignTo     = fullfile(path2anat, 't1.nii.gz');
resample_type = [];
system(['mri_convert  --out_orientation RAS --reslice_like ' alignTo ...
         ' ' subjID ' ' outfile]);
system(['mri_convert  --out_orientation RAS --reslice_like ' alignTo ...
         ' ' subjID2009 ' ' outfile2009]);
% And now do the ac-pc using the xform we already had
aparc = matchfiles(fullfile(path2anat, 't1_aparcaseg.nii.gz'));
aparc2009 = matchfiles(fullfile(path2anat, 't1_aparcaseg2009.nii.gz'));
acpcMatrix = load(fullfile(path2anat, 'xform2acpc.mat'));
NameOfAcpcFile = 't1_class_std_acpc';
acpcOutfile     = fullfile(path2anat, 't1_aparcaseg_std_acpc.nii.gz');
acpcOutfile2009     = fullfile(path2anat, 't1_aparcaseg2009_std_acpc.nii.gz');
T1acpc = fullfile(path2anat, 't1_std_acpc.nii.gz');
mrAnatAverageAcpcNifti(aparc, acpcOutfile, acpcMatrix.alignLandmarks);
mrAnatAverageAcpcNifti(aparc2009, acpcOutfile2009, acpcMatrix.alignLandmarks);
close all;
% It is not working, the rotation changes the labels names, talk to Eugenio
path2anat= '/Users/gari/Documents/BCBL_PROJECTS/MINI/ANALYSIS/ret/S011/anat';
fsIn   = fullfile(path2anat, 'aparc+aseg.mgz');
outDir = fullfile(path2anat,'aparcRoi');
if ~exist(outDir, 'dir'), mkdir(outDir), end
type   = 'nifti';
refT1  = fullfile(path2anat, 'T1.mgz');
fs_roisFromAllLabels(fsIn,outDir,type,refT1);

wholebrainfgPath= '/Users/gari/Documents/BCBL_PROJECTS/MINI/ANALYSIS/DWI/S011/dmri60/dti60trilin/fibers/WholeBrainFG.mat'; 
fgMori = dtiReadFibers(fullfile(wholebrainfgPath, 'MoriGroups.mat'));
L_arcuate= fgMori(19);
R_arcuate= fgMori(20);
fsROIdir= outDir;
outdir = fullfile(wholebrainfgPath,'VOF');
if ~exist(outdir, 'dir'), mkdir(outdir), end
thresh= [];
v_crit= [];
dt= dtiLoadDt6('/Users/gari/Documents/BCBL_PROJECTS/MINI/ANALYSIS/DWI/S011/dmri60/dti60trilin/dt6.mat');
savefiles= true;
arcThresh= [];
parcThresh= [];
[L_VOF, R_VOF, L_pArc, R_pArc, L_pArc_vot, R_pArc_vot] = ...
                                AFQ_FindVOF(wholebrainfgPath,...
                                            L_arcuate,...
                                            R_arcuate,...
                                            fsROIdir,...
                                            outdir,...
                                            thresh,...
                                            v_crit, ...
                                            dt, ...
                                            savefiles, ...
                                            arcThresh, ...
                                            parcThresh)


% TO DO: integrate parallel version:
% afq = AFQ_run_sge_LH(afq, 2, 3); %tmp



