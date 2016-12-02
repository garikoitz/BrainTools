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


                             
                                        
                                        
% TO DO: integrate parallel version:
% afq = AFQ_run_sge_LH(afq, 2, 3); %tmp





%% Lanzar los sujetos que me faltan para el batch_fslpre...  con un parfor a ver si se quita la mierda del Image_Toolbox license error
OnlyDoAfqRun = {'S095','S050','S069','S010','S016','S032','S087','S018','S075','S057','S055','S024','S071','S034','S059','S077','S091','S035','S090','S079','S005','S060','S041'};
% myPool = parpool(length(OnlyDoAfqRun))
% addAttachedFiles(myPool, ...
%                  '/opt/matlab/R2014b/toolbox/matlab/lang/@char/exist');
A = {'S034','S035','S041','S055','S057','S059','S060','S071','S075','S077','S079','S090','S091','S099'}
A = {'S024'};
parfor ii =1:length(A)
    AnalysisDir = '/bcbl/home/public/Gari/MINI/ANALYSIS';
    subName = A{ii};
    doPreProc = 0;
    doDtiInit = 0;
    doAfqCreate = 0;
    doAfqRun = 1;
    batch_fslpreprocessdiffusion(subName, AnalysisDir, doPreProc, ...
                                 doDtiInit, doAfqCreate, doAfqRun);
end



