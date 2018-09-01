function batch_fslpreprocessdiffusion(projectsDir, projectName, subName, shell, ...
                                      doPreProc, doBias,...
                                      doDtiInit, doAfqCreate, ...
                                      doAfqRun, doCreateProfiles)
dataWhere = 'HCPblack';  % 'HCPbcbl', 'HCPblack', 'MINIbcbl'                                      
%% GLU MINI project adapted from: BDE lab preprocessing for diffusion data
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
% 
% subName = 'S041'
% % AnalysisDir = '/bcbl/home/public/Gari/MINI/ANALYSIS'
% % AnalysisDir = '/Users/gari/Documents/BCBL_PROJECTS/MINI/ANALYSIS'
% AnalysisDir = '/bcbl/home/home_g-m/glerma/00local/PROYECTOS/dr/ANALYSIS'
% shell = '1000'
%
%{
subName = 'SC0131'
projectName = 'ILLITERATE'
projectsDir = '/bcbl/home/public/Gari'
shell = '1000'

doPreProc = 1
doBias   = 0
doDtiInit = 0
doAfqCreate = 0
doAfqRun = 0
doCreateProfiles = 0
%}
% Edited by GLU on June 2016
% update instructions



%% Set up directories and find files
projectDir  = fullfile(projectsDir, projectName);
DataDir     = fullfile(projectDir, 'DATA');
dicomDir    = fullfile(DataDir, 'dicoms');
AnalysisDir = fullfile(projectDir, 'ANALYSIS');
DWIdir  = fullfile(AnalysisDir, 'DWI');
retdir  = fullfile(AnalysisDir, 'ret');

fsp = filesep;
qMRIdir = fullfile(AnalysisDir, 'qMRI_acpc');
qMRIsubPATH = fullfile(qMRIdir, subName, 'OutPutFiles_1','BrainMaps');
basedir = fullfile(DWIdir,subName);
% rawdir = fullfile(DWIdir,subName, 'raw');
rawdir = fullfile(dicomDir, subName);
if ~exist(rawdir, 'dir'), mkdir(rawdir), end

%%% EDIT %%%
% t1dir = fullfile(retdir, subName, 'anat');

dmridir = fullfile(basedir,['dmri' num2str(shell)]);  % 'mrtrix3']);
% copyfile(fullfile(basedir,'t1_std_acpc.nii.gz'), ...
%          fullfile(dmridir,'t1_std_acpc.nii.gz')); % Otherwise it won't find it downstream, since in the dt6.mat files.t1 only the name is stores
% t1 = fullfile(dmridir, 't1_std_acpc.nii.gz');
% FSdir   = fullfile(AnalysisDir, 'freesurferacpc');
FSdir   = fullfile(AnalysisDir, 'freesurfer');
aparcAsegDir = fullfile(FSdir, subName, 'mri');
% t1        = fullfile(aparcAsegDir, 'T1.nii.gz');
t1        = fullfile(dicomDir, subName, 'T1w.nii');
t1dir   = aparcAsegDir;
aparcAseg = fullfile(aparcAsegDir, 'aparc+aseg.mgz');
% copyfile(aparcAseg, dmridir);
%%%      %%%

if ~exist(dmridir)
    mkdir(dmridir)
end

% Note glu: conversion from dicom in the ipython notebook file with qsubs
switch projectName
    case {'MINI'}
        d30 = dir(fullfile(rawdir,'*d35b1000*.nii'));
        d60 = dir(fullfile(rawdir,'*d65b2500*.nii'));
        b0 = dir(fullfile(rawdir,'*d6b0*.nii')); % grab post-anterior encoded file 
        dMRIFiles = {};
        % The dwell time is exactly the same of the example = 0.095 
        % (Calculated from echo spacing 0.75 and EPI factor 128)
        % Phase encode matrix. This denotes, for each volume, which direction is
        % the phase encode
        % Edit GLU: See here the info used for MINI data in Siemens TRIO
        % http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/TOPUP/Faq
        % We have 5xb0 + 30 dirs & 5xb0 + 60 dirs with one A>>P (this is y=-1)
        % and 6xb0 with P>>A (this is y=1)
        switch shell
            case {'1000'}  % Use singleShell b = 1000
                dMRIFiles{1}=fullfile(rawdir,b0(1).name);
                dMRIFiles{2}=fullfile(rawdir,d30(1).name);
                pe_mat = [0 1 0; 0 -1 0];
                dwellTime = 0.095;
                runTopup = 1;
            case {'2500'}  % Use singleShell b = 2500
                dMRIFiles{1}=fullfile(rawdir,b0(1).name);
                dMRIFiles{2}=fullfile(rawdir,d60(1).name);
                pe_mat = [0 1 0; 0 -1 0];
                dwellTime = 0.095;            
                runTopup = 1;
            case {'MS'}  % Use multishell
                dMRIFiles{1}=fullfile(rawdir,b0(1).name);
                dMRIFiles{2}=fullfile(rawdir,d30(1).name);
                dMRIFiles{3}=fullfile(rawdir,d60(1).name);
                pe_mat = [0 1 0; 0 -1 0; 0 -1 0];
                dwellTime = 0.095; 
                runTopup = 1;
            otherwise
                error('Unknown shell, valid values are MS, 1000 and 2500')
        end
        % Bvals and Bvecs files
        for ii = 1:length(dMRIFiles)
            bvals{ii} = [prefix(prefix(dMRIFiles{ii})) '_bvals'];
            bvecs{ii} = [prefix(prefix(dMRIFiles{ii})) '_bvecs'];
        end
    case {'dr'}
        disp('to be checked, update this file with info in black')
    case {'ILLITERATE'}
        % d = dir(fullfile(rawdir,'DWI.nii'));
        d = dir(fullfile(dicomDir, subName, 'DWI.nii'));
        dMRIFiles = {};
        switch shell
            case {'1000'}  % There is only one shell
                dMRIFiles{1}=fullfile(rawdir,d(1).name);
                pe_mat = []; % we don't have reverse encoded file
                dwellTime = [];  
                runTopup = 0;
            otherwise
                error('Unknown shell, valid value is 1000')
        end
        % Bvals and Bvecs files
        for ii = 1:length(dMRIFiles)
            % bvals{ii} = [prefix(prefix(dMRIFiles{ii})) '.bval'];
            % bvecs{ii} = [prefix(prefix(dMRIFiles{ii})) '.bvec'];
            bvals{ii} = fullfile(rawdir, 'bval');
            bvecs{ii} = fullfile(rawdir, 'bvec');
        end
    otherwise
        error('This project has not been created yet')

end

% Directory to save everything
outdir = fullfile(dmridir);
% mrtrixdiir = fullfile(outdir, ['dti' shell 'trilin'],'mrtrix');

%% doPreProc: This is mostly done with command line calls to FSL
if doPreProc
    fsl_preprocess(dMRIFiles, bvecs, bvals, pe_mat, outdir, ...
                   dwellTime, shell, runTopup);
end

%% doBias: anadido a posteriori como parte de mrTrix quantitative analysis para obtener AFD
% ESTO NO ESTA HECHO EN MINI, LO HICE EN UNO EN LOCAL PARA VER QUE FUNCIONA
% DWI bias field correction is perfomed by first estimating a correction field 
% from the DWI b=0 image, then applying the field to correct all DW volumes. 
% This can be done in a single step using the dwibiascorrect script in MRtrix. 
% The script uses bias field correction algorthims available in ANTS or FSL. 
% In our experience the N4 algorithm in ANTS gives superiour results. 
% To install N4 install the ANTS package, then run perform bias field correction 
% on DW images using:

% dmridir = '/Users/gari/Documents/BCBL_PROJECTS/MINI/ANALYSIS/DWI/S005/dmri'

% Para esto hace falta mrtrix y ANTS instalados
if doBias
    bvecs = fullfile(dmridir,'eddy','bvecs');
    bvals = fullfile(dmridir,'eddy','bvals');
    input_brain_mask = fullfile(dmridir,'eddy','nodif_brain_mask.nii.gz');
    input_dwi = fullfile(dmridir,'eddy','data.nii.gz');
    output_corrected_dwi = fullfile(dmridir,'eddy','biasdata.nii.gz');
    cmd = ['dwibiascorrect -ants ' ...
           '-fslgrad ' bvecs ' ' bvals ' ' ...
           '-mask ' input_brain_mask ' ' ...
            input_dwi ' ' ...
            output_corrected_dwi];   
    system(cmd);
end

%% Global intensity normalisation across subjects
% ESTO NO ESTA HECHO EN MINI, LO HICE EN UNO EN LOCAL PARA VER QUE FUNCIONA
% Ahora hay que hacer esto pero en paralelo para todos los sujetos, o sea que
% hay que llamar primero al batch con las opciones de arriba, hacer esto, y
% luego seguir de doDtiInit para abajo. 

% Copiar todos los archivos a una carpeta. Pero habra que cambiar el nombre
    % para poder identificarlos luego. 
% if(0) % Generar el comando pero lanzarlo por command line
%     DWIdir = '/bcbl/home/public/Gari/MINI/ANALYSIS/DWI';
%     DWI2dir = '/export/home/glerma/glerma/00local/PROYECTOS/MINI/ANALYSIS/DWI';
%     input_dwi_folder = [DWI2dir fsp 'input_dwi_folder'];
%     input_brain_mask_folder = [DWI2dir fsp 'input_brain_mask_folder'];
%     output_normalised_dwi_folder = [DWI2dir fsp 'output_normalised_dwi_folder'];
%     output_fa_template = [DWI2dir fsp 'output_fa_template'];
%     output_template_wm_mask = [DWI2dir fsp 'output_template_wm_mask'];
%     tempdir = '/scratch/glerma';
%     cmd = ['dwiintensitynorm -force -nocleanup -nthreads 16 ' ...
%            '-tempdir ' tempdir ' ' ...
%             input_dwi_folder ' ' input_brain_mask_folder ' ' ...
%             output_normalised_dwi_folder ' ' output_fa_template ' ' ...
%             output_template_wm_mask]
%     % Los dos ultimos son mriconvert de output fa template y template brain
%     % mask. Pide por archivo, no folder. Lo he hecho a mano, y he tanto
%     % .mif como .nii.gz
%         
% end

%% Computing a group average response function
% idem arriba
% average_response <input_response_files (mulitple inputs accepted)> 
%                  <output_group_average_response>
% if (0)
%     wmAll = [];
%     csAll = [];
%     gmAll = [];
%     gm_output_group_average_response = [DWI2dir fsp 'output_group_average_response' fsp ...
%                                         'gm_output_group_average_response.txt'];
%     wm_output_group_average_response = [DWI2dir fsp 'output_group_average_response' fsp ...
%                                         'wm_output_group_average_response.txt'];
%     cs_output_group_average_response = [DWI2dir fsp 'output_group_average_response' fsp ...
%                                         'cs_output_group_average_response.txt'];                               
%     for ns = 1 : length(subs)
%         subname = subs(ns).name
%         % Folders
%         DWIdir = '/bcbl/home/public/Gari/MINI/ANALYSIS/DWI';
%         dmridir = fullfile(DWIdir, subname, 'dmri*');
%         mrtrixdir = fullfile(dmridir, 'noNorm_dti90trilin','mrtrix');
%         wmResponse = fullfile(mrtrixdir, 'data_aligned_trilin_noMEC_wmResponse.txt');
%         gmResponse = fullfile(mrtrixdir, 'data_aligned_trilin_noMEC_gmResponse.txt');
%         csResponse = fullfile(mrtrixdir, 'data_aligned_trilin_noMEC_csfResponse.txt');
%         wmAll = [wmResponse ' ' wmAll]; 
%         csAll = [csResponse ' ' csAll]; 
%         gmAll = [gmResponse ' ' gmAll]; 
%     end
%     wm_cmd = ['average_response ' wmAll ' ' wm_output_group_average_response];
%     gm_cmd = ['average_response ' gmAll ' ' gm_output_group_average_response];
%     cs_cmd = ['average_response ' csAll ' ' cs_output_group_average_response];
%     system(wm_cmd)
%     system(gm_cmd)
%     system(cs_cmd)
% end

%% doDtiInit: to fit tensor model (ACTUALIZAR 'dtEddy' CON LA SALIDA DE LO ANTERIOR
% Turn off motion and eddy current correction, that was taken care of by FSL

% GLU: I shouldn't fit the tensor model since I have the multishell data

params = dtiInitParams; % Set up parameters for controlling dtiInit
switch projectName
    case {'MINI'}
        params.eddyCorrect=-1; % This turns off eddy current and motion correction
        params.rotateBvecsWithCanXform=1; % Siemens data requires this to be 1
        params.phaseEncodeDir=2; % AP phase encode, 1 is for R>>L (2 = A/P 'col')
        params.clobber=1; % Overwrite anything previously done
        params.fitMethod='ls'; % 'ls, or 'rt' for robust tensor fitting (longer)
        params.dt6BaseName = ''
        dtEddy = fullfile(dmridir,'eddy','data.nii.gz'); % Path to the data
        params.bvalsFile = fullfile(dmridir,'eddy','bvals'); % Path to bvals
        params.bvecsFile = fullfile(dmridir,'eddy','bvecs'); % Path to the bvecs
    case  {'ILLITERATE'}
        params.eddyCorrect=1; % In this case with no reversed PE, no topup, no fslpreproc, all in mrDiff
        params.rotateBvecsWithCanXform=1; % Philips data requires this to be 1
        params.phaseEncodeDir=2;  % AP phase encode
        params.clobber=1; % Overwrite anything previously done
        params.fitMethod='ls'; % 'ls, or 'rt' for robust tensor fitting (longer)
        params.dt6BaseName = '';
        mkdir([dmridir filesep 'eddy'])
        copyfile(dMRIFiles{1}, [dmridir filesep 'eddy'])
        copyfile(bvals{1}    , [dmridir filesep 'eddy'])
        copyfile(bvecs{1}    , [dmridir filesep 'eddy'])
        dtEddy = fullfile(dmridir,'eddy','dwi.nii.gz'); % Path to the data
        params.bvalsFile = fullfile(dmridir,'eddy','dwi.bval'); % Path to bvals
        params.bvecsFile = fullfile(dmridir,'eddy','dwi.bvec'); % Path to the bvecs
    otherwise
        error('Project unknown')
end

if doDtiInit
    dt6FileName = dtiInit(dtEddy,t1,params); % Run dtiInit to preprocess data
else
    % dtiDir = dir(fullfile(dmridir,'dti*trilin*'));  % for MINI
    % if dtiDir.isdir && exist(fullfile(dmridir, dtiDir.name,'dt6.mat'),'file' )
    %     dt6FileName = {fullfile(dmridir, dtiDir.name,'dt6.mat')};
    % else
    %     error('Cannot find dt6.mat file')
    % end
    dt6FileName = {fullfile(dmridir,'dti','dt6.mat')};    
end

%% doAfqCreate
% Cell array with paths to the dt6 directories
dt6dirs = horzcat({fileparts(dt6FileName{1})});

% afq = AFQ_Create('sub_dirs',dt6dirs,'sub_group',[0 0],'clip2rois', 0);
% To run AFQ in test mode so it will go quickly
% afq = AFQ_Create('sub_dirs',dt6dirs,'sub_group',[0 0],'run_mode','test');


% To run AFQ using mrtrix for tractography
if doAfqCreate
    cd(dmridir)  
    % Delete old version if exists
    if exist(fullfile(dmridir, [subName '_b' shell '_afqOutAfqCreate.mat']))
        delete(fullfile(dmridir, [subName '_b' shell '_afqOutAfqCreate.mat']))
    end
    % Create a new one
    disp(['\n\n This is the dt6dirs{1}: ' dt6dirs{1}])
    afq = AFQ_Create('sub_dirs',dt6dirs,...
                     'sub_group', [0], ...
                     'sub_names', [subName],...
                     'computeCSD',1);
    switch dataWhere
        case {'MINIbcbl'}
            % FOR MINI: Create paths to qMRI images. Only 1 subject, so it is always 1.
            % Removing it for HCP in blackn 
            % t1Path{1}  = fullfile(qMRIsubPATH, 'T1_map_Wlin.nii.gz');
            % mtvPath{1} = fullfile(qMRIsubPATH, 'TV_map.nii.gz');
            % afq = AFQ_set(afq, 'images', t1Path);
            % afq = AFQ_set(afq, 'images', mtvPath);
        otherwise
            disp('Do not add qMRI files to create profiles')
    end
    save(fullfile(dmridir, [subName '_b' shell '_afqOutAfqCreate']), 'afq')
end

%% doAfqRun
if doAfqRun
    load(fullfile(dmridir, [subName '_b' shell '_afqOutAfqCreate.mat']))
    % Now add the new FA calculated by mrtrix and create a new measurement
    % This depends on the project, again
    switch projectName
        case {'MINI'} 
            % Create paths to qMRI images. Only 1 subject, so it is always 1.
            t1Path{1}  = fullfile(qMRIsubPATH, 'T1_map_Wlin.nii.gz');
            mtvPath{1} = fullfile(qMRIsubPATH, 'TV_map.nii.gz');
            if strcmp(shell,'1000'); ndirs='30';end;
            if strcmp(shell,'2500'); ndirs='60';end;
            mrtrixPath = fullfile(dmridir, ['dti' ndirs 'trilin'], 'mrtrix');
            faPath{1}  = fullfile(mrtrixPath, 'data_aligned_trilin_noMEC_fa.nii.gz');
            
            afq = AFQ_set(afq, 'images', t1Path);
            afq = AFQ_set(afq, 'images', mtvPath);
            afq = AFQ_set(afq, 'images', faPath);
        case {'ILLITERATE'} 
            ndirs = '64';
            mrtrixPath = fullfile(dmridir, ['dti' ndirs 'trilin'], 'mrtrix');
            famif = fullfile(mrtrixPath, 'dwi_aligned_trilin_fa.mif');
            fanii = fullfile(mrtrixPath, 'dwi_aligned_trilin_fa.nii.gz');
            AFQ_mrtrix_mrconvert(famif, fanii, 0, 0, '3');
            faPath{1}  = fanii;
            afq = AFQ_set(afq, 'images', faPath);
        case {'MINIbcbl'}
            if strcmp(shell,'1000'); ndirs='30';end;
            if strcmp(shell,'2500'); ndirs='60';end;
            mrtrixPath = fullfile(dmridir, ['dti' ndirs 'trilin'], 'mrtrix');
            faPath{1}  = fullfile(mrtrixPath, 'data_aligned_trilin_noMEC_fa.nii.gz');
            afq = AFQ_set(afq, 'images', faPath);
            afq = AFQ_run([],[],afq);
            save(fullfile(dmridir, [subName '_b' shell '_afqOutAfqRun']), 'afq')
        case {'HCPbcbl', 'HCPblack'}
            mrtrixPath = fullfile(dmridir, 'dti', 'mrtrix');
            faPath{1}  = fullfile(mrtrixPath, 'dwi_aligned_trilin_noMEC_fa.nii.gz')
            afq = AFQ_set(afq, 'images', faPath);
            % Dont do tracking, I would have done it separately and then run LiFE before this
            afq = AFQ_run([],[],afq);
            save(fullfile(dmridir, [subName '_b' shell '_afqOutAfqRun']), 'afq')
        otherwise
            disp('Unknown dataWhere')
    end
end

%% doCreateProfiles
% if doCreateProfiles
%    disp('todo')
% end
    




