%% Desde tractografia a ROI

%% Inicializar
% Despues de haber hecho el ROI analisis en funcional, mi objetivo en este
% caso es el:
% 1.- Encontrar las fibras del VOF y del poserior arcuate y del arcuate
% 2.- HAcer conteo de fibras, caracteristicas, y crear ROIs de los tractos
% 3.- Crear ROI-s individuales de donde estan llegando estos tractos

clear all; close all; 
fsp = filesep;

% Folder Names
% LOCAL
%{
MINIDIR = '/Users/gari/Documents/BCBL_PROJECTS/MINI';
fsbin = '/Applications/freesurfer6/bin';
fshome = '/Applications/freesurfer6'; 
%}

% SERVER
%
MINIDIR = '/bcbl/home/public/Gari/MINI';
fsbin = '/opt/freesurfer-5.3.0/freesurfer/bin';
fshome = '/opt/freesurfer-5.3.0/freesurfer'; 
%}



AnalysisDir = [MINIDIR fsp 'ANALYSIS'];
fs_SUBJECTS_DIR = fullfile(AnalysisDir, 'freesurferacpc');
qmridir = fullfile(AnalysisDir, 'qMRI_acpc');
DWIdir  = fullfile(AnalysisDir, 'DWI');
cd(DWIdir);
subs = dir('S*');
retDIR = fullfile(AnalysisDir, 'ret');
fMRIDIR = fullfile(AnalysisDir, 'fMRI_SPM', 'block', 'data');

      





%%%%%   CLUSTER PARPOOL    %%%%%%
% myclusterLocal = parcluster('local');
% myclusterLocal.NumWorkers
%{
[st, re] = system('qstat -g c | grep matlab.q');
[Tok, Rem] = strtok(re);
[Tok, Rem] = strtok(Rem);
[Tok, Rem] = strtok(Rem);
[Tok, Rem] = strtok(Rem);
[available] = strtok(Rem)
parpool('ips_base', str2num(available))
%}
%%%%% END CLUSTER PARPOOL  %%%%%%

%% Encontrar VOF PARC y pasarlos a surfaces en fsaverage
if(0)
  for ns = 1 : length(subs)
    %% Find the VOF per every subject
    subname = subs(ns).name

    setenv('FREESURFER_HOME', fshome); 


    path2anat = fullfile(retDIR, subname, 'anat');
    path2fMRIanat = fullfile(fMRIDIR, subname, 'anat');
      dmridir = fullfile(DWIdir, subname, 'dmri');
    cd(dmridir)




    % Create labels from freesurfer for every subject
    fsIn   = fullfile(dmridir, 'aparc+aseg.mgz');
    outDir = fullfile(dmridir,'ROIs');
    if ~exist(outDir, 'dir'), mkdir(outDir), end
    type   = 'mat';
    refT1  = fullfile(dmridir, 't1_std_acpc.nii.gz');

  
%   fs_roisFromAllLabels(fsIn,outDir,type,refT1);

    % %%%%%%%%% NOTE  %%%%%%%%%%%%    
    %{
    'dti90trilin' has been edited to 'noNorm_dti90trilin'.
    The noNorm_ version is the original one. Afterwards I created the 
    the normalized version because it was a requirement for mrtrix's SIFT2
    I tried to reutilize all the files that where useful. 
    Be careful with the not normalized version, as there was no space in
    public, I had to create symbolic links to my folder and store the new
    analyses there. It can be a little bit messy.  
    
    
    %}
    % %%%%%%%%% NOTE  %%%%%%%%%%%%        



    % wholebrainfgPath= fullfile(dmridir, 'dti90trilin', 'fibers');
    wholebrainfgPath= fullfile(dmridir, 'noNorm_dti90trilin', 'fibers'); 
    wholebrainfg= fullfile(wholebrainfgPath, 'WholeBrainFG.mat'); 
    fgMori = dtiReadFibers(fullfile(wholebrainfgPath, 'MoriGroups.mat'));
    L_arcuate= fgMori(19);
    R_arcuate= fgMori(20);
    fsROIdir= outDir;
    outdir = fullfile(wholebrainfgPath,'VOF');
    if ~exist(outdir, 'dir'), mkdir(outdir), end
    thresh= [];
    v_crit= [];
    % GLU dt= dtiLoadDt6(fullfile(dmridir, 'dti90trilin', 'dt6.mat'));
    % load(fullfile(dmridir, 'dti90trilin', 'dt6.mat'));
    % dt.dataFile = fullfile(dmridir, 'dti90trilin');
    load(fullfile(dmridir, 'noNorm_dti90trilin', 'dt6.mat'));
    dt.dataFile = fullfile(dmridir, 'noNorm_dti90trilin');
    savefiles= true;
    arcThresh= [];
    parcThresh= [];
    
    % LEave this commented, I want to obtain only ARCUATE
    %{
    % Obtain the tracts of interest
    [L_VOF, R_VOF, L_pArc, R_pArc, L_pArc_vot, R_pArc_vot] = ...
                                    AFQ_FindVOF(wholebrainfg,...
                                                L_arcuate,...
                                                R_arcuate,...
                                                fsROIdir,...
                                                outdir,...
                                                thresh,...
                                                v_crit, ...
                                                dt, ...
                                                savefiles, ...
                                                arcThresh, ...
                                                parcThresh);    
    % % and save them
    save(fullfile(outdir, 'VOF_all.mat'), ...
        'L_VOF', 'R_VOF', 'L_pArc', 'R_pArc', 'L_pArc_vot', 'R_pArc_vot');
    % % and now load them, just to be shure they work fine.
    % load( fullfile(outdir, 'VOF_all.mat'));
    % 
    
    
    
    %% RENDER
    % % Read the ROIs in the cortex
%     roi_L_fusiform = dtiReadRoi(fullfile(dmridir,'ROIs', ...
%                                 '1007_ctx-lh-fusiform.mat'));
%     roi_L_inferiortemporal = dtiReadRoi(fullfile(dmridir,'ROIs', ...
%                                              '1009_ctx-lh-inferiortemporal.mat'));
%     roi_L_lateraloccipital = dtiReadRoi(fullfile(dmridir,'ROIs', ...
%                                              '1011_ctx-lh-lateraloccipital.mat'));                                     
%     
%     
%     % Rnder them the fibers and the cortex rois
%     AFQ_RenderFibers(L_VOF , 'color',  [158 47 88]/256, ...
%                     'tubes',[0]); % Render the fibers
%     AFQ_RenderFibers(L_pArc , 'color',  [237 139 140]/256, 'tubes',[0], ...
%                      'newfig',false); % Render the fibers
%     % Render the roi in orange
%     AFQ_RenderRoi(roi_L_fusiform, [241 217 201]/256, 'mesh'); 
%     AFQ_RenderRoi(roi_L_inferiortemporal, [241 217 181]/256, 'mesh'); 
%     AFQ_RenderRoi(roi_L_lateraloccipital, [241 217 161]/256, 'mesh'); 
% 
%     % con este comando binarizo aparc+aseg y ademas solo me quedo con el GM
    FSLDIR = '/opt/fsl/fsl-5.0.9/fsl';
    setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
    setenv('FSLDIR', FSLDIR);
    system([FSLDIR fsp 'bin' fsp 'fslmaths ' dmridir fsp 'aparc+aseg.nii.gz '...
            '-thr 1000 -bin ' dmridir fsp 'segmentation.nii.gz']);
%     % Y ahora creo el mesh
%}     
    
     im = [dmridir fsp 'segmentation.nii.gz'];
   % msh = AFQ_meshCreate(im, 'color', [.8 .7 .6])
% 
% 
% 

    
    %% Obtener el file con el intersect entre fiber y cortex
    % por ahora solo hago el posterior arcuate y el VOF
    % fiberRois = {L_VOF, L_pArc, L_pArc_vot};
    fiberRois = {L_arcuate};
    for fr =1:length(fiberRois)
        fg = fiberRois{fr};
        segmentation = niftiRead(im);
        fdImg = zeros([size(segmentation.data) length(fg)]);
        % Extraido de AFQ_RenderFibersOnCortex
        % Check if the segmentation is binary or is mrVista format
        if length(unique(segmentation.data(:)))>2
            segmentation.data = uint8(segmentation.data==3 | segmentation.data==4);
        end
        for ii = 1:length(fg)
            fdImg(:,:,:,ii) = smooth3(dtiComputeFiberDensityNoGUI(...
                                      fg(ii), ... % Fibras en .mat
                                      segmentation.qto_xyz, ... % matriz del segmentation.nii.gz
                                      size(segmentation.data), ... % tamano en voxels
                                      1, ... % = 1, Normalize to 1. =0, fiber count 
                                      [],... % FibreGroupNum: si quieres elegir solo alguna fibra concreta
                                      0), ...% endptFlag=1, solo usar fiber endpoints. LO CAMBIO!!
                               'gaussian', ...
                               5); 
        end

        % Tack on an extra volume that will mark voxels with no fibers
        fdImg = cat(4,zeros(size(fdImg(:,:,:,1)))+.000001,fdImg);
        % Find the volume with the highest fiber density in each voxel
        [~,fdMax] = max(fdImg,[],4);
        % clear fdImg; % Lo inicializo arriba con zeros a ver si arregla el
        % parfor

        % Zero out voxels with no fibers
        fdMax = fdMax-1;
        % Make into a nifti volume
        fdNii = segmentation;
        fdNii.data = fdMax;

        % niftiWrite(fdNii, fdNii.fname)

        % Render it
        % [p, msh, lightH] =  AFQ_RenderCorticalSurface(segmentation, ...
        %                         'overlay',fdNii, ...
        %                         'boxfilter',1, ...
        %                         'thresh',[1 20], ...
        %                         'interp','nearest', ...
        %                         'cmap',colormap);


        % Al archivo anterior le he dicho que escriba el nifti con el overlap entre
        % los tractos y la corteza, que esta dada por el archivo de aparc+aseg. 
        % Prueba 1. visualizarlo a ver que tal se ve y ver si podre hacer overlay al
        % espacio individual.
        % Prueba 2. Podria crear ya los ROIs metidos un par de mm hacia dentro, o
        % sea, estaran en volumen, y luego podria salvar los tractos en nifti tb y
        % ver el overlap, luego inflar y buscar el overlap con el white matter...

        % niftiRead-Write y MRIread-write hacen cosas diferentes e inservibles en
        % freeview, aunque en mrview de mrtrix se vieran bien.

        % fdNii.fname = [dmridir fsp fg.name '_overlayGM_vista.nii.gz'];
        % niftiWrite(fdNii, fdNii.fname)

        % No hace falta escribirlo en formato mrVista, ya que estos no se
        % ven vien en freeview, solo se ven bien en mrview de mrtrix, pero
        % los de MRIwrite si he conseguido que se vean igual tanto en uno
        % como en otro.

        % en fs ahora
        segRead = MRIread([dmridir fsp 'segmentation.nii.gz']);
        segRead.vol = permute(fdNii.data, [2 1 3]);  % mierdas de x,y en Matlab
        MRIwrite(segRead, [dmridir fsp strrep(fg.name,' ','_') '_tracts.nii.gz']);
        % MRIwrite(segRead, [dmridir fsp 'NORM_' fg.name '_tracts.nii.gz']);


        % Hay que pensar si hago el paso a la superficie con todos los
        % voxeles que pertenecen a los tractos, o solo me quedo con
        % aquellos voxeles que coinciden con los ROI de interes y luego
        % hago el paso a la superifice. >> He pasado todos los voxeles,
        % luego con aparc podre elegir los voxeles que me interesen para
        % los rois. Lo de los tractos tiene que ser bidireccional. 

        % Y ahora los convertimos a superficie usando fs
        movname    = fullfile(dmridir, [strrep(fg.name,' ','_') '_tracts.nii.gz']);
        oname      = fullfile(dmridir, [strrep(fg.name,' ','_') '_tracts.mgh']);
        oname305   = fullfile(dmridir, [strrep(fg.name,' ','_') '_tracts305.mgh']);
        % movname    = fullfile(dmridir, ['NORM_' fg.name '_tracts.nii.gz']);
        % oname      = fullfile(dmridir, ['NORM_' fg.name '_tracts.mgh']);
        % oname305   = fullfile(dmridir, ['NORM_' fg.name '_tracts305.mgh']);

        % fshomecajal02 = '/usr/local/freesurfer';
        % fsbincajal02 = '/usr/local/freesurfer/bin';

        % setenv('FREESURFER_HOME', fshome);       
        % Uso --projfrac -1 para meterlo un poco dentro del cortex, si no
        % se ve mucho mas cuarteado. He probado con -2 y -3 pero casi no
        % hay mejora. Al final el problema es que a los gyrus no llegan las
        % fibras.
        cmd2 =  [fsbin fsp 'mri_vol2surf ' ...
                   '--srcsubject '  subname  ' ' ...
                   '--projdist -1 ' ... % '--projfrac 0.5 ' ... %  
                   '--interp trilinear ' ...
                   '--hemi lh ' ...
                   '--regheader '  subname  ' ' ...
                   '--mov '  movname  ' ' ...
                   '--o '  oname ...
                   ];
        cmd3 = [fsbin fsp 'mri_surf2surf ' ...
                   '--srcsubject '  subname  ' ' ...
                   '--srchemi lh ' ...
                   '--srcsurfreg sphere.reg ' ...
                   '--sval '  oname   ' ' ...
                   '--trgsubject fsaverage ' ...
                   '--trghemi lh ' ...
                   '--trgsurfreg sphere.reg ' ...
                   '--tval '  oname305  ' ' ...
                   '--sfmt ' ...
                   '--curv ' ...
                   '--noreshape ' ...
                   '--no-cortex ' ...
                   ];

        system(cmd2);
        system(cmd3);
% 
% 
% %         cortex = fullfile(dmridir, 'segmentation.nii.gz');
% %         % overlay = fullfile(AFQdata,'mesh','Left_Arcuate_Endpoints.nii.gz');
% %         thresh = .01; % Threshold for the overlay image
% %         crange = [.01 .8]; % Color range of the overlay image
% %         % Render the cortical surface colored by the arcuate endpoint density 
% %         [p, msh, lightH] = AFQ_RenderCorticalSurface(cortex, 'overlay' , overlay, 'crange', crange, 'thresh', thresh)
% % 
% %         msh = AFQ_meshCreate(cortex, 'color', [.8 .7 .6])
% %         AFQ_RenderCorticalSurface(msh)
% % 
% 
% 
% 
% % 
% %         %% Ahora voy a ir con la siguiente solucion en mrtrix para freesurfer
% %         % % If you use the read_mrtrix_tracks.m matlab function you can load
% %         % .tck files into a matlab structure. Then run a simple loop to keep 
% %         % the first and last coordinates of each streamline in the .data structure.
% %         % % The streamline coordinates should be in mm space which you can then 
% %         % match to freesurfer vertices as follows...
% %         % % Load a freesurfer surface (e.g. lh.white) into matlab using the 
% %         % read_surf.m function provided in the set of freesurfer matlab functions. 
% %         % The vertex_coords variable gives mm coordinates of each vertex. 
% %         % You can then find the Euclidean distance between an end point and the 
% %         % vertices to find the nearest vertex for a fiber termination.
% %         % % Freesurfer then has a bunch of matlab functions to write surface 
% %         % overlays or annotation files depending on your desired outcome
% %         % (e.g. save_mgh).
% %         fname = 'WordHighVsPhaseScrambledWords_Sphere4.tck';
% %         fname = 'WordHighVsFF_Sphere5.tck';
% % 
% % 
% %         data =  read_mrtrix_tracks(fullfile(dmridir, 'dti90trilin','mrtrix',fname));
% %         endPoints = zeros(2*length(data.data), 3);
% %         for ii =1:(2*length(data.data))
% %             tractNo = ceil(ii/2);
% %             if mod(ii,2)
% %                 endPoints(ii,:) = data.data{tractNo}(1,:);
% %             else
% %                 endPoints(ii,:) = data.data{tractNo}(end,:);
% %             end
% %         end
% %         WhiteSurf = read_surf(fullfile(fs_SUBJECTS_DIR,subname,'surf','lh.white'));
% % 
% %         % Find the index and coordinate of closest vertex
% %         vertexIndex = knnsearch(WhiteSurf, endPoints);
% %         vertexPoints = WhiteSurf(knnsearch(WhiteSurf, endPoints),:);
% % 
% %         % Write it
% %         ok = write_label(vertexIndex,[], [], ...
% %                      fullfile(dmridir, 'dti90trilin','mrtrix',[fname '.label']));
% 
% 
    end

  end
end

%% Convertir VOF y PARC a tcks para usar en mrtrix (se hizo a posteriori)
if(0)
   for ns = 1 : length(subs) % no funcional parfor pero es rapido
        subname = subs(ns).name

        setenv('FREESURFER_HOME', fshome); 

        % Read the VOF per every subject
        dmridir = fullfile(DWIdir, subname, 'dmri');
        wholebrainfgPath= fullfile(dmridir, 'dti90trilin', 'fibers'); 
        MRtrixPath= fullfile(dmridir, 'dti90trilin', 'mrtrix'); 
        outdir = fullfile(wholebrainfgPath,'VOF');
        cd(dmridir)
        load( fullfile(outdir, 'VOF_all.mat'));

        % Create empty struct to put the tract data
        tractData = struct(...
                        'act', ['/bcbl/home/public/Gari/MINI/ANALYSIS/DWI/' subname '/dmri/dti90trilin/mrtrix/data_aligned_trilin_noMEC_5tt.mif'], ...
                  'backtrack', '0', ...
          'downsample_factor', '3', ...
                  'fod_power', '0.25', ...
             'init_threshold', '0.100000001', ...
                       'lmax', '8', ...
                  'max_angle', '45', ...
           'max_num_attempts', '50000000', ...
             'max_num_tracks', '500000', ...
          'max_seed_attempts', '1', ...
                 'max_trials', '1000', ...
                     'method', 'iFOD2', ...
             'mrtrix_version', '0.3.15-65-gaeb862d2', ...
           'output_step_size', '1', ...
                        'rk4', '0', ...
           'samples_per_step', '4', ...
             'sh_precomputed', '1', ...
                     'source', ['/bcbl/home/public/Gari/MINI/ANALYSIS/DWI/' subname '/dmri/dti90trilin/mrtrix/data_aligned_trilin_noMEC_wmCsd_lmax4.mif'], ...
                  'step_size', '1', ...
        'stop_on_all_include', '0', ...
                  'threshold', '0.100000001', ...
                  'timestamp', '1488215560.9647302628', ...
             'unidirectional', '0', ...
                   'datatype', 'Float32LE', ...
                      'count', num2str(size(cellfun(@transpose,L_VOF.fibers,'un',0), 1)), ...
                'total_count', '500000', ...
                       'data', {cellfun(@transpose,L_VOF.fibers,'un',0)'} );
          tractData.data = cellfun(@transpose,L_VOF.fibers,'un',0)';
          write_mrtrix_tracks(tractData, [MRtrixPath fsp 'afq_L_vOF.tck']);
          % Now write the pAF
          tractData.count = num2str(size(cellfun(@transpose,L_pArc.fibers,'un',0), 1));
          tractData.data  = cellfun(@transpose,L_pArc.fibers,'un',0)';
          write_mrtrix_tracks(tractData, [MRtrixPath fsp 'afq_L_pAF.tck']);
    end      
end

%% Convertir los avg sem y perc ROIs de VOT IFG PPC a individual space
if(0)
    
    parfor ns = 1 : length(subs)
        ROIs = {'PPC_perc_averages', 'PPC_sem_averages', ...
                'VOT_perc_averages', 'VOT_sem_averages', ...
                'IFG_perc_averages', 'IFG_sem_averages'};
        dilateLabelBy = '1';
        subname = subs(ns).name

          setenv('FREESURFER_HOME', fshome); 
          setenv('SUBJECTS_DIR', fs_SUBJECTS_DIR);
        
        
        for ROI = ROIs
            roi = ROI{:};

            iname = fullfile(fs_SUBJECTS_DIR, 'fsaverage', 'label', ...
                                                ['lh.' roi dilateLabelBy '.label']);
            oname = fullfile(fs_SUBJECTS_DIR, subname, 'label', ...
                                                ['lh.' roi dilateLabelBy '.label']);


           cmd = [fsbin fsp 'mri_label2label ' ...
                   '--srcsubject fsaverage ' ...
                   '--hemi lh ' ...
                   '--srclabel '  iname   ' ' ...
                   '--trgsubject '  subname  ' ' ...
                   '--trglabel '  oname  ' ' ...
                   '--regmethod surface '];
            system(cmd)
        end   
    end   
    
end
      
%% Crear los seis pares de tractos creando esferas en esos puntos
% No he conseguido hacer funcionar mrtrix dentro de parfor...
if(0)    

    for ns = 1 : length(subs)
        ROIs = {'PPC_perc_averages', 'PPC_sem_averages', ...
                 'VOT_perc_averages', 'VOT_sem_averages', ...
                 'IFG_perc_averages', 'IFG_sem_averages'};

%         tcktype       = {'sem', 'perc'};
%         connections1   = {'VOT2IFG', 'VOT2PPC', 'PPC2IFG'}
%         connections2   = {'VOT2VOT', 'PPC2PPC'};
%         dilateLabelBy = '1';
%         setenv('FREESURFER_HOME', fshome); 
%         setenv('SUBJECTS_DIR', fs_SUBJECTS_DIR);
%         oldPath = getenv('PATH');

        ROIs = {'VOT_sem_averages', 'IFG_sem_averages'};
        tcktype       = {'sem'};  %, 'perc'};
        connections   = {'VOT2IFG'}; %, 'VOT2PPC', 'PPC2IFG'};
        dilateLabelBy = '1';
        setenv('FREESURFER_HOME', fshome); 
        setenv('SUBJECTS_DIR', fs_SUBJECTS_DIR);
%         oldPath = getenv('PATH');
%         setenv('PATH', ...
%                 ['/bcbl/home/home_g-m/glerma/GIT/mrtrix3/release/bin:' ...
%                 '/bcbl/home/home_g-m/glerma/GIT/mrtrix3/scripts:' ...
%                 oldPath]);

%         ROISphereRadius = [8,10,12];
        ROISphereRadius = [12];
        % WholeTractogramName = 'data_aligned_trilin_noMEC_wmCsd_lmax4_data_aligned_trilin_noMEC_wmMask_data_aligned_trilin_noMEC_wmMask_iFOD2-500000.tck';        
        WholeTractogramName = 'data_alignedNorm_trilin_noMEC.tck';
        subname = subs(ns).name
        dmridir = fullfile(DWIdir, subname, 'dmri');
        mrtrixdir = fullfile(dmridir, 'dti90trilin','mrtrix');
        % T1std = MRIread([dmridir fsp  't1_std_acpc.nii.gz']);
        T1 = MRIread([fs_SUBJECTS_DIR fsp  subname fsp 'mri' fsp 'T1.mgz']);
        %% Convertimos los one vertex voxels al espacio individual
        coords = struct();
        for ROI = ROIs
            roi = ROI{:};
            label = read_label(subname, ['lh.' roi dilateLabelBy]);
            surfRAS =  label(1, 2:4);
            % Convertir desde surfaceRAS a scannerRas
            scanRAS  =  T1.vox2ras  * inv(T1.tkrvox2ras) *  [surfRAS';1];
            coords.([roi dilateLabelBy]) = scanRAS;
        end   
        %% Creamos los tractos con las esferas
        % Primero hicimos el triangulo sem y perc por separado
        for tckt = tcktype;for conn = connections1;for sphR=ROISphereRadius
            tctname = ['L_' tckt{:} '_' conn{:} '_R' num2str(sphR) '.tck'];
            con1 = conn{:}(1:3);
            coord1 = round(coords.([con1 '_' tckt{:} '_averages' dilateLabelBy])');
            coord1(4) = sphR;
            roi1 = strjoin(arrayfun(@(x) num2str(x),coord1,'UniformOutput',false),',');
            
            con2 = conn{:}(5:7);
            coord2 = round(coords.([con2 '_' tckt{:} '_averages' dilateLabelBy])');
            coord2(4) = sphR;
            roi2 = strjoin(arrayfun(@(x) num2str(x),coord2,'UniformOutput',false),',');
            
%             cmd_str = ['tckedit ' ...
%                        '-include ' roi1  ' ' ...
%                        '-include ' roi2  ' ' ... %'-ends_only ' ... 
%                          [mrtrixdir filesep WholeTractogramName] ' ' ...
%                          [mrtrixdir filesep 'sinEndsOnly_' tctname]];
            cmd_str = ['tckgen -algorithm iFOD2  ' ...
                       '-seed_sphere ' roi1  ' ' ...
                       '-include ' roi2  ' ' ...
                       '-act ' mrtrixdir filesep 'data_aligned_trilin_noMEC_5tt.mif ' ...
                       '-number 5000 ' ...
                       mrtrixdir filesep 'data_aligned_trilin_noMEC_wmCsd_lmax4.mif ' ...
                       mrtrixdir filesep 'vMC_ILF_UNC_v01.tck'];
            bkgrnd = false;
            verbose = false;
            mrtrixVersion = 3;
            AFQ_mrtrix_cmd(cmd_str, bkgrnd, verbose,mrtrixVersion)
            
        end;end;end
        % Luego VOT2VOT y PPC2PPC
        for conn = connections2;for sphR=ROISphereRadius
            tctname = ['L_percsem_' conn{:} '_R' num2str(sphR) '.tck'];
            con1 = conn{:}(1:3);
            coord1 = round(coords.([con1 '_sem_averages' dilateLabelBy])');
            coord1(4) = sphR;
            roi1 = strjoin(arrayfun(@(x) num2str(x),coord1,'UniformOutput',false),',');
            
            con2 = conn{:}(5:7);
            coord2 = round(coords.([con2 '_perc_averages' dilateLabelBy])');
            coord2(4) = sphR;
            roi2 = strjoin(arrayfun(@(x) num2str(x),coord2,'UniformOutput',false),',');
            
            cmd_str = ['tckedit ' ...
                       '-include ' roi1  ' ' ...
                       '-include ' roi2  ' ' ...
                       '-ends_only ' ... 
                         [mrtrixdir filesep WholeTractogramName] ' ' ...
                         [mrtrixdir filesep tctname]];
            bkgrnd = false;
            verbose = false;
            mrtrixVersion = 3;
            AFQ_mrtrix_cmd(cmd_str, bkgrnd, verbose,mrtrixVersion)
        end;end
    
    end
    
end

%% Obtener las estadisticas sin SIFT

% Lo iba a hacer fuera, pero tengo todas las variables aqui ya...
if(0)    
    resultsFldr = fullfile(DWIdir, 'ResultsCSV');
    for ns = 1 : length(subs)
        ROIs = {'PPC_perc_averages', 'PPC_sem_averages', ...
                 'VOT_perc_averages', 'VOT_sem_averages', ...
                 'IFG_perc_averages', 'IFG_sem_averages'};
        tcktype       = {'sem', 'perc'};
        connections   = {'VOT2IFG', 'VOT2PPC', 'PPC2IFG'};
        dilateLabelBy = '1';
        setenv('FREESURFER_HOME', fshome); 
        setenv('SUBJECTS_DIR', fs_SUBJECTS_DIR);
        oldPath = getenv('PATH');
        setenv('PATH', ...
                ['/bcbl/home/home_g-m/glerma/GIT/mrtrix3/release/bin:' ...
                '/bcbl/home/home_g-m/glerma/GIT/mrtrix3/scripts:' ...
                oldPath]);

        subname = subs(ns).name
        dmridir = fullfile(DWIdir, subname, 'dmri');
        mrtrixdir = fullfile(dmridir, 'dti90trilin','mrtrix');
        faFile = fullfile(mrtrixdir, 'data_aligned_trilin_noMEC_fa.mif');
        MTVfile = fullfile(qmridir,subname,'OutPutFiles_1','BrainMaps','TV_map.nii.gz');
        T1qfile = fullfile(qmridir,subname,'OutPutFiles_1','BrainMaps','T1_map_Wlin.nii.gz');

        
           
        %% Por cada tracto que hemos leido antes, escribimos stats
        ROISphereRadius = [12]; % 8 y 12 son muy pequeos, igual hasta meter 15
        for tckt = tcktype;for conn = connections;for sphR=ROISphereRadius
            tctname = ['L_' tckt{:} '_' conn{:} '_R' num2str(sphR) '.tck'];
            
            tckFile = [mrtrixdir filesep tctname];
            csvFA  = fullfile(resultsFldr,['FA_' tctname '_' subname '.csv']);
            csvMTV = fullfile(resultsFldr,['MTV_' tctname '_' subname '.csv']);
            csvT1q = fullfile(resultsFldr,['T1q_' tctname '_' subname '.csv']);
            
            cmdFA  = ['tcksample -stat_tck mean ' tckFile ' ' faFile ' ' csvFA];
            cmdMTV = ['tcksample -stat_tck mean ' tckFile ' ' MTVfile ' ' csvMTV];
            cmdT1q = ['tcksample -stat_tck mean ' tckFile ' ' T1qfile ' ' csvT1q];

            bkgrnd = false;
            verbose = false;
            mrtrixVersion = 3;
            AFQ_mrtrix_cmd(cmdFA,  bkgrnd, verbose, mrtrixVersion)
            AFQ_mrtrix_cmd(cmdMTV, bkgrnd, verbose, mrtrixVersion)
            AFQ_mrtrix_cmd(cmdT1q, bkgrnd, verbose, mrtrixVersion)
            
        end;end;end
        tctname = ['afq_L_vOF.tck'];
            tckFile = [mrtrixdir filesep tctname];
            csvFA  = fullfile(resultsFldr,['FA_' tctname '_' subname '.csv']);
            csvMTV = fullfile(resultsFldr,['MTV_' tctname '_' subname '.csv']);
            csvT1q = fullfile(resultsFldr,['T1q_' tctname '_' subname '.csv']);
            
            cmdFA  = ['tcksample -stat_tck mean ' tckFile ' ' faFile ' ' csvFA];
            cmdMTV = ['tcksample -stat_tck mean ' tckFile ' ' MTVfile ' ' csvMTV];
            cmdT1q = ['tcksample -stat_tck mean ' tckFile ' ' T1qfile ' ' csvT1q];

            bkgrnd = false;
            verbose = false;
            mrtrixVersion = 3;
            AFQ_mrtrix_cmd(cmdFA,  bkgrnd, verbose, mrtrixVersion)
            AFQ_mrtrix_cmd(cmdMTV, bkgrnd, verbose, mrtrixVersion)
            AFQ_mrtrix_cmd(cmdT1q, bkgrnd, verbose, mrtrixVersion)
        
        tctname = 'afq_L_pAF.tck';
            tckFile = [mrtrixdir filesep tctname];
            csvFA  = fullfile(resultsFldr,['FA_' tctname '_' subname '.csv']);
            csvMTV = fullfile(resultsFldr,['MTV_' tctname '_' subname '.csv']);
            csvT1q = fullfile(resultsFldr,['T1q_' tctname '_' subname '.csv']);
            
            cmdFA  = ['tcksample -stat_tck mean ' tckFile ' ' faFile ' ' csvFA];
            cmdMTV = ['tcksample -stat_tck mean ' tckFile ' ' MTVfile ' ' csvMTV];
            cmdT1q = ['tcksample -stat_tck mean ' tckFile ' ' T1qfile ' ' csvT1q];

            bkgrnd = false;
            verbose = false;
            mrtrixVersion = 3;
            AFQ_mrtrix_cmd(cmdFA,  bkgrnd, verbose, mrtrixVersion)
            AFQ_mrtrix_cmd(cmdMTV, bkgrnd, verbose, mrtrixVersion)
            AFQ_mrtrix_cmd(cmdT1q, bkgrnd, verbose, mrtrixVersion)        
        
    end
    
end

%% doBias: En vez de iPython, lanzando por aqui los comandos de batch_fslpreproc
% Lo iba a hacer fuera, pero tengo todas las variables aqui ya...
if(1)
    for ns = 1 : length(subs)
        subname = subs(ns).name

        % Folders
        dmridir = fullfile(DWIdir, subname, 'dmri');
        mrtrixdir = fullfile(dmridir, 'dti90trilin','mrtrix');
        faFile = fullfile(mrtrixdir, 'data_aligned_trilin_noMEC_fa.mif');

        % Options
        doPreProc = 0
        doBias = 1
        doDtiInit = 0
        doAfqCreate = 0
        doAfqRun = 0
        batch_fslpreprocessdiffusion(subname, AnalysisDir, ...
                                     doPreProc, doBias,...
                                     doDtiInit, doAfqCreate, doAfqRun)
    end
end

%% Tengo que crear los simbolic links a los DWI
% Lo iba a hacer fuera, pero tengo todas las variables aqui ya...
if(0)
    for ns = 1 : length(subs)
        subname = subs(ns).name

        % Folders
        dmridir = fullfile(DWIdir, subname, 'dmri');
        mrtrixdir = fullfile(dmridir, 'dti90trilin','mrtrix');
        faFile = fullfile(mrtrixdir, 'data_aligned_trilin_noMEC_fa.mif');

        
        
        % No   funciona pq antes tengo que meter bvacs y    bvals
%         
%         cmd0 = ['mrconvert -fslgrad ' dmridir fsp 'eddy' fsp 'bvecs ' ...
%                  dmridir fsp 'eddy' fsp 'bvals ' ...
%                  dmridir fsp 'eddy' fsp 'biasdata.nii.gz ' ...
%                  dmridir fsp 'eddy' fsp 'biasdata.mif'];
%             bkgrnd = false;
%             verbose = false;
%             mrtrixVersion = 3;
%         AFQ_mrtrix_cmd(cmd0, bkgrnd, verbose,mrtrixVersion)
        DWI2dir = '/export/home/glerma/glerma/00local/PROYECTOS/MINI/ANALYSIS/DWI';
        cmd1 = ['ln -s ' dmridir fsp 'eddy' fsp 'biasdata.mif ' ...
                         DWI2dir fsp 'input_dwi_folder' fsp subname '.mif'];
        cmd2 = ['ln -s ' dmridir fsp 'eddy' fsp 'nodif_brain_mask.nii.gz ' ...
                         DWI2dir fsp 'input_brain_mask_folder' fsp subname '.nii.gz'];
        system(cmd1)
        system(cmd2)

    end
end

%% Los archivos normalizados son .mif, convertirlos a .nii.gz para que lea mrDiffussion y seguir con el pipeline de antes
% Lo iba a hacer fuera, pero tengo todas las variables aqui ya...
if(0)
    for ns = 1 : length(subs)
        %setenv('PATH','/bcbl/home/home_g-m/glerma/GIT/mrtrix3/release/bin:/bcbl/home/home_g-m/glerma/GIT/mrtrix3/scripts');
        subname = subs(ns).name
        DWI2dir = '/export/home/glerma/glerma/00local/PROYECTOS/MINI/ANALYSIS/DWI';
        output_normalised_dwi_folder = [DWI2dir fsp 'output_normalised_dwi_folder'];
        cmd_str = ['mrconvert ' output_normalised_dwi_folder fsp subname '.mif ' ...
               output_normalised_dwi_folder fsp subname '.nii.gz'];
        bkgrnd = false;
        verbose = false;
        mrtrixVersion = 3;
        AFQ_mrtrix_cmd(cmd_str, bkgrnd, verbose,mrtrixVersion)
    end
end

%% No hay sitio en public, crear carpeta en home folder y link desde public
% Lo iba a hacer fuera, pero tengo todas las variables aqui ya...
if(0)
      oldDWI = '/bcbl/home/public/Gari/MINI/ANALYSIS/DWI';
      newDWI = '/bcbl/home/home_g-m/glerma/00local/PROYECTOS/MINI/ANALYSIS/DWI';
      
      for ns = 4 : length(subs)
        subname = subs(ns).name

        %Rename old dti90trilin directories
        olddmri = fullfile(oldDWI, subname, 'dmri','dti90trilin');
        noNorm_olddmri = fullfile(oldDWI, subname, 'dmri','noNorm_dti90trilin');
%         system(['mv ' olddmri ' '  noNorm_olddmri])
%         
%         % Rename old data.nii.gz too noNorm
%           oldData = fullfile(oldDWI, subname, 'dmri','eddy','data.nii.gz');
%           noNormOldData = fullfile(oldDWI, subname, 'dmri','eddy','noNorm_data.nii.gz');
%           system(['mv ' oldData ' '  noNormOldData])
%           % New Folders
%         subjectdir = fullfile(newDWI, subname);
%         dmridir = fullfile(subjectdir, 'dmri');
%         trilindir = fullfile(dmridir, 'dti90trilin');
%         
%         % Create new directories
%         mkdir(subjectdir)
%         mkdir(dmridir)
%         mkdir(trilindir)
%         
%         % Create symbolic links in public
%         cmd1 = ['ln -s ' trilindir ' ' olddmri];
%         system(cmd1)
% 
%         % Create symbolic links in public
%         output_normalised_dwi_folder = [newDWI fsp 'output_normalised_dwi_folder'];
%         normData = fullfile(output_normalised_dwi_folder,[subname '.nii.gz']); % Path to the data
%         cmd1 = ['ln -s ' normData ' ' oldData];
%         system(cmd1)
        
%         old5tt = 'data_aligned_trilin_noMEC_5tt.mif';
%         nuevo5tt = 'data_alignedNorm_trilin_noMEC_5tt.mif';
%         zahar  = fullfile(noNorm_olddmri,'mrtrix' , old5tt  );
%         berria = fullfile(olddmri, 'mrtrix', nuevo5tt);
%         mkdir(fullfile(olddmri, 'mrtrix'))
%         cmd  = ['ln -s ' zahar ' ' berria]
%         system(cmd)
        
        % Cambiar el nombre a los afqOut.mat que ya existen, para el S001 y
        % el S002 no los he podido salvar, pero los otros guardar por si
        % acaso
%         zahar = fullfile(oldDWI,subname,'afqOut.mat');
%         berria = fullfile(oldDWI,subname,'afqOut_noNorm.mat');
%         cmd  = ['mv ' zahar ' ' berria]
%         system(cmd)
    end
end

%% doDtiInit: En vez de iPython, lanzando por aqui los comandos de batch_fslpreproc
% Lo iba a hacer fuera, pero tengo todas las variables aqui ya...
if(0)
    clear all; close all; 
    fsp = filesep;
    subs = dir('S*');
    for ns =2 1 : length(subs)
        subname = subs(ns).name
        subName = subname;
        % Folders
        AnalysisDir = '/bcbl/home/public/Gari/MINI/ANALYSIS';
        
        % Options
        doPreProc = 0;
        doBias = 0;
        doDtiInit = 1;
        doAfqCreate = 1;
        doAfqRun = 1;
        tic
        batch_fslpreprocessdiffusion(subname, AnalysisDir, ...
                                     doPreProc, doBias,...
                                     doDtiInit, doAfqCreate, doAfqRun);
        toc                         
    end
end
% batch_fslpreprocessdiffusion('S002', '/bcbl/home/public/Gari/MINI/ANALYSIS', 0, 0,1, 1, 1);
% batch_fslpreprocessdiffusion('S002',  '/bcbl/home/public/Gari/MINI/ANALYSIS', 0, 0,0, 0, 1);
     
%% dwi2fod, tckgen, SIFT2 y Crear los seis pares de tractos creando esferas en esos puntos
% No he conseguido hacer funcionar mrtrix dentro de parfor...
if(0)    

pDWI = '/bcbl/home/public/Gari/MINI/ANALYSIS/DWI/';
hDWI = '/bcbl/home/home_g-m/glerma/00local/PROYECTOS/MINI/ANALYSIS/DWI';
fsp = filesep;
    % for ns = 1 : length(subs)
    for ns = 1  : 12
    % for ns = 13  : 24
    % for ns = 25 : 36 
    % for ns = 37 : 48
    % for ns = 49  : 60
    % for ns = 61 : 72 
    % for ns = 73 : 84
    % for ns = 85 : 97

        tic
        numFiberNum = '5000000';
        numFiberName = '5M';
        
        subname = subs(ns).name
        hacerCSD = 0
        hacerTCK = 1
        hacerSIFT2 = 1
        
        if hacerCSD 
            tic
            system(['dwi2fod msmt_csd -force -nthreads 16 ' ...
            '-grad ' pDWI fsp  subname '/dmri/dti90trilin/mrtrix/data_alignedNorm_trilin_noMEC.b ' ...
            '-mask  /bcbl/home/public/Gari/MINI/ANALYSIS/DWI/' subname '/dmri/dti90trilin/mrtrix/data_alignedNorm_trilin_noMEC_brainmask.mif ' ...
            pDWI fsp subname '/dmri/dti90trilin/mrtrix/data_alignedNorm_trilin_noMEC_dwi.mif ' ...
            hDWI fsp 'output_group_average_response/wm_output_group_average_response.txt ' ...
            pDWI fsp subname '/dmri/dti90trilin/mrtrix/data_alignedNorm_trilin_noMEC_wmCsd_lmax4.mif ' ...
            hDWI fsp 'output_group_average_response/gm_output_group_average_response.txt ' ...
            pDWI fsp subname '/dmri/dti90trilin/mrtrix/data_alignedNorm_trilin_noMEC_gmCsd_lmax4.mif ' ...
            hDWI fsp 'output_group_average_response/cs_output_group_average_response.txt ' ...'
            pDWI fsp subname '/dmri/dti90trilin/mrtrix/data_alignedNorm_trilin_noMEC_csfCsd_lmax4.mif'])
            toc
       end
               
        if hacerTCK
%             tic
            system(['tckgen /bcbl/home/public/Gari/MINI/ANALYSIS/DWI/' subname '/dmri/dti90trilin/mrtrix/data_alignedNorm_trilin_noMEC_wmCsd_lmax4.mif '  ...
                       '-algo iFOD2 ' ...
                       '-seed_image /bcbl/home/public/Gari/MINI/ANALYSIS/DWI/' subname '/dmri/dti90trilin/mrtrix/data_alignedNorm_trilin_noMEC_wmCsd_lmax4.mif ' ...
                       '-act /bcbl/home/public/Gari/MINI/ANALYSIS/DWI/' subname '/dmri/dti90trilin/mrtrix/data_alignedNorm_trilin_noMEC_5tt.mif ' ...
                       '-num ' numFiberNum ' ' ...
                       '/bcbl/home/public/Gari/MINI/ANALYSIS/DWI/' subname '/dmri/dti90trilin/mrtrix/' numFiberName 'data_alignedNorm_trilin_noMEC.tck']);
%             toc
        end
        
        if hacerSIFT2
%           tic
          system(['tcksift2 -force -nthreads 16 ' ...
                    '-act ' pDWI fsp subname '/dmri/dti90trilin/mrtrix/data_alignedNorm_trilin_noMEC_5tt.mif ' ...
                    pDWI fsp subname '/dmri/dti90trilin/mrtrix/' numFiberName 'data_alignedNorm_trilin_noMEC.tck ' ...                    
                    pDWI fsp subname '/dmri/dti90trilin/mrtrix/data_alignedNorm_trilin_noMEC_wmCsd_lmax4.mif ' ...
                    pDWI fsp subname '/dmri/dti90trilin/mrtrix/' numFiberName 'data_alignedNorm_trilin_noMEC.sift2.weight']);
%           toc
        end
        
        ROIs = {'PPC_perc_averages', 'PPC_sem_averages', ...
                 'VOT_perc_averages', 'VOT_sem_averages', ...
                 'IFG_perc_averages', 'IFG_sem_averages'};
        tcktype       = {'sem', 'perc'};
        connections1   = {'VOT2IFG', 'VOT2PPC', 'PPC2IFG'}
        connections2   = {'VOT2VOT', 'PPC2PPC'};
        dilateLabelBy = '1';
        setenv('FREESURFER_HOME', fshome); 
        setenv('SUBJECTS_DIR', fs_SUBJECTS_DIR);
        oldPath = getenv('PATH');

%         ROISphereRadius = [8,10,12];
        ROISphereRadius = [12];
        % WholeTractogramName = 'data_aligned_trilin_noMEC_wmCsd_lmax4_data_aligned_trilin_noMEC_wmMask_data_aligned_trilin_noMEC_wmMask_iFOD2-500000.tck';        
        WholeTractogramName =  [numFiberName 'data_alignedNorm_trilin_noMEC.tck'];
        dmridir = fullfile(DWIdir, subname, 'dmri');
        mrtrixdir = fullfile(dmridir, 'dti90trilin','mrtrix');
        % T1std = MRIread([dmridir fsp  't1_std_acpc.nii.gz']);
        T1 = MRIread([fs_SUBJECTS_DIR fsp  subname fsp 'mri' fsp 'T1.mgz']);
        %% Convertimos los one vertex voxels al espacio individual
        coords = struct();
        for ROI = ROIs
            roi = ROI{:};
            label = read_label(subname, ['lh.' roi dilateLabelBy]);
            surfRAS =  label(1, 2:4);
            % Convertir desde surfaceRAS a scannerRas
            scanRAS  =  T1.vox2ras  * inv(T1.tkrvox2ras) *  [surfRAS';1];
            coords.([roi dilateLabelBy]) = scanRAS;
        end   
        %% Creamos los tractos con las esferas
        % Primero hicimos el triangulo sem y perc por separado
        for tckt = tcktype;for conn = connections1;for sphR=ROISphereRadius
            tctname = [numFiberName 'L_' tckt{:} '_' conn{:} '_R' num2str(sphR) '.tck'];
            con1 = conn{:}(1:3);
            coord1 = round(coords.([con1 '_' tckt{:} '_averages' dilateLabelBy])');
            coord1(4) = sphR;
            roi1 = strjoin(arrayfun(@(x) num2str(x),coord1,'UniformOutput',false),',');
            
            con2 = conn{:}(5:7);
            coord2 = round(coords.([con2 '_' tckt{:} '_averages' dilateLabelBy])');
            coord2(4) = sphR;
            roi2 = strjoin(arrayfun(@(x) num2str(x),coord2,'UniformOutput',false),',');
            
            cmd_str = ['tckedit -force -nthreads 16 ' ...
                       '-tck_weights_in ' pDWI fsp subname '/dmri/dti90trilin/mrtrix/' numFiberName 'data_alignedNorm_trilin_noMEC.sift2.weight ' ...
                       '-tck_weights_out ' [mrtrixdir filesep tctname '.sift2.weight'] ' ' ...
                       '-include ' roi1  ' ' ...
                       '-include ' roi2  ' ' ...
                       '-ends_only ' ... 
                         [mrtrixdir filesep WholeTractogramName] ' ' ...
                         [mrtrixdir filesep tctname]];
            bkgrnd = false;
            verbose = false;
            mrtrixVersion = 3;
            AFQ_mrtrix_cmd(cmd_str, bkgrnd, verbose,mrtrixVersion)
            
        end;end;end
        % Luego VOT2VOT y PPC2PPC
        for conn = connections2;for sphR=ROISphereRadius
            tctname = [numFiberName 'L_percsem_' conn{:} '_R' num2str(sphR) '.tck'];
            con1 = conn{:}(1:3);
            coord1 = round(coords.([con1 '_sem_averages' dilateLabelBy])');
            coord1(4) = sphR;
            roi1 = strjoin(arrayfun(@(x) num2str(x),coord1,'UniformOutput',false),',');
            
            con2 = conn{:}(5:7);
            coord2 = round(coords.([con2 '_perc_averages' dilateLabelBy])');
            coord2(4) = sphR;
            roi2 = strjoin(arrayfun(@(x) num2str(x),coord2,'UniformOutput',false),',');
            
            cmd_str = ['tckedit -force -nthreads 16 ' ...
                       '-tck_weights_in ' pDWI fsp subname '/dmri/dti90trilin/mrtrix/' numFiberName 'data_alignedNorm_trilin_noMEC.sift2.weight ' ...
                       '-tck_weights_out ' [mrtrixdir filesep tctname '.sift2.weight'] ' ' ...
                       '-include ' roi1  ' ' ...
                       '-include ' roi2  ' ' ...
                       '-ends_only ' ... 
                         [mrtrixdir filesep WholeTractogramName] ' ' ...
                         [mrtrixdir filesep tctname]];
            bkgrnd = false;
            verbose = false;
            mrtrixVersion = 3;
            AFQ_mrtrix_cmd(cmd_str, bkgrnd, verbose,mrtrixVersion)
        end;end
        toc
    end
    
    
end


%  find -maxdepth 1 -type d -readable -exec sh -c 'echo "$1"; find "$1"/dmri/dti90trilin/mrtrix/2M*.weigth  | wc -l' sh {} ';'

%% Obtener las estadisticas para tractos con SIFT2 
% Lo iba a hacer fuera, pero tengo todas las variables aqui ya...
if(0)    
    resultsFldr = fullfile(DWIdir, 'ResultsCSVSIFT2');
    % mkdir(fullfile(DWIdir, 'ResultsCSVSIFT2'))
    for ns = 1 : length(subs)
        TCTs = { 'L_perc_VOT2PPC_R12.tck',    'L_sem_VOT2PPC_R12.tck', ...
                 'L_perc_VOT2IFG_R12.tck',    'L_sem_VOT2IFG_R12.tck', ...
                 'L_perc_PPC2IFG_R12.tck',    'L_sem_PPC2IFG_R12.tck', ...
                 'L_percsem_VOT2VOT_R12.tck', 'L_percsem_PPC2PPC_R12.tck'};
        numFibers = {'','2M'}
        
        setenv('FREESURFER_HOME', fshome); 
        setenv('SUBJECTS_DIR', fs_SUBJECTS_DIR);
        
        subname = subs(ns).name
        dmridir = fullfile(DWIdir, subname, 'dmri');
        mrtrixdir = fullfile(dmridir, 'dti90trilin','mrtrix');
        % faFile = fullfile(mrtrixdir, 'data_aligned_trilin_noMEC_fa.mif');
        MTVfile = fullfile(qmridir,subname,'OutPutFiles_1','BrainMaps','TV_map.nii.gz');
        T1qfile = fullfile(qmridir,subname,'OutPutFiles_1','BrainMaps','T1_map_Wlin.nii.gz');

        
           
        %% Por cada tracto que hemos leido antes, escribimos stats
        ROISphereRadius = [12]; % 8 y 12 son muy pequeos, igual hasta meter 15
        for nF = numFibers;for tct = TCTs;for sphR=ROISphereRadius
            tctname = [nF{:} tct{:}]
            
            tckFile = [mrtrixdir filesep tctname];
            % csvFA  = fullfile(resultsFldr,['FA_' tctname '_' subname '.csv']);
            csvMTV = fullfile(resultsFldr,['MTV_' tctname '_' subname '.csv']);
            csvT1q = fullfile(resultsFldr,['T1q_' tctname '_' subname '.csv']);
            
            % cmdFA  = ['tcksample -stat_tck mean ' tckFile ' ' faFile ' ' csvFA];
            cmdMTV = ['tcksample -stat_tck mean ' tckFile ' ' MTVfile ' ' csvMTV];
            cmdT1q = ['tcksample -stat_tck mean ' tckFile ' ' T1qfile ' ' csvT1q];

            bkgrnd = false;
            verbose = false;
            mrtrixVersion = 3;
            % AFQ_mrtrix_cmd(cmdFA,  bkgrnd, verbose, mrtrixVersion)
            AFQ_mrtrix_cmd(cmdMTV, bkgrnd, verbose, mrtrixVersion)
            AFQ_mrtrix_cmd(cmdT1q, bkgrnd, verbose, mrtrixVersion)
            
        end;end;end
%         tctname = ['afq_L_vOF.tck'];
%             tckFile = [mrtrixdir filesep tctname];
%             csvFA  = fullfile(resultsFldr,['FA_' tctname '_' subname '.csv']);
%             csvMTV = fullfile(resultsFldr,['MTV_' tctname '_' subname '.csv']);
%             csvT1q = fullfile(resultsFldr,['T1q_' tctname '_' subname '.csv']);
%             
%             cmdFA  = ['tcksample -stat_tck mean ' tckFile ' ' faFile ' ' csvFA];
%             cmdMTV = ['tcksample -stat_tck mean ' tckFile ' ' MTVfile ' ' csvMTV];
%             cmdT1q = ['tcksample -stat_tck mean ' tckFile ' ' T1qfile ' ' csvT1q];
% 
%             bkgrnd = false;
%             verbose = false;
%             mrtrixVersion = 3;
%             AFQ_mrtrix_cmd(cmdFA,  bkgrnd, verbose, mrtrixVersion)
%             AFQ_mrtrix_cmd(cmdMTV, bkgrnd, verbose, mrtrixVersion)
%             AFQ_mrtrix_cmd(cmdT1q, bkgrnd, verbose, mrtrixVersion)
%         
%         tctname = 'afq_L_pAF.tck';
%             tckFile = [mrtrixdir filesep tctname];
%             csvFA  = fullfile(resultsFldr,['FA_' tctname '_' subname '.csv']);
%             csvMTV = fullfile(resultsFldr,['MTV_' tctname '_' subname '.csv']);
%             csvT1q = fullfile(resultsFldr,['T1q_' tctname '_' subname '.csv']);
%             
%             cmdFA  = ['tcksample -stat_tck mean ' tckFile ' ' faFile ' ' csvFA];
%             cmdMTV = ['tcksample -stat_tck mean ' tckFile ' ' MTVfile ' ' csvMTV];
%             cmdT1q = ['tcksample -stat_tck mean ' tckFile ' ' T1qfile ' ' csvT1q];
% 
%             bkgrnd = false;
%             verbose = false;
%             mrtrixVersion = 3;
%             AFQ_mrtrix_cmd(cmdFA,  bkgrnd, verbose, mrtrixVersion)
%             AFQ_mrtrix_cmd(cmdMTV, bkgrnd, verbose, mrtrixVersion)
%             AFQ_mrtrix_cmd(cmdT1q, bkgrnd, verbose, mrtrixVersion)        
        
    end
    
end



% Encontrar y tar: 
% find ./*/dmri/dti90trilin/mrtrix -type f -name *R12.tck -print0 | tar -cvzf  AllSIFT2_Tracts.tar.gz --null -T -

% Econtrar para validar
% find -maxdepth 1 -type d -readable -exec sh -c 'echo "$1"; find "$1"/OutPutFiles_1/BrainMaps/*WM_305.mgh  | wc -l' sh {} ';

