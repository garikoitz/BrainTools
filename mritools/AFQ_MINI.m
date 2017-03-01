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
MINIDIR = '/Users/gari/Documents/BCBL_PROJECTS/MINI';
fsbin = '/Applications/freesurfer/bin';
fshome = '/Applications/freesurfer'; 
% SERVER
% MINIDIR = '/bcbl/home/public/Gari/MINI';
% fsbin = '/opt/freesurfer-5.3.0/freesurfer/bin';
% fshome = '/opt/freesurfer-5.3.0/freesurfer'; 

AnalysisDir = [MINIDIR fsp 'ANALYSIS'];
fs_SUBJECTS_DIR = fullfile(AnalysisDir, 'freesurferacpc');
DWIdir  = fullfile(AnalysisDir, 'DWI');
cd(DWIdir);
subs = dir('S*');
retDIR = fullfile(AnalysisDir, 'ret');
fMRIDIR = fullfile(AnalysisDir, 'fMRI_SPM', 'block', 'data');

      





%%%%%   CLUSTER PARPOOL    %%%%%%
% % myclusterLocal = parcluster('local');
% % myclusterLocal.NumWorkers
% [st, re] = system('qstat -g c | grep matlab.q');
% [Tok, Rem] = strtok(re);
% [Tok, Rem] = strtok(Rem);
% [Tok, Rem] = strtok(Rem);
% [Tok, Rem] = strtok(Rem);
% [available] = strtok(Rem)
% parpool('ips_base', str2num(available))
%%%%% END CLUSTER PARPOOL  %%%%%%



%% Encontrar VOF PARC y pasarlos a surfaces en fsaverage
if(0)
    for ns = 1 : length(subs)
    subname = subs(ns).name

    setenv('FREESURFER_HOME', fshome); 

    %% Find the VOF per every subject

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

    wholebrainfgPath= fullfile(dmridir, 'dti90trilin', 'fibers'); 
    wholebrainfg= fullfile(wholebrainfgPath, 'WholeBrainFG.mat'); 
    fgMori = dtiReadFibers(fullfile(wholebrainfgPath, 'MoriGroups.mat'));
    L_arcuate= fgMori(19);
    R_arcuate= fgMori(20);
    fsROIdir= outDir;
    outdir = fullfile(wholebrainfgPath,'VOF');
    if ~exist(outdir, 'dir'), mkdir(outdir), end
    thresh= [];
    v_crit= [];
    dt= dtiLoadDt6(fullfile(dmridir, 'dti90trilin', 'dt6.mat'));
    savefiles= true;
    arcThresh= [];
    parcThresh= [];
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
    % %% RENDER
    % % Read the ROIs in the cortex
    roi_L_fusiform = dtiReadRoi(fullfile(dmridir,'ROIs', ...
                                '1007_ctx-lh-fusiform.mat'));
    roi_L_inferiortemporal = dtiReadRoi(fullfile(dmridir,'ROIs', ...
                                             '1009_ctx-lh-inferiortemporal.mat'));
    roi_L_lateraloccipital = dtiReadRoi(fullfile(dmridir,'ROIs', ...
                                             '1011_ctx-lh-lateraloccipital.mat'));                                     
    
    
    % Rnder them the fibers and the cortex rois
    AFQ_RenderFibers(L_VOF , 'color',  [158 47 88]/256, ...
                    'tubes',[0]); % Render the fibers
    AFQ_RenderFibers(L_pArc , 'color',  [237 139 140]/256, 'tubes',[0], ...
                     'newfig',false); % Render the fibers
    % Render the roi in orange
    AFQ_RenderRoi(roi_L_fusiform, [241 217 201]/256, 'mesh'); 
    AFQ_RenderRoi(roi_L_inferiortemporal, [241 217 181]/256, 'mesh'); 
    AFQ_RenderRoi(roi_L_lateraloccipital, [241 217 161]/256, 'mesh'); 
% 
%     % con este comando binarizo aparc+aseg y ademas solo me quedo con el GM
    FSLDIR = '/opt/fsl/fsl-5.0.9/fsl';
    setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
    setenv('FSLDIR', FSLDIR);
    system([FSLDIR fsp 'bin' fsp 'fslmaths ' dmridir fsp 'aparc+aseg.nii.gz '...
            '-thr 1000 -bin ' dmridir fsp 'segmentation.nii.gz']);
%     % Y ahora creo el mesh
%     im = [dmridir fsp 'segmentation.nii.gz'];
%     % msh = AFQ_meshCreate(im, 'color', [.8 .7 .6])
% 
% 
% 
% 
    % Obtener el file con el intersect entre fiber y cortex
    % por ahora solo hago el posterior arcuate y el VOF
    fiberRois = {L_VOF, L_pArc, L_pArc_vot};
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
        MRIwrite(segRead, [dmridir fsp fg.name '_tracts.nii.gz']);


        % Hay que pensar si hago el paso a la superficie con todos los
        % voxeles que pertenecen a los tractos, o solo me quedo con
        % aquellos voxeles que coinciden con los ROI de interes y luego
        % hago el paso a la superifice. >> He pasado todos los voxeles,
        % luego con aparc podre elegir los voxeles que me interesen para
        % los rois. Lo de los tractos tiene que ser bidireccional. 

        % Y ahora los convertimos a superficie usando fs
        movname    = fullfile(dmridir, [fg.name '_tracts.nii.gz']);
        oname      = fullfile(dmridir, [fg.name '_tracts.mgh']);
        oname305   = fullfile(dmridir, [fg.name '_tracts305.mgh']);

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
if(1)
    for ns = 1 : length(subs)
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
if(1)
    ROIs = {'PPC_perc_averages', 'PPC_sem_averages', ...
             'VOT_perc_averages', 'VOT_sem_averages', ...
             'IFG_perc_averages', 'IFG_sem_averages'};
    dilateLabelBy = '1';
    setenv('FREESURFER_HOME', fshome); 
    setenv('SUBJECTS_DIR', fs_SUBJECTS_DIR);
    for ns = 1 : length(subs)
        subname = subs(ns).name

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
if(1)    
    ROIs = {'PPC_perc_averages', 'PPC_sem_averages', ...
             'VOT_perc_averages', 'VOT_sem_averages', ...
             'IFG_perc_averages', 'IFG_sem_averages'};
    tcktype       = {'sem', 'perc'};
    connections   = {'VOT2IFG', 'VOT2PPC', 'PPC2IFG'};
    dilateLabelBy = '1';
    setenv('FREESURFER_HOME', fshome); 
    setenv('SUBJECTS_DIR', fs_SUBJECTS_DIR);
    ROISphereRadius = [8,10,12];
    WholeTractogramName = 'data_aligned_trilin_noMEC_wmCsd_lmax4_data_aligned_trilin_noMEC_wmMask_data_aligned_trilin_noMEC_wmMask_iFOD2-500000.tck';
    for ns = 1 : length(subs)
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
        for tckt = tcktype;for conn = connections;for sphR=ROISphereRadius
            tctname = ['L_' tckt{:} '_' conn{:} '_R' num2str(sphR) '.tck'];
            con1 = conn{:}(1:3);
            coord1 = round(coords.([con1 '_' tckt{:} '_averages' dilateLabelBy])');
            coord1(4) = sphR;
            roi1 = strjoin(arrayfun(@(x) num2str(x),coord1,'UniformOutput',false),',');
            
            con2 = conn{:}(5:7);
            coord2 = round(coords.([con2 '_' tckt{:} '_averages' dilateLabelBy])');
            coord2(4) = sphR;
            roi2 = strjoin(arrayfun(@(x) num2str(x),coord2,'UniformOutput',false),',');
            
            cmd = ['tckedit ' ...
                   '-include ' roi1  ' ' ...
                   '-include ' roi2  ' ' ...
                   '-ends_only ' ... 
                   [mrtrixdir filesep WholeTractogramName] ' ' ...
                   [mrtrixdir filesep tctname]];
            system(cmd)
        end;end;end
    end      
end



