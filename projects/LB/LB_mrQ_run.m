function LB_mrQ_run(session,startsub,proclus)
if notDefined('session')
    session=0;
end
if notDefined('startsub')
    startsub = 1;
end
if notDefined('proclus')
    proclus = 1;
end
if session ==1
    rawDirs = {'/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session1/LB1_20130630/20130630_1437'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session1/LB2_20130628/20130628_1812'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session1/LB11_20130709/20130709_1008'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session1/LB17_20130728/20130728_1345'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session1/LB18_20130805/20130805_1025'};
    refIms = {'/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB1/t1.nii.gz'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB2/t1.nii.gz'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB11/t1.nii.gz'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB17/t1.nii.gz'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB18/t1.nii.gz'};
elseif session ==2
    rawDirs = {'/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session2/LB1_20130716/20130716_1606'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session2/LB2_20130715/20130715_1805'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session2/LB11_20130731/20130731_1805'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session2/LB17_20130813/20130813_1910'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session2/LB18_20130820/20130820_1438'};
    refIms = {'/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB1/t1.nii.gz'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB2/t1.nii.gz'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB11/t1.nii.gz'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB17/t1.nii.gz'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB18/t1.nii.gz'};
elseif session ==3
    rawDirs = {'/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session3/LB1_20130730/20130730_1004'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session3/LB2_20130729/20130729_1738'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session3/LB11_20130819/20130819_1012'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session3/LB17_20130827/20130827_1835'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session3/LB18_20130904/20130904_1540'};
    refIms = {'/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB1/t1.nii.gz'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB2/t1.nii.gz'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB11/t1.nii.gz'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB17/t1.nii.gz'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB18/t1.nii.gz'};
elseif session ==4
    rawDirs = {'/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session4/LB1_20130818/20130818_1602'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session4/LB2_20130808/20130808_1746'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session4/LB11_20130909/20130909_1617'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session4/LB17_20130910/20130910_1907'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session4/LB18_20131016/20131016_1443'}
    refIms = {'/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB1/t1.nii.gz'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB2/t1.nii.gz'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB11/t1.nii.gz'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB17/t1.nii.gz'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB18/t1.nii.gz'};
elseif session == 0
    rawDirs = {'/biac4/wandell/data/Lindamood_Bell/MRI/child/control/session1/LB7_20130712/20130712_1717'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/child/control/session1/LB8_20130723/20130723_1304'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/child/control/session1/LB9_20130711/20130711_1616'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/child/control/session1/LB10_20130703/20130703_1148'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/child/control/session1/LB12_20130717/20130717_1743'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/child/control/session1/LB15_20130729/20130729_1325'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/child/control/session1/LB16_20130722/20130722_1744'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/child/control/session2/LB7_20130802/20130802_1028'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/child/control/session2/LB8_20130813/20130813_1214'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/child/control/session2/LB9_20130807/20130807_1547'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/child/control/session2/LB10_20130717/20130717_1334'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/child/control/session2/LB12_20130806/20130806_1315'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/child/control/session2/LB15_20131017/20131017_1635'...
        }
    refIms = {'/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB7/t1.nii.gz'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB8/t1.nii.gz'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB9/t1.nii.gz'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB10/t1.nii.gz'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB12/t1.nii.gz'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB15/t1.nii.gz'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB16/t1.nii.gz'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB7/t1.nii.gz'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB8/t1.nii.gz'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB9/t1.nii.gz'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB10/t1.nii.gz'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB12/t1.nii.gz'...
        '/biac4/wandell/data/Lindamood_Bell/MRI/anatomy/LB15/t1.nii.gz'...
        }
end

for ii = 1:length(rawDirs)
    outdirs{ii} = fullfile(rawDirs{ii},'mrQ_revised');
end
for ii = startsub:length(rawDirs)
    if ~exist(outdirs{ii},'dir'); mkdir(outdirs{ii}); end
    mrQ = mrQ_Create(rawDirs{ii},[],outdirs{ii});
    
    % Set other parameters
    if ~isfield(mrQ,'sub') || isempty(mrQ.sub)
        mrQ = mrQ_Set(mrQ,'sub',[num2str(session) '_' num2str(ii)]);
    end
    mrQ = mrQ_Set(mrQ,'proclus',proclus);
    %mrQ = mrQ_Set(mrQ,'sungrid',1);
    mrQ = mrQ_Set(mrQ,'fieldstrength',3);
    
    if ~isfield(mrQ,'refIm') || isempty(mrQ.refIm)
        % New input to automatically acpc align
        mrQ = mrQ_Set(mrQ,'refim',refIms{ii});
        
        % Specific arrange function for nimsfs
        mrQ = mrQ_arrangeData_nimsfs(mrQ);
    end
    % RUN IT
    mrQ_run(mrQ.name);
end

return

%% register maps to dti
sub_dirs = {'/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session1/LB1_20130630/20130630_1437/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session2/LB1_20130716/20130716_1606/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session3/LB1_20130730/20130730_1004/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session4/LB1_20130818/20130818_1602/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session1/LB2_20130628/20130628_1812/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session2/LB2_20130715/20130715_1805/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session3/LB2_20130729/20130729_1738/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session4/LB2_20130808/20130808_1746/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session1/LB4_20130807/20130807_1120/raw/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session2/LB4_20130906/20130906_1538/raw/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session3/LB4_20130927/20130927_1512/raw/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session4/LB4_20131120/20131120_1054/raw/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session1/LB11_20130709/20130709_1008/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session2/LB11_20130731/20130731_1805/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session3/LB11_20130819/20130819_1012/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session4/LB11_20130909/20130909_1617/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session1/LB17_20130728/20130728_1345/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session2/LB17_20130813/20130813_1910/dti70trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session3/LB17_20130827/20130827_1835/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session4/LB17_20130910/20130910_1907/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session1/LB18_20130805/20130805_1025/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session2/LB18_20130820/20130820_1438/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session3/LB18_20130904/20130904_1540/dti111trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session4/LB18_20131016/20131016_1443/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/control/session1/LB7_20130712/20130712_1717/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/control/session2/LB7_20130802/20130802_1028/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/control/session1/LB8_20130723/20130723_1304/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/control/session2/LB8_20130813/20130813_1214/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/control/session1/LB9_20130711/20130711_1616/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/control/session2/LB9_20130807/20130807_1547/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/control/session1/LB10_20130703/20130703_1148/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/control/session2/LB10_20130717/20130717_1334/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/control/session1/LB12_20130717/20130717_1743/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/control/session2/LB12_20130806/20130806_1315/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/control/session1/LB15_20130729/20130729_1325/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/control/session2/LB15_20131017/20131017_1635/dti80trilin'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/control/session1/LB16_20130722/20130722_1744/dti80trilin'};
mapDirs = {'/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session1/LB1_20130630/20130630_1437/mrQ_revised/OutPutFiles_1/BrainMaps'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session2/LB1_20130716/20130716_1606/mrQ_revised/OutPutFiles_1/BrainMaps'...
    []...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session4/LB1_20130818/20130818_1602/mrQ_revised/OutPutFiles_1/BrainMaps'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session1/LB2_20130628/20130628_1812/mrQ_revised/OutPutFiles_1/BrainMaps'...
    []...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session3/LB2_20130729/20130729_1738/mrQ_revised/OutPutFiles_1/BrainMaps'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session4/LB2_20130808/20130808_1746/mrQ_revised/OutPutFiles_1/BrainMaps'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session1/LB4_20130807/20130807_1120/SPGR_2/Align_0.9375_0.9375_1.5/maps'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session2/LB4_20130906/20130906_1538/SPGR_2/Align_0.9375_0.9375_1.5/maps'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session3/LB4_20130927/20130927_1512/SPGR_1/Align_0.9375_0.9375_1.5/maps'...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session4/LB4_20131120/20131120_1054/SPGR_2/Align_0.9375_0.9375_1.5/maps'...
    []...
    []...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session3/LB11_20130819/20130819_1012/mrQ_revised/OutPutFiles_1/BrainMaps'...
    []...
    []...
    []...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session3/LB17_20130827/20130827_1835/mrQ_revised/OutPutFiles_1/BrainMaps'...
    []...
    []...
    []...
    '/biac4/wandell/data/Lindamood_Bell/MRI/child/intervention/session3/LB18_20130904/20130904_1540/mrQ_revised/OutPutFiles_1/BrainMaps'...
    []}
% Warp maps to dti
for ii = 1:length(mapDirs)
    if ~exist(fullfile(sub_dirs{ii},'bin','T1_map_lsq_2DTI.nii.gz'),'file');
        t1File = fullfile(mapDirs{ii},'T1_map_lsq.nii.gz');
        b0File = fullfile(sub_dirs{ii},'bin', 'b0.nii.gz');
        otherMaps = {fullfile(mapDirs{ii},'SIR_map.nii.gz')...
            fullfile(mapDirs{ii},'TV_map.nii.gz')...
            fullfile(mapDirs{ii},'VIP_map.nii.gz')...
            fullfile(mapDirs{ii},'WF_map.nii.gz')};
        if exist(t1File,'file')
            mrQ_registerMap2DTI(b0File,t1File,otherMaps,fullfile(sub_dirs{ii},'bin'));
            % make an R1 map from the t1 map
            t1 = readFileNifti(fullfile(sub_dirs{ii},'bin','T1_map_lsq_2DTI.nii.gz'));
            t1.data(t1.data>4) = 4;
            t1.data(t1.data<.4) = .4;
            t1.data=1./t1.data;
            t1.fname = fullfile(sub_dirs{ii},'bin','R1_map_lsq_2DTI.nii.gz');
            writeFileNifti(t1); clear t1;
        else
            fprintf('\nT1 map does not exist')
        end
    end
end