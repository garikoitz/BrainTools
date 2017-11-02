function [] = mainSplitShells(DWIdir, subName, shell)
% normalizes the bvals and splits the bvecs
if ~exist([DWIdir filesep subName filesep 'raw' num2str(shell)])
	        mkdir([DWIdir filesep subName filesep 'raw' num2str(shell)])
end
cd ([DWIdir filesep subName filesep 'raw' num2str(shell)])

% if isempty(getenv('SERVICE_DIR'))
%     disp('setting SERVICE_DIR to pwd')
%     setenv('SERVICE_DIR', pwd)
% end
% 
% switch getenv('ENV')
%     case 'IUHPC'
%         disp('loading paths (HCP)')
%         addpath(genpath('/N/u/hayashis/BigRed2/git/jsonlab'));
%         addpath(genpath('/N/u/hayashis/BigRed2/git/vistasoft'));
%     case 'VM'
%         disp('loading paths (VM)')
%         addpath(genpath('/usr/local/jsonlab'))
%         addpath(genpath('/usr/local/vistasoft'))
% end

% load config.json
% config = loadjson('config.json');
% (GLU: I prefer to have it here copy-pasted, the json is for the Docker file)
%{
%        "shell": 2000,
%        "bvals": "/input/dwi.bvals",
%        "bvecs": "/input/dwi.bvecs",
%        "dwi": "/input/dwi.nii.gz"
%}
config.shell  = shell;
config.bvals  = '../raw/dwi.bvals';
config.bvecs  = '../raw/dwi.bvecs';
config.dwi    = '../raw/dwi.nii.gz';


% Parameters used for normalization
params.shells       = config.shell;
bvals = dlmread(config.bvals);
bvecs = dlmread(config.bvecs);
dwi = niftiRead(config.dwi);

for i = 1:length(params.shells)
    index = (bvals == params.shells(i));
    index0 = (bvals == 0);
    all_index = or(index, index0);
    assertEqual(sum(all_index), sum(index0) + sum(index));
    
    %write files
    %dlmwrite('dwi.bvals',bvals.valnorm,'delimiter',' ');
    dlmwrite(sprintf('dwi.bvals', params.shells(i)), bvals(all_index), 'delimiter',' ');
    dlmwrite(sprintf('dwi.bvecs', params.shells(i)), bvecs(:,all_index));
    dwi_oneshell = dwi;
    dwi_oneshell.fname = sprintf('dwi.nii.gz', params.shells(i));
    size(dwi.data)
    size(all_index)
    dwi_oneshell.data = dwi.data(:,:,:,all_index);
    dwi_oneshell.dim(4) = size(dwi_oneshell.data,4);
    niftiWrite(dwi_oneshell);
end


