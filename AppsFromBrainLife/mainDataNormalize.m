function [] = mainDataNormalize(DWIdir, subName)

if ~exist([DWIdir filesep subName filesep 'raw'])
	mkdir([DWIdir filesep subName filesep 'raw'])
end
cd ([DWIdir filesep subName filesep 'raw'])
% if isempty(getenv('SERVICE_DIR'))
%     disp('setting SERVICE_DIR to pwd')
%     setenv('SERVICE_DIR', pwd)
% end
% 
% switch getenv('ENV')
%     case 'IUHPC'
%         disp('loading paths (HPC)')
%         addpath(genpath('/N/u/hayashis/BigRed2/git/jsonlab'))
%     case 'VM'
%         disp('loading paths (VM)')
%         addpath(genpath('/usr/local/jsonlab'))
% end
% 
% normalizes the bvals and flips the bvecs

% load config.json
% config = loadjson('config.json');
% GLU: instead of reading it, copy it here.
% {
%         "bvals": "/input/dwi.bvals",
%         "bvecs": "/input/dwi.bvecs",
%         "dwi": "/input/dwi.nii.gz",
%         "xflip": true,
%         "yflip": false,
%         "zflip": false
% }
config.bvals = [DWIdir filesep subName filesep 'T1w/Diffusion_7T/bvals'];
config.bvecs = [DWIdir filesep subName filesep 'T1w/Diffusion_7T/bvecs'];
config.dwi   = [DWIdir filesep subName filesep 'T1w/Diffusion_7T/data.nii.gz'];
config.xflip = true; 
config.yflip = false;
config.zflip = false;



system(sprintf('ln -s %s %s', config.dwi, 'dwi.nii.gz'))

% Parameters used for normalization
params.thresholds.b0_normalize    = 200;
params.thresholds.bvals_normalize = 100;

%% Normalize HCP files to the VISTASOFT environment

bvals.val = dlmread(config.bvals);

% Round the numbers to the closest thousand 
% This is necessary because the VISTASOFT software does not handle the B0
% when they are not rounded.
[bvals.unique, ~, bvals.uindex] = unique(bvals.val);

bvals.unique(bvals.unique <= params.thresholds.b0_normalize) = 0;
bvals.unique  = round(bvals.unique./params.thresholds.bvals_normalize) ...
    *params.thresholds.bvals_normalize;
bvals.valnorm = bvals.unique( bvals.uindex );
dlmwrite('dwi.bvals',bvals.valnorm,'delimiter',' ');

%% Flip the Bvecs on chosen axis

%load bvecs
bvecs = dlmread(config.bvecs);

%params
params.x_flip = config.xflip;
params.y_flip = config.yflip;
params.z_flip = config.zflip;

if ~(size(bvecs,1) == 3), bvecs = bvecs'; end

if params.x_flip
    bvecs(1,:) = -bvecs(1,:);
end
if params.y_flip
    bvecs(2,:) = -bvecs(2,:);
end
if params.z_flip
    bvecs(3,:) = -bvecs(3,:);
end

%savejson('', config,'products.json')

dlmwrite('dwi.bvecs',bvecs,'delimiter',' ');

end


