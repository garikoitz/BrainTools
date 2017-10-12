function png2avi(pngdir,filestr, outname)
% Edited by GLU to be able to sort by filename



if ~exist('filestr','var') || isempty(filestr)
    filestr = '*.png';
end
d = dir(fullfile(pngdir,filestr));


% Sort files by creation date
[tmp ind]=sort({d.date});
d=d(ind);

d = d(1:20);

d(:).name

for ii = 1:length(d)
   im(:,:,:,ii) = imread(fullfile(pngdir,d(ii).name));
end
v = VideoWriter(outname);
open(v);
writeVideo(v,im);
close(v);



%% Creando videos para la defensa de la tesis

% LOTS
% pngdir = '/Users/gari/Documents/BCBL_PROJECTS/MINI/ANALYSIS/myGLMFIT/vMC3_sem_MTV_LOTS_behav_fmri_block_RWvsPW/clusterResult13PNGsoloVOT';
% filestr = 'CON_DWI-Rsq_fMRIPPC*.fhmw5.png';
% outname = '/Users/gari/gDrive/BCBL/PROYECTOS/MINI/_PUBLISH_/Defensa/vMC3_sem_MTV_LOTS_behav_fmri_block_RWvsPW.avi';
% 
% pngdir = '/Users/gari/Documents/BCBL_PROJECTS/MINI/ANALYSIS/myGLMFIT/vMC2_perc_MTV_LOTS_behav_fmri_block_RWvsCB/clusterResult13PNGsoloVOT';
% filestr = '*.png';
% outname = '/Users/gari/gDrive/BCBL/PROYECTOS/MINI/_PUBLISH_/Defensa/vMC2_perc_MTV_LOTS_behav_fmri_block_RWvsCB.avi';
% 
% % IFG
% pngdir = '/Users/gari/Documents/BCBL_PROJECTS/MINI/ANALYSIS/myGLMFIT/vMC2_perc_MTV_IFG_behav_fmri_block_RWvsCB/clusterResult13PNGsoloVOT';
% filestr = '*.png';
% outname = '/Users/gari/gDrive/BCBL/PROYECTOS/MINI/_PUBLISH_/Defensa/vMC2_perc_MTV_IFG_behav_fmri_block_RWvsCB.avi';
% 
% pngdir = '/Users/gari/Documents/BCBL_PROJECTS/MINI/ANALYSIS/myGLMFIT/vMC3_sem_MTV_IFG_behav_fmri_block_RWvsPW/clusterResult13PNGsoloVOT';
% filestr = '*.png';
% outname = '/Users/gari/gDrive/BCBL/PROYECTOS/MINI/_PUBLISH_/Defensa/vMC3_sem_MTV_IFG_behav_fmri_block_RWvsPW.avi';
