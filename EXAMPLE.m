% Example meshkit script


% add tools
addpath('~/Downloads/meshkit/');
cd ~/faces/ ;           % wd


f   = loadarrayspm(F)'; % F is a cell array of MEEG filenames


global inflate ; inflate = .04;   % inflate mesh


% Calculate some group responses - here there are conditions called
% 'NeutralAll' and 'HappyAll' in D.condlist which have been source
% localised in spm

cfg.woi  = [0 .2]; % cfg optional - see help mesh_pca1
cfg.neig = 20; 

[Neutm,Neut]   = plotmesh_fo_grp_pca(f,{'NeutralAll'},cfg);
[Happm,Happ]   = plotmesh_fo_grp_pca(f,{'HappyAll'},cfg);
[Angrm,Angr]   = plotmesh_fo_grp_pca(f,{'AngryAll'},cfg);
[FTAllm,FTAll] = plotmesh_fo_grp_pca(f,{'FTAll'},cfg);

[Allm,All]     = plotmesh_fo_grp_pca(f,{'NeutralAll',...
                              'HappyAll',...
                              'AngryAll',...
                              'FTAll'},cfg);

crit  = 5.5025e-06; % if using t-image, use a critical t on colbar
alph  = .2;         % alpha value for overlay

% New: optionally write these [group aves] g/nifti meshs:
writeims(f{1},abs(All),'AllFcGam');


% Plots
%-----------------------------------
plotmeshfo(F{1},abs(All),.05,alph); % plot on mesh with some alpha
[x,y,z] = get_xyz(1,F{1});          % click to get peak coords

% plot on mesh with click coords marked
plotmeshfo(F{1},abs(All),.05,alph,[x' y' z']); 

% convert these vertex coords to subject MNIs
for i = 1:length(F)
    MNI(i,:,:) = mesh2mni(F{i},[x' y' z']);
end

% More mesh plots
%-----------------
figure,
subplot(121),plotmeshfo(F{1},abs(Neut),crit,alph);
title('Neutral');set(findall(gca, 'type', 'text'), 'visible', 'on');

subplot(122),plotmeshfo(F{1},abs(Happ),crit,alph);
title('Happy');set(findall(gca, 'type', 'text'), 'visible', 'on');

