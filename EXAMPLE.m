% Example meshkit script


% add tools
addpath('~/Downloads/meshkit/');
cd ~/faces/ ;           % wd


f   = loadarrayspm(F)'; % F is a cell array of MEEG filenames
toi = [0 .35];          % time window of interest
foi = [];               % freqyency window


% plot options
%------------------------------------------------
global inflate ; inflate = .04;   % inflate
global thr     ; thr     = [];    % threshold for making t(+/- n) = 0
global trs     ; trs     = .07;   % rescale theshold

trans = .8;       % overlay alpha (0 - 1)
typ   = 'trials'; % trials or evoked


% Calculate some group responses - here there are conditions called
% 'NeutralAll' and 'HappyAll' in D.condlist which have been source
% localised in spm

Neut  = plotmesh_fo_grp(f,[],{'NeutralAll'},toi,foi,typ,[],trans);
Happ  = plotmesh_fo_grp(f,[],{'HappyAll'},toi,foi,typ,[],trans);

% I also want to average over these 2 conditions:
All   = plotmesh_fo_grp(f,[],{'NeutralAll',...
                              'HappyAll'},toi,foi,typ,[],trans);

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

