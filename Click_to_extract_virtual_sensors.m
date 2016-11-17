% Plot SPM MEEG source localised data on a glass brain.
% Includes ability to also project MNI coordinates (M) for visualisation 
% Can be used to select [click] region of interest over a set of subject,
% find local corresponding vertices and extract full datasets from here
% 
% see also plotmesh and plotmesh_fo
% AS2016

% input cell array of spm object filenames
f = dir('ica_*.mat'); f = {f.name}';


t = 1;               % trial / condition number
r = 5;               % radius to search

M = [...             % MMN MNI regions to project for guidance
   -46    20     8;
   -61   -32     8;
   -42   -14     7;
    46    20     8;
    59   -25     8;
    46   -14     8];

M = [-7 -103 -4]; %  ~ V1


woi  = [];
foi  = [];
type = [];

plotmesh_fo_grp(f,M,t,woi,foi,type); % projection with mean 'activations' over datsets
[rx,ry,rz] = get_xyz();              % get cursor click location [turn off rotation first]

for i = 1:length(f)
    
    [rx ry rz]            = ind2mni(f{i},rx,ry,rz);      % find nearest vertex to selected point for each subject
    [dD{i},C,T,XYZ{i},FT] = MNI2TS({f{i}},[rx ry rz],r); % extract full dataset [samples by trials] for each identified vertex

end




[-46    20     8]
[-61   -32     8]
[-42   -14     7]
[46    20     8]
[59   -25     8]
[46   -14     8]
