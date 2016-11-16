function plotmesh_fo(D,o,t)
% Plot glass mesh brain with MNI coordinates marked
% from source localised SPM MEEG object
%
%
% AS2016

% options
%---------------------------------------------------------
dS    = 100; % dot size for functional overlay of trial{t}
s     = 500; % patch size for MNI coordinates
CL    = 'r'; % colour of MNI patch

woi   = [0 .3];   % time window if interest for source data in trial{t}
foi   = [];       % freq window if interest for source data in trial{t}
type  = 'evoked'; % 'evoked', 'induced' or 'trial'


%for i = 1:2; subplot(1,2,i),plotmesh(D,MNI);end



% verts & faces for brain
%---------------------------------------------------------
vert  = D.inv{end}.forward(end).mesh.vert;
x     = vert(:,1);
y     = vert(:,2);
z     = vert(:,3);
face  = D.inv{end}.forward(end).mesh.face;

% glass brain
%---------------------------------------------------------
h = patch('faces',face,'vertices',[x(:) y(:) z(:)]);
set(h,'FaceColor',[.4 .4 .4]);
box off;
grid off;
whitebg(1,'w'); 
camlight('right')
axis tight
set(h,'EdgeColor','none')
material dull
alpha(.2);
lighting phong
set(gca,'visible','off');
set(gcf,'inverthardcopy','off');
hold on;

% call function for projecting into source space
%---------------------------------------------------------
FO   = rebuild(D,woi,type,foi);
t    = FO.JW{t};
scatter3(x,y,z,[],t,'filled');alpha(.3)


% find vertices corresponding to provided MNIs
%---------------------------------------------------------
XYZ   = o;
inv   = D.inv{end};
rad   = 1;   

Ns    = size(XYZ, 1); % n points to plot
svert = {};
for i = 1:Ns
    dist = sqrt(sum([vert(:,1) - XYZ(i,1), ...
                     vert(:,2) - XYZ(i,2), ...
                     vert(:,3) - XYZ(i,3)].^2, 2));
    if rad > 0
        for j = 1:rad
            [junk,svert{i,j}] = min(dist);
            dist(svert{i,j}) = NaN;
        end
    else
        [junk,svert{i}] = min(dist);
        XYZ(i, :) = vert(svert{i}, :);
    end
end

% add selected points to plot
%---------------------------------------------------------
for i = 1:Ns
    for j = 1:size(svert,2)
        scatter3(vert(svert{i,j},1),...
                 vert(svert{i,j},2),...
                 vert(svert{i,j},3),...
                 s,CL,'filled');
                    alpha(.2)
    end
end





% MNI=[-46 20 8;
%  -61 -32 8;
%  -42 -14 7; %-22
%   46 20 8;
%   59 -25 8;
%   46 -14 8];






