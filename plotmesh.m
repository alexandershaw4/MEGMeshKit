function plotmesh(D,o)
% Plot glass mesh brain with MNI coordinates marked
% from source localised SPM MEEG object
%
%
% AS2016

%for i = 1:2; subplot(1,2,i),plotmesh(D,MNI);end

% verts & faces for brain
vert = D.inv{end}.forward(end).mesh.vert;
x = vert(:,1);
y = vert(:,2);
z = vert(:,3);

face  = D.inv{end}.forward(end).mesh.face;

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

%camorbit(280,270,'data');

set(gcf,'inverthardcopy','off');

hold on;


% find vertices corresponding to MNIs
XYZ   = o;
inv   = D.inv{end};
rad   = 1;
s     = 500; % size of marker

Ns    = size(XYZ, 1); % n points to plot
svert = {};
for i = 1:Ns
    dist = sqrt(sum([vert(:,1) - XYZ(i,1), ...
                     vert(:,2) - XYZ(i,2), ...
                     vert(:,3) - XYZ(i,3)].^2, 2));
    if rad > 0
        %svert{i} = find(dist < rad);
        for j = 1:rad
            [junk,svert{i,j}] = min(dist);
            dist(svert{i,j}) = NaN;
            %XYZ(i,j,:) = vert(svert{i,j});
        end
    else
        [junk,svert{i}] = min(dist);
        XYZ(i, :) = vert(svert{i}, :);
    end
end

% add thse to plot
for i = 1:Ns
    for j = 1:size(svert,2)
        scatter3(vert(svert{i,j},1),vert(svert{i,j},2),vert(svert{i,j},3),s,i,'r','filled');
    end
end


% MNI=[-46 20 8;
%  -61 -32 8;
%  -42 -22 7;
%   46 20 8;
%   59 -25 8;
%   46 -14 8];






