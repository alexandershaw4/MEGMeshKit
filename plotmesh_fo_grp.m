function plotmesh_fo_grp(D,o,t,woi,foi,type,CL)
% Plot glass mesh brain with MNI coordinates marked
% from source localised SPM MEEG objects
%
% D is a cell array of subjects' D filenames
%
% AS2016

% options
%---------------------------------------------------------
dS    = 100; % dot size for functional overlay of trial{t}
s     = 500; % patch size for MNI coordinates

try CL; catch; CL    = 'r'; end%[.4 .4 .4];%'r'; % colour of MNI patch

try woi;   catch; woi = [-.1 .3]; end  % time window if interest for source data in trial{t}
try foi;   catch; foi = [];       end  % freq window if interest for source data in trial{t}
try type;  catch; type  = 'evoked';end % 'evoked', 'induced' or 'trial'

%for i = 1:2; subplot(1,2,i),plotmesh(D,MNI);end

if iscell(D); D = loadarrayspm(D);
else return;
end

% verts & faces for brain
%---------------------------------------------------------
vert  = D{1}.inv{end}.forward(end).mesh.vert;
x     = vert(:,1);
y     = vert(:,2);
z     = vert(:,3);
face  = D{1}.inv{end}.forward(end).mesh.face;

% glass brain
%---------------------------------------------------------
h = patch('faces',face,'vertices',[x(:) y(:) z(:)]);
set(h,'FaceColor',[.4 .4 .4]);
box off;
grid off;
%whitebg(1,'w'); 
camlight('right')
%axis tight
set(h,'EdgeColor','none')
material dull
alpha(.2);
lighting phong
set(gca,'visible','off');
set(gcf,'inverthardcopy','off');
hold on;

% call function for projecting into source space
%---------------------------------------------------------
if isempty(woi);  woi = [0 .3]; end
if isempty(type); type = 'evoked'; end
for SUB = 1:length(D)
    if SUB > 1; fprintf(repmat('\b',[size(str),1])); end
    str = sprintf('Fetching projections for %d of %d datasets\n',SUB,length(D));
    fprintf(str);

    FO        = rebuild(D{SUB},woi,type,foi);
    it        = FO.JW{t};
    st(SUB,:) = it;
end

%st = PEig(full(st'));
mst = mean(st,1);
%cbr = round(TSNorm(st,6,2));
%cbr = sum(cbr,1);
%cbr = cbr-min(cbr);

%scatter3(x,y,z,[],cbr,'filled');alpha(.3);
scatter3(x,y,z,[],mst(:),'filled');alpha(.3);


% Discard other datasets now:
D = D{1};

if ~isempty(o) % MNIs
    
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

end




% MNI=[-46 20 8;
%  -61 -32 8;
%  -42 -14 7; %-22
%   46 20 8;
%   59 -25 8;
%   46 -14 8];






