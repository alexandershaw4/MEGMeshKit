function plotmesh_fo_tmap(D,tmap,thr,trs,trans)
% Plot glass mesh brain with MNI coordinates marked
% from source localised SPM12 MEEG object
%
% tmap - from contrast
% thr  - threshold for overlay [e.g. critical t val]
% trs  - blanking values [remove colour in overlay]
% trans - transparency
%
% AS2016

% options
%---------------------------------------------------------
dS    = 100; % dot size for functional overlay of trial{t}
s     = 500; % patch size for MNI coordinates
CL    = 'r'; % colour of MNI patch

try CC;   catch;CL    = 'r'; end
try trans;catch;trans = .4; end; if isempty(trans); trans=.4; end

warning off
if isempty(thr); global thr ; end
if isempty(trs); global trs ; end
global flt
global inflate

% verts & faces for brain
%---------------------------------------------------------
vert  = D.inv{end}.forward(end).mesh.vert;
x     = vert(:,1);
y     = vert(:,2);
z     = vert(:,3);
face  = D.inv{end}.forward(end).mesh.face;

% fix unspecified params
if isempty(flt);    flt = []; end
if isempty(inflate);inflate = [];end


% enable inflation
if ~isempty(inflate) 
    face = double(face);
    if inflate == 0 ; inflate = 0.02; end
    if inflate ~= 1;
        fprintf('inflating...\n');
        vert = vsmooth([x y z], face, inflate);
        fprintf('done\n');
        x = vert(:,1);
        y = vert(:,2);
        z = vert(:,3);
    end
end


% glass brain
%---------------------------------------------------------
h = patch('faces',face,'vertices',[x(:) y(:) z(:)]);
set(h,'FaceColor',[.4 .4 .4]);
box off;
grid off;
%whitebg(1,'w'); 
camlight('right')
axis tight
set(h,'EdgeColor','none')
material dull
alpha(.2);
lighting phong
set(gca,'visible','off');
set(gcf,'inverthardcopy','off');
hold on;

% functional overlay [ from contrast]

% enable thresholding [slider]
if ~isempty(thr);
    tmap = sparse(tmap);
    tmap(tmap<thr) = 0;
end

% enable blanking of overlay at crit t [trs]
if ~isempty(trs);
    cmap = jet();
    tm   = max(tmap);
    
    % Make values 0-5 gray:
    sz  = length(cmap);
    sz  = round(sz*(.5));    
    mid = sz;
    sz  = ceil(sz/tm*trs);
    
    zs  = zeros(size(cmap(mid-sz:mid+sz,:)));
    cmap(mid-sz:mid+sz,:)=zs+.4;
    colormap(cmap);
end


trisurf(face,x,y,z,tmap); 
alpha(trans)

%set(h,'EdgeColor','interp')
%set(h,'FaceVertexCData',tmap);
shading interp
camlight headlight

%scatter3(x,y,z,[],st,'filled');
%alpha(.3)


% % find vertices corresponding to provided MNIs
% %---------------------------------------------------------
% XYZ   = o;
% inv   = D.inv{end};
% rad   = 1;   
% 
% Ns    = size(XYZ, 1); % n points to plot
% svert = {};
% for i = 1:Ns
%     dist = sqrt(sum([vert(:,1) - XYZ(i,1), ...
%                      vert(:,2) - XYZ(i,2), ...
%                      vert(:,3) - XYZ(i,3)].^2, 2));
%     if rad > 0
%         for j = 1:rad
%             [junk,svert{i,j}] = min(dist);
%             dist(svert{i,j}) = NaN;
%         end
%     else
%         [junk,svert{i}] = min(dist);
%         XYZ(i, :) = vert(svert{i}, :);
%     end
% end
% 
% % add selected points to plot
% %---------------------------------------------------------
% for i = 1:Ns
%     for j = 1:size(svert,2)
%         scatter3(vert(svert{i,j},1),...
%                  vert(svert{i,j},2),...
%                  vert(svert{i,j},3),...
%                  s,CL,'filled');
%                     alpha(.2)
%     end
% end
% 
% 



% MNI=[-46 20 8;
%  -61 -32 8;
%  -42 -14 7; %-22
%   46 20 8;
%   59 -25 8;
%   46 -14 8];






