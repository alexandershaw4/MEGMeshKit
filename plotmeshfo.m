function MNI = plotmeshfo(D,ol,tval,trans,varargin)

global inflate
trs = tval;

% blob coords
try;varargin{1};
      o = varargin{1};
catch o = [];
end; 

try trans;catch;trans = .4; end
D = spm_eeg_load(D);

% fix unspecified params
if isempty(inflate);inflate = [];end

% verts & faces for brain
vert  = D.inv{end}.forward(end).mesh.vert;
x     = vert(:,1);
y     = vert(:,2);
z     = vert(:,3);
face  = D.inv{end}.forward(end).mesh.face;


% vert = D.inv{end}.mesh.tess_mni.vert;
% face = D.inv{end}.mesh.tess_mni.face;
% x     = vert(:,1);
% y     = vert(:,2);
% z     = vert(:,3);

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
h = patch('faces',face,'vertices',[x(:) y(:) z(:)]);
set(h,'FaceColor',[.9 .9 .9]);
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


% enable blanking of overlay at crit t [trs]
if ~isempty(trs);
    mapping = max(abs(ol));
    cmap = jet(length(ol));
    tm   = max(ol);
    
    % Make values 0-5 gray:
    sz  = length(cmap);
    sz  = round(sz*(.5));    
    mid = sz;
    sz  = ceil(sz/tm*trs);
    
    try    zs  = zeros(size(cmap(mid-sz:mid+sz,:)));
    catch  sz  = sz/length(cmap)*100;
           zs  = zeros(size(cmap(mid-sz:mid+sz,:)));
    end
    cmap(mid-sz:mid+sz,:)=zs+.4;
    colormap(cmap);
else
   mapping = max(abs(ol));
 
end


%overlay
trisurf(face,x,y,z,ol); 
alpha(trans)
set(h,'EdgeColor','interp')
set(h,'FaceVertexCData',ol');
shading interp
%camlight headlight
lighting none
brighten(.2);

caxis([-mapping mapping]);



if ~isempty(o) % MNIs
    
% find vertices corresponding to provided MNIs
%---------------------------------------------------------
XYZ   = o;
inv   = D.inv{end};
rad   = 1;   
s     = 20;
CL    = 'b';

Ns    = size(XYZ, 1); % n points to plot
svert = {};

% convert to MNI
nvert = D.inv{end}.mesh.tess_mni.vert;


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

for i = 1:length(svert)
    MNI(i,:) = nvert(svert{i},:);
end

% add selected points to plot
%---------------------------------------------------------
for i = 1:Ns
    for j = 1:size(svert,2)
        scatter3(vert(svert{i,j},1),...
                 vert(svert{i,j},2),...
                 vert(svert{i,j},3),...
                 s,CL,'filled');
                    alpha(trans)
    end
end

else MNI = [];
    
end
