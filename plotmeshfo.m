function plotmeshfo(D,ol,tval,trans)

global inflate
trs = tval;

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
    
    zs  = zeros(size(cmap(mid-sz:mid+sz,:)));
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