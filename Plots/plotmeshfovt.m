function plotmeshfovt(D,ol,tval,trans,sname,angs,t)
%
% Like plotmeshfo but generates a video if ol is a matrix
%
% D is a meeg object with [template] mesh
% ol is matrix of overlay by time
% tval is colorbar value at which to blank
% trans is alpha value
% sname is save name.avi
% angs is view(angs(1),angs(2))
% t is t(1):t(2) times, displayed as title
%
% AS

global inflate


checkexist(sname);

fc  = .1; % face col [.9]
frt = 20;
D   = spm_eeg_load(D);

% fix unspecified params
if  isempty(inflate);inflate = [];end
try trans; catch; trans = .4;     end



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

% enable blanking of overlay at crit t [tval]
if ~isempty(tval);
    if any(ol(:)<0)
        mapping = max(abs(mean(ol,1)));
        cmap = jet(length(ol));
        tm   = max(mean(ol,1));
        
        % Make values 0-5 gray:
        sz  = length(cmap);
        sz  = round(sz*(.5));
        mid = sz;
        sz  = ceil(sz/tm*tval);
        
        try    zs  = zeros(size(cmap(mid-sz:mid+sz,:)));
        catch  sz  = sz/length(cmap)*100;
            zs  = zeros(size(cmap(mid-sz:mid+sz,:)));
        end
        cmap(mid-sz:mid+sz,:)=zs+.4;
        colormap(cmap);
    
    else
        % if is abs(y) then thresh lower colbar by n%
        mapping = max(abs(mean(ol,1)));
        %cmap = jet(length(ol));
        cmap = jet(length(ol));
        tm   = max(mean(ol,1));
        
        % Make values 0-5 gray:
        sz  = length(cmap);
        pos = round( tval*(sz/100));
        cmap(1:pos,:) = .4;
        colormap(cmap)
    end
    
    
else
   mapping = max(abs(mean(ol,1)));
 
end

% increase frames by adding averages? [lin interp]
% nfram = 500;
% d     = [];
% 
% while length(a) < nfram
%     mergim   = squeeze(inner3d(a));
%     d        = [];
%     for i    = 1:(length(a)-1)
%         b    = [mergim(i+1,:,:,:);mergim(i,:,:,:)];
%         c    = (squeeze(mean(b,1)));
%         d    = [d [a(i) c]];
%     end
%     a = d;
%     fprintf('Num images %d\n',length(a));
% end



% VIDEO overlay
%-----------------------------------------------------------------
writerObj = VideoWriter([sname '.avi']);
writerObj.FrameRate = size(ol,1)/frt;
writerObj.Quality   = 100;
open(writerObj);
fprintf('\n');

tvec = linspace(t(1),t(2),size(ol,1));
fh   = gcf;

for k = 1:size(ol,1)
    strg = sprintf('Generating frame %d of %d',k,size(ol,1));
    if k > 1; fprintf(repmat('\b',size(strg))); end
    fprintf(strg);

    % glass brain
    cla();
    h = patch('faces',face,'vertices',[x(:) y(:) z(:)]);
    set(h,'FaceColor',[fc fc fc]);
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
    
    mesh    = D.inv{1}.forward(1).mesh;
    ol(k,:) = spm_mesh_smooth(export(gifti(mesh)), (ol(k,:))', 8);

    trisurf(face,x,y,z,ol(k,:)); 
    alpha(trans)
    set(h,'EdgeColor','none')
    %set(h,'FaceVertexCData',ol(k,:)');
    shading interp
    %camlight headlight
    lighting none
    %brighten(.2);
    view(angs(1),angs(2));
    
    title(num2str(tvec(k)),'fontsize',20);
    set(findall(gca, 'type', 'text'), 'visible', 'on');
    
    hold off;

    %caxis([-mapping mapping]);

   frame = getframe(fh);
   writeVideo(writerObj,frame);
end

close(writerObj);

end

function checkexist(name)

if any(exist([name '.avi']));
    disp('File exists! - return to overwrite');pause
end

end
