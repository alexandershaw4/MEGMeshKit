function [mst,MNI] = plotmesh_fo_grp(D,o,t,woi,foi,type,CL,trans)
% Plot glass mesh brain with MNI coordinates marked
% from source localised SPM MEEG objects
%
% D is a cell array of subjects' D filenames
% o is a optional MNI coords to project [n x 3]     [optional]
% t is trial(s) types [conditions] to include - e.g. 1 or [1 5 6] or 'On'
%
% woi  is time window of interest                   [optional]
% foi  is freq window of interest                   [optional]
% type is 'evoked' 'induced' or 'trials'            [optional]
% CL   is colour of projected MNI blobs             [optional]
% trans is transparency of overlay                  [optional]
%
% AS2016

% options
%---------------------------------------------------------
robust_cluster = 1;

dS    = 100; % dot size for functional overlay of trial{t}
s     = 500; % patch size for MNI coordinates

try CL;    catch; CL    = 'r';    end  %[.4 .4 .4];%'r'; % colour of MNI patch
try woi;   catch; woi = [-.1 .3]; end  % time window if interest for source data in trial{t}
try foi;   catch; foi = [];       end  % freq window if interest for source data in trial{t}
try type;  catch; type = 'evoked';end  % 'evoked', 'induced' or 'trial'
try trans; catch; trans= .4;      end


if iscell(D) && ~isobject(D{1}); D = loadarrayspm(D);
elseif ~iscell(D); return;
end

warning off %
global subnorm
global thr
global trs
global flt
global npca
global inflate

if thr == 0; thr = []; end

% verts & faces for brain
%---------------------------------------------------------
vert  = D{1}.inv{end}.forward(end).mesh.vert;
x     = vert(:,1);
y     = vert(:,2);
z     = vert(:,3);
face  = D{1}.inv{end}.forward(end).mesh.face;


% fix unspecified parameters
if isempty(woi);   woi   = [0 .3]; end
if isempty(type);  type  = 'evoked'; end
if isempty(trans); trans = .4; end
if isempty(subnorm); subnorm = []; end
if isempty(flt);    flt = []; end
if isempty(npca);   npca = []; end
if isempty(inflate);inflate = [];end


% enable inflation
if ~isempty(inflate) 
    face = double(face);
    if inflate == 0 ; inflate = 0.02; end
    fprintf('inflating...\n');
    vert = vsmooth([x y z], face, inflate);
    fprintf('done\n');
    x = vert(:,1);
    y = vert(:,2);
    z = vert(:,3);
end

% glass brain
%---------------------------------------------------------
h = patch('faces',face,'vertices',[x(:) y(:) z(:)]);
set(h,'FaceColor',[.3 .3 .3]);
box off;
grid off;
%whitebg(1,'w'); 
%camlight('headlight')
%axis tight
set(h,'EdgeColor','none')
material dull
alpha(.2);
%lighting flat %phong
set(gca,'visible','off');
set(gcf,'inverthardcopy','off');
hold on;


% call function for projecting into source space
%---------------------------------------------------------
for SUB = 1:length(D)
    if isempty(foi)
        if SUB > 1; fprintf(repmat('\b',[size(str),1])); end
        str = sprintf('Fetching projections for %d of %d datasets\n',SUB,length(D));
        fprintf(str);
    end

    FO        = rebuild(D{SUB},woi,type,foi);
    if isnumeric(t) && length(t) == 1
        % just get this trial / type
        it        = FO.JW{t};
        st(SUB,:) = it;     
    elseif isnumeric(t)
        % get trials of vector 
        it        = spm_cat({FO.JW{t}});
        st(SUB,:) = mean(it,2);
    elseif iscell(t)
        % get indices of this condition[s] name [spm]
        L        = D{1}.condlist;
        for cond = 1:length(t)
            this = find(strcmp(t{cond},L));
            
            if ~isempty(this) 
                % found condition specified
                fprintf('Averaging projections for %d condition labels\n',length(t));
                it      = FO.JW{this};
                if cond == 1
                    st(SUB,:) = [mean(it,2)*(1/length(t))];
                else
                    st(SUB,:) = [mean(it,2)*(1/length(t))]' + st(SUB,:);
                end
                
            elseif isempty(this)
                % find any partial matches
                this(SUB,:) = strmatch(t{cond},L);
                if ~isempty(this)
                   fprintf('Averaging projections for %d condition labels\n',length(this));

                for k = 1:length(this)
                    it        = FO.JW{this(k)};
                    if cond == 1
                        st(SUB,:) = [mean(it,2)*(1/length(t))];
                    else
                        st(SUB,:) = [mean(it,2)*(1/length(t))]' + st(SUB,:);
                    end
                end
                
                else
                    fprintf('no conditions found!\n');
                end
            end
        end
        
    end
end




if ~iscell(st);
    %st = PEig(full(st'));
    if ~isempty(subnorm)
        fprintf('normalising subjects prior to averaging using %s\n',func2str(subnorm.f));
        %st = TSNorm(st',6,1,1)';
        st = feval(subnorm.f,st',subnorm.type,subnorm.range)';
    elseif ~isempty(flt)
        if any(size(st)==1);
             st  = HighResMeanFilt(diag(st),1,round(flt));
             st  = diag(st)';
        else st  = HighResMeanFilt(st,1,round(flt));
        end
    end
    
    % otherwise
    %mst = mean(st,1);
    %st  = smoother(st);
    st  = spm_conv(st,1.5);
    %mst = spm_robust_average(st);
    
    %st = smoother(st);
    
    if robust_cluster 
        st = clust(st,D{1}.inv{end}.forward(end).mesh);
    end
    
    
    mst  = spm_robust_average(st);
    %mst = spm_conv(st,1.5);
    
    
elseif ~isempty(flt)
    st  = HighResMeanFilt(st,1,round(flt));
    mst = squeeze(cat(2,st{:}));
    %mst = mean(mst,2);
    mst = spm_robust_average(st);
else
    mst = squeeze(cat(2,st{:}));
    %mst = mean(mst,2);
    mst = spm_robust_average(st);
end

% enable pca
if ~isempty(npca)
    mst = PEig(mst',npca)';
end


% enable thresholding [slider]
if ~isempty(thr);
    mst = sparse(mst);
    mst(mst<thr) = 0;
end


%scatter3(x,y,z,[],mst(:),'filled');alpha(.3);
%h = patch('faces',face,'vertices',[x(:) y(:) z(:)],'facevertexcdata',mst(:));%,...

if ~isempty(trs);
    mapping = max(abs(mst));
    cmap = jet(length(mst));
    tm   = max(mst);
    tm   = mapping;
    
    % Make values 0-5 gray:
    sz  = length(cmap);
    sz  = round(sz*(.5));    
    mid = sz;
    %sz  = ceil(sz/tm*trs);
    %sz  = ceil(tm/sz*trs);
    
    sz = ceil(length(cmap)/(trs*100));
    
    
    zs  = zeros(size(cmap(mid-sz:mid+sz,:)));
    cmap(mid-sz:mid+sz,:)=zs+.4;
    colormap(cmap);
end

trisurf(face,x,y,z,mst); 
alpha(trans)
set(h,'EdgeColor','interp')
set(h,'FaceVertexCData',mst');
shading interp
axis tight
caxis([-mapping mapping]);
%camlight headlight



% enable blanking of overlay
% if ~isempty(trs);
%     cmap = jet();
%     % Make values 0-5 black:
%     %cmap(1:6,:) = zeros(6,3)+.4;
%     sz = length(cmap);
%     sz = round(sz*(.5*trs));
%     cmap(1:sz,:)=zeros(sz,3)+.4;
%     colormap(cmap);
% end

  
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

end

end

function y = smoother(x)
fprintf('Smoothing..\n');
for s = 1:size(x,1)
    np    = round(.2 * size(x,2) );
    %fprintf('(a.) finding potential outliers\n');
    [v,i] = findpeaks(x(s,:),'SortStr','descend','NPeaks',np);
    thr   = findthenearest(mean(v),v);
    
    v = v(thr:end);
    i = i(thr:end);
    
    %v = v(1:thr);
    %i = i(1:thr);
    
    %fprintf('(b.) weighting %d potential outliers\n',length(i));
    
    for k = 1:length(i)
        d       = x(s,:);
        d(i(k)) = 0;
        m       = mean(d);
        
        try        win = [ x(s,i(k)-1) x(s,i(k)+1) ];
        catch try  win = [ m           x(s,i(k)+1) ];
            catch  win = [ x(s,i(k)-1) m           ];
            end
        end
            
        x(s,i(k)) = mean(win);
    end
    
end
y = x;
end

function o = clust(in,xyz);
fprintf('Spatial smoothing...\n');
A = spm_mesh_adjacency(xyz.face); fprintf('(a.) Finding neighbours\n');
N = spm_mesh_neighbours(A);       fprintf('(b.) Creating radius function\n');
N = moreN(N);

mn = mean(in,1);  fprintf('(c.) Finding local maxima\n');
L  = spm_mesh_get_lm(xyz.face,mn');

% sort local maxima
fprintf('(d.) Sorting maxima\n');
[thev,thein]=sort(mean(in(:,L),1),'descend');
NBLOB = 12;
L = L(thein(1:NBLOB));

% N is a function: nearest neighbours of peak pos L(1) == N(L(1),:)


% Align peaks
fprintf('(e.) Aligning [peaks]\n');
for v = 1:length(L)
    LR(v,:,:) = [in(:,L(v)) in(:,N(L(v),:))];
end

out = in;

for s = 1:size(in,1)
    for l = 1:length(L)
        out(s,L(l)) = max(squeeze(LR(l,s,:)));
    end
end

N(find(N==0)) = 1;

% Smooth all vertices
fprintf('(f.) Robust Averaging\n'); nn = 0;
    for i = 1:length(N)
        nn = nn + 1;
        if nn > 1; fprintf(repmat('\b',size(str))); end
        str = sprintf('%d / %d',nn,(length(N)));
        fprintf(str);
        
        out(:,i) = spm_robust_average( [out(:,i) out(:,N(i,:)) ]');
    end
    disp('');
% for s = 1:size(in,1)
%     for i = 1:length(N)
%         nn = nn + 1;
%         if nn > 1; fprintf(repmat('\b',size(str))); end
%         str = sprintf('%d / %d',nn,(size(in,1)*length(N)));
%         fprintf(str);
%         
%         out(s,i) = spm_robust_average( [out(s,i) out(s,N(i,:)) ]');
%     end
% end

o = out;

end

function y = moreN(N)

for v = 1:size(N,1)
    colla = [];
    for i = 1:6
        t = N(v,i);
        
        if t ~= 0            
            colla = [colla N(t,:)];
        end

    end
    try   NEW(v,:) = [N(v,:) colla];
    catch NEW(v,:) = ( NEW(v-1,:)*0 ) + N(v,i); % if we're at an edge
    end

end

y = NEW;
end

function y = smoother(x)
fprintf('Smoothing..\n');
for s = 1:size(x,1)
    np    = round(.6 * size(x,2) );
    [v,i] = findpeaks(x(s,:),'SortStr','descend','NPeaks',np);
    thr   = findthenearest(mean(v),v);
    
    %v = v(thr:end);
    %i = i(thr:end);
    
    v = v(1:thr);
    i = i(1:thr);
    
    %fprintf('smoothing %d peaks\n',length(i));
    
    for k = 1:length(i)
        d       = x(s,:);
        d(i(k)) = 0;
        m       = mean(d);
        
        try        win = [ x(s,i(k)-1) x(s,i(k)+1) ];
        catch try  win = [ m           x(s,i(k)+1) ];
            catch  win = [ x(s,i(k)-1) m           ];
            end
        end
            
        x(s,i(k)) = mean(win);
    end
end
y = x;
end

function y = inner(x)

[S1,S2] = size(x);
[S3,S4] = size(x{1});
y       = zeros(S1,S2,S4,S3);

for i = 1:S1
    for j = 1:S2
        y(i,j,:,:) = x{i,j}';
    end
end

end
    



% MNI=[-46 20 8;
%  -61 -32 8;
%  -42 -14 7; %-22
%   46 20 8;
%   59 -25 8;
%   46 -14 8];
