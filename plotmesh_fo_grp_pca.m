function [mst,st] = plotmesh_fo_grp_pca(D,t,cfg)
% Plot glass mesh brain with MNI coordinates marked
% from source localised SPM MEEG objects
%
% D is a cell array of subjects' D filenames
% t is trial(s) types [conditions] to include - e.g. 1 or [1 5 6] or 'On'
%
% See mesh_pca1.m for options
%
%
% AS2016

% options
%---------------------------------------------------------
robust_cluster = 0;


if iscell(D) && ~isobject(D{1}); D = loadarrayspm(D);
elseif ~iscell(D); return;
end

warning off 

% verts & faces for brain
%---------------------------------------------------------
vert  = D{1}.inv{end}.forward(end).mesh.vert;
face  = D{1}.inv{end}.forward(end).mesh.face;


% call function for projecting into source space
for s = 1:length(D)
    
    if isnumeric(t)
        for nt = 1:length(t)
            tt{nt,:} = D{s}.inv{end}.inverse.trials{t(nt)};
        end
    else
        tt = t;
    end
    
    clear y
    for nt = 1:length(tt)
        try    y{nt} = mesh_pca1(D{s},tt(nt),cfg);
        catch; y{nt} = mesh_pca1(D{s},tt(nt));
        end
    end
    
    yy     = cat(2,y{:});
    yy     = mean(yy,2);
    
    st(s,:) = yy;
end

% between subject smoothing
st  = spm_conv(st,1.5);


% Experimental clustering [no]
if robust_cluster
    st = clust(st,D{1}.inv{end}.forward(end).mesh);
end

mst  = spm_robust_average(st);


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
fprintf('(e.) Smoothing\n');
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
