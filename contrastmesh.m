function tmap = contrastmesh(DD,t1,t2,T,F,val)
% constrast conditions from meshs using t-stat
%
% DD = SPM12 MEEG objects [source localised]
% t1 = trial type 1 (to contrast with..)
% t2 = trial type 2
% T  = window of interest 
%
% eg.
% tmap = contrastmesh(f,'Neutral','Neutral',toi,foi); 
% where toi = [-1 0;0.3]
%
% All = contrastmesh(f,{'Neutral','Happy','Angry'},'FT',[0 .3],foi);
% for multiple triggers vs 1
%
% Note that the actual condition names could be Neutral_1,...Neutral_n 
% however the function does partial string matching and includes all
% [partial] matches
% 
% If input array 'DD' is 2xN an independent t-test is performed between
% group DD(1,:) and DD(2,:)
%
% - conditions can be in diff order in diff datasets [read per person]
% - evoked or induced
% AS2016

if nargin < 6; val = []; end
if nargin < 5; foi = []; else foi = F; end
if nargin < 4; woi = []; else woi = T; end
if ~isobject(DD{1}); DD = loadarrayspm(DD); end

type = 'evoked';

D = DD{1};
v = @spm_vec;

if isempty(woi); woi = [D.time(1) D.time(end)]; end
js = 1:size(D.inv{end}.inverse.M,1);

if size(DD,1) > 1 && size(DD,2) > 1
     DD2 = DD(2,:);
     DD  = DD(1,:);
     ttt = @ttest2;
else ttt = @ttest;
end

% Get projections
if all(size(woi,1)>1 && size(woi,2)>1); 
     fprintf('Using different time windows for contrasts\n');
     SepTime = 1; 
     woi2    = woi(2,:);
     woi     = woi(1,:);
else SepTime = 0; 
end

% Trials / conditions 1 [@time1, Group0]
%------------------------------------------------------------------------
DD = DD(~cellfun(@isempty,DD));
for s = 1:length(DD)
    if s == 1; fprintf('Fetching mesh projections; please wait...\n'); end
    D = DD{s};
    switch type
        case 'evoked';  out = rebuild(D,woi,'evoked',foi,val);
        case 'induced'; out = rebuild(D,woi,'induced',foi,val);
    end
    JW(s,:) = out.JW;
end

% TIME 2 [IF REQ]
%------------------------------------------------------------------------
if SepTime
    for s = 1:length(DD)
        if s == 1; fprintf('Fetching mesh projections; please wait...\n'); end
        D = DD{s};
        switch type
            case 'evoked';  out = rebuild(D,woi2,'evoked',foi);
            case 'induced'; out = rebuild(D,woi2,'induced',foi);
        end
        JW2(s,:) = out.JW;
    end
end

% GROUP 2 [IF REQ]
%------------------------------------------------------------------------
if exist('DD2','var')
    for s = 1:length(DD2)
        if s == 1; fprintf('Fetching mesh projections [group 2]; please wait...\n'); end
        D2 = DD2{s};
        switch type
            case 'evoked';  out = rebuild(D2,woi,'evoked',foi);
            case 'induced'; out = rebuild(D2,woi,'induced',foi);
        end
        JW2(s,:) = out.JW;
    end
end


try t1 = eval(t1); end
for s  = 1:length(DD)
    D  = DD{s};
    
    if ischar(t1)
        L        = D.condlist;
        this1{s} = strmatch(t1,L);
    elseif iscell(t1)
        L        = D.condlist;
        for cond = 1:length(t1)
            if cond == 1
                 this1{s} = strmatch(t1{cond},L);
            else this1{s} = [v(this1{:}); v(strmatch(t1{cond},L))]; 
            end
        end
        this1{s} = unique(spm_vec(this1{s}));
        
    elseif isnumeric(t1)
        this1 = t1;
    end
end
try this1 = squeeze(inner(this1));
catch
    try
    this1 = squeeze(inner3d(this1));
    catch
    end
end
fprintf('Found %d trial types for condition 1\n',size(this1,2));


try t2 = eval(t2); end
if exist('DD2','var')
    DD = DD2;
end
for s  = 1:length(DD)
    % if there is a gorup 2, switch to them
    % othwerwise just find the other conditions in same subs
%     if exist('DD2','var');
%          D  = DD2{s};
%     else D  = DD {s};
%     end
    
    if ischar(t2)
        L        = D.condlist;
        this2{s} = strmatch(t2,L);
    elseif iscell(t2);
        L        = D.condlist;
        for cond = 1:length(t2)
            if cond == 1
                % collect indices
                 this2{s} = strmatch(t2{cond},L);
            else this2{s} = [v(this2{:}); v(strmatch(t2{cond},L))]; 
            end
        end
        this2{s} = unique(spm_vec(this2{s}));
        
    elseif isnumeric(t2)
        this2 = t2;
        L = D.condlist;
    end
end
try this2 = squeeze(inner(this2));
catch
    try
    this2 = squeeze(inner3d(this2));
    catch
    end
end
fprintf('Found %d trial types for condition 2\n',size(this2,2));



% Get relevant projections
if ~isvector(this1)
    J1 = JW(this1);
else
    ATR = zeros(length(JW),length(L));
    %ATR(:,this1) = 1;
    for i = 1:size(ATR,1)
        J1(i,1) = JW(i,this1(i));
        %J1(i,1) = JW(i,find(ATR(i,:)));
    end  
end

if SepTime || exist('DD2','var')
    % if different wois or groups
    if ~isvector(this1)
        %J1 = JW(this1);
         J2 = JW2(this2);
    else
        ATR = zeros(length(JW),length(L));
        %ATR(:,this2) = 1;
        for i = 1:size(ATR,1)
            J2(i,1) = JW2(i,this2(i));
            %J2(i,1) = JW2(i,find(ATR(i,:)));
        end  
    end
 
     %J2 = JW2(this2);
else
    if ~isvector(this2)
        J2 = JW(this2);
    else
        ATR = zeros(length(JW),length(L));
        %ATR(:,this2) = 1;
        for i = 1:size(ATR,1)
            J2(i,1) = JW(i,this2(i));
            %J2(i,1) = JW(i,find(ATR(i,:)));
        end  
    end        
end

J1 = squeeze(inner(J1));
J2 = squeeze(inner(J2));

% Average this trial types
if ndims(J1) > 2 
    J1 = squeeze(mean(J1,2));
    J2 = squeeze(mean(J2,2));
end

try
    J1 = smoother(J1); 
    J2 = smoother(J2); 
end

J1 = clust((J1),DD{1}.inv{1}.forward.mesh);
J2 = clust((J2),DD{1}.inv{1}.forward.mesh);



J1 = abs(J1);
J2 = abs(J2);

% t-tests
for s = 1:size(J1,2)
    if s > 1; fprintf(repmat('\b',size(str))); end
    str = sprintf('Computing %d of %d',s,size(J1,2));
    fprintf(str);
    
    [H(s,:),P(s,:),CI,ST] = ttt(J1(:,s),J2(:,s));
    Tst(:,s) = ST.tstat;
end
fprintf('\n');

tmap = Tst;

% do plot
% tmap = Tst;
% figure,
% vert = D.inv{end}.forward(end).mesh.vert;
% face = D.inv{end}.forward(end).mesh.face;
% x = vert(:,1);
% y = vert(:,2);
% z = vert(:,3);
% %b = HighResMeanFilt(tmap,1,2);
% 
% Brain = trisurf(face,x,y,z,tmap,'EdgeColor','none');
% alpha(.7); grid off;set(gca,'visible','off');
% colorbar
end

function o = clust(in,xyz);
fprintf('Spatial smoothing (robust averaging)\n');
A = spm_mesh_adjacency(xyz.face);
N = spm_mesh_neighbours(A);
N = moreN(N);

mn = mean(in,1);
L  = spm_mesh_get_lm(xyz.face,mn');

% sort local maxima
[thev,thein]=sort(mean(in(:,L),1),'descend');
NBLOB = 12;
L = L(thein(1:NBLOB));

% N is a function: nearest neighbours of peak pos L(1) == N(L(1),:)


% Align peaks
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
for s = 1:size(in,1)
    for i = 1:length(N)
        out(s,i) = spm_robust_average( [out(s,i) out(s,N(i,:)) ]');
    end
end

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
    catch NEW(v,:) = ( NEW(v-1,:)*0 ) + N(v,i);
          %NEW(v,:) = N(v,:)  
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