function tmap = contrastmesh(DD,t1,t2,T,F)
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
        case 'evoked';  out = rebuild(D,woi,'evoked',foi);
        case 'induced'; out = rebuild(D,woi,'induced',foi);
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
    ATR(:,this1) = 1;
    for i = 1:size(ATR,1)
        J1(i,1) = JW(i,find(ATR(i,:)));
    end  
end

if SepTime || exist('DD2','var')
    % if different wois or groups
     J2 = JW2(this2);
else J2 = JW(this2);
end

J1 = squeeze(inner(J1));
J2 = squeeze(inner(J2));

% Average this trial types
if ndims(J1) > 2
    J1 = squeeze(mean(J1,2));
    J2 = squeeze(mean(J2,2));
end

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