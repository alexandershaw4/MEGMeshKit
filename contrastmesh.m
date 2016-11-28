function tmap = contrastmesh(DD,t1,t2,T,F)
% constrast conditions from mesh vertices using t-stat
%
% DD = SPM12 MEEG objects [source localised]
% t1 = trial type 1 (to contrast with..)
% t2 = trial type 2
% T  = window of interest 
%
% AS2016

if nargin < 5; foi = []; else foi = F; end
if nargin < 4; woi = []; else woi = T; end
if ~isobject(DD{1}); DD = loadarrayspm(DD); end

type = 'induced';

D = DD{1};

if isempty(woi); woi = [D.time(1) D.time(end)]; end
js = 1:size(D.inv{end}.inverse.M,1);

% Get projections
for s = 1:length(DD)
    D = DD{s};
    switch type
        case 'evoked'
            out = rebuild(D,woi,'evoked',foi);
        case 'induced';
            out = rebuild(D,woi,'induced',foi);
    end
    JW(s,:) = out.JW;
end

if ischar(t1);
    L     = D.condlist;
    this1 = strmatch(t1,L);
elseif isnumeric(t1)
    this1 = t1;
end
fprintf('Found %d trial types for condition 1\n',length(this1));

if ischar(t2)
    L     = D.condlist;
    this2 = strmatch(t2,L);
elseif isnumeric(t2)
    this2 = t2;
end    
fprintf('Found %d trial types for condition 2\n',length(this2));

% Get relevant projections
J1 = JW(:,this1);
J2 = JW(:,this2);
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
    
    [H(s,:),P(s,:),CI,ST] = ttest(J1(:,s),J2(:,s));
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