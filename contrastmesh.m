function tmap = contrastmesh(D,t1,t2,T)
% constrast conditions from mesh vertices using t-stat
%
% D  = SPM12 MEEG object [source localised]
% t1 = trial type 1 (to contrast with..)
% t2 = trial type 2
% T  = window of interest 
%
% AS2016

if nargin < 5; foi = []; else foi = F; end
if nargin < 4; woi = []; else woi = T; end

if isempty(woi); woi = [D.time(1) D.time(end)]; end
js = 1:size(D.inv{end}.inverse.M,1);

for i = 1:length(js)
     tmp = GetSS(D,js(i));
     SS(i,:,:) = full(tmp);
end

% extract [averaged] trials of interest
T1   = squeeze(SS(:,:,t1));
T2   = squeeze(SS(:,:,t2));

% woi
time = D.time;
ti   = [findthenearest(woi(1),time):findthenearest(woi(2),time)];
T1   = T1(:,ti);
T2   = T2(:,ti);

% t-tests
for s = 1:size(T1,1)
    [H(s,:),P(s,:),CI,ST] = ttest(T1(s,:),T2(s,:));
    Tst(s,:) = ST.tstat;
end

% do plot
tmap = Tst;
figure,
vert = D.inv{end}.forward(end).mesh.vert;
face = D.inv{end}.forward(end).mesh.face;
x = vert(:,1);
y = vert(:,2);
z = vert(:,3);
%b = HighResMeanFilt(tmap,1,2);

Brain = trisurf(face,x,y,z,tmap,'EdgeColor','none');
alpha(.7); grid off;set(gca,'visible','off');
colorbar
end

function SS = GetSS(D,js)
% Convert from D(sens ,samps,trials)
%         to   D(verts,samps,trials)

inv     = D.inv{end};
%js      = 1:length(inv.inverse.M);
scale   = inv.inverse.scale;
J       = inv.inverse.J;  
U       = inv.inverse.U; % spatial projector 
T       = inv.inverse.T; %
TT      = T*T';
M       = inv.inverse.M(js, :);
Ic      = inv.inverse.Ic;
It      = inv.inverse.It;
Np      = length(It);
trial   = D.condlist;
clabell = {};

for i = 1:numel(trial)
    c       = D.indtrial(trial{i}, 'GOOD');
    clabell = [clabell D.conditions(c)];
    
    Nt    = length(c);
    fprintf('reconstructing trial type %d\n',i);
    
    for j = 1:Nt
        for k = 1:length(U)
            Y       = D(Ic{k},It,c(j));
            UY{k,1} = U{k}*Y*scale(k);
        end
        Y = spm_cat(UY);
        
        if j > 1
            %MY{i} = MY{i} + M*Y;
            %MY{i}(:,c(j)) = mean(M*Y,1);
            MY(:,c(j)) = mean(M*Y,1);
        else
            %MY{i} = mean(M*Y,1);
            %MY{i}(:,c(j)) = mean(M*Y,1);
            MY(:,c(j)) = mean(M*Y,1);
        end
    end
    %MY{i} = MY{i} / j;
end


%MY = squeeze(inner(MY));
%SS = permute(MY,[3 2 1]);
%SS =        (mean(SS,1));

SS = MY;
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