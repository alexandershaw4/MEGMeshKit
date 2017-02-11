function y = mesh_pca1(D,cond)
% mesh-based pca [with plots] for spm source reconstructed m/eegs
%
% AS

doplot = 0;       % plot 
doimg  = 0;       % save images [g/nifti]
invi   = 1;       % inversion index (see D.val)
woi    = [0 .35]; % window of interest (secs)
foi    = [];      % opt. foi [but reduces rank]
k      = 8;       % smoothing kern size
neig   = 20;      % num eigens 

coni   = strmatch(cond,D.inv{1}.inverse.trials);

%-----------------------------------------------
inv  = D.inv{invi};
vert = inv.forward(end).mesh.vert;
face = inv.forward(end).mesh.face;
pst  = inv.inverse.pst;
T    = inv.inverse.T;
Ja   = inv.inverse.J;
J    = Ja{coni};       
time = D.inv{1}.inverse.pst;
mesh = inv.forward(end).mesh;

% get temp projr for woi
w    = 1; 
fwhm = max(diff(woi(w,:)),8);
t    = exp(-4*log(2)*(pst(:) - mean(woi(w,:))).^2/(fwhm^2));
t    = t/sum(t);
W    = t(:);

if ~isempty(foi)
    wt = 2*pi*pst(:)/1000;
    W  = [];
    for f = foi(1):foi(end)
        W = [W sin(f*wt) cos(f*wt)];
    end
    W  = diag(t)*W;
    W  = spm_svd(W,1);
else
    W  = t(:);
end

TW   = T'*W;
TTW  = T*TW;

% full pst:
if isempty(foi); wtime = J*T';
else             wtime = J*(TW*TTW');
end

wi    = [findthenearest(time,woi(1)*1000):findthenearest(time,woi(2)*1000)];
TOI   = abs(wtime(:,wi));
y     = reduce_eig_mesh(mesh,TOI,neig,k);


if doimg
        writeims(D,mean(y,2),['comp_' cond]);
end
