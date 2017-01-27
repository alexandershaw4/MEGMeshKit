% mesh-based ica [with plots] for spm source reconstructed m/eegs
%
% AS

doplot = 1;       % plot 
invi   = 1;       % inversion index (see D.val)
coni   = 1;       % condition index (see D.condlist)
woi    = [0 .35]; % window of interest (secs)

%-----------------------------------------------
inv  = D.inv{invi};
vert = inv.forward(end).mesh.vert;
face = inv.forward(end).mesh.face;
pst  = inv.inverse.pst;
T    = inv.inverse.T;
Ja   = inv.inverse.J;


fwhm = max(diff(woi(1,:)),8);
t    = exp(-4*log(2)*(pst(:) - mean(woi(1,:))).^2/(fwhm^2));
t    = t/sum(t);
W    = t(:);

TW   = T'*W;     % temporal projector
J    = Ja{coni}; % condition specific ica

[~,n]     = PEig90(J);               % prop energ > 90%
[C, A, W] = fastica(J','numOfIC',n); % ica

nc = size(C,1);
W  = pinv(W);

for i = 1:nc
    this{i} = ( C(i,:)'*W(:,i)' )';
    proj{i} = this{i}'*TW(:,1);
end

for i = 1:nc
    o{i} = clusterr(proj{i}',inv.forward(end).mesh);
end

if doplot
    nv = 2;
    nh = ceil(nc/nv);

    for i = 1:nc
        subplot(nv,nh,i),plotmeshfo(D,o{i}',[]); 
    end
end

mo = cat(1,o{:});