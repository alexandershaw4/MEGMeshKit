% mesh-based ica [with plots] for spm source reconstructed m/eegs
%
% AS

doplot = 1;       % plot 
doimg  = 1;       % save images [g/nifti]
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
J    = Ja{coni};       
time = D.inv{1}.inverse.pst;

% get temp projr for woi
w    = 1; 
fwhm = max(diff(woi(w,:)),8);
t    = exp(-4*log(2)*(pst(:) - mean(woi(w,:))).^2/(fwhm^2));
t    = t/sum(t);
W    = t(:);
TW   = T'*W;
TTW  = T*TW;

% full pst:
wtime = J*T';
wi    = [findthenearest(time,woi(1)*1000):findthenearest(time,woi(2)*1000)];
TOI   = wtime(:,wi);

[~,n]     = PEig90(TOI);              % prop energ > 90%
[C, A, W] = fastica(TOI,'numOfIC',n); % ica

nc = size(C,1);
W  = pinv(W);

for i = 1:nc
    this{i} = ( C(i,:)'*W(:,i)' )';   % spatio-temporal component projections
end

for i = 1:nc
    o{i} = clusterr(this{i}',inv.forward(end).mesh); % robust average
end

for i = 1:nc
    o{i} = o{i}'*TTW;                 % temporal projector TTW
end

if doplot
    nv = 2;
    nh = ceil(nc/nv);

    for i = 1:nc
        subplot(nv,nh,i),plotmeshfo(D,o{i}',[]); 
    end
end

mo = cat(1,o{:});

if doimg
    for i = 1:nc
        writeims(D,o{i},['comp_' num2str(i)]);
    end
end



