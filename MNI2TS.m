function [dD,XYZ] = MNI2TS(D,MNI,rad)
% From a source localised MEEG object, extract timeseries from MNI
% coordinates
%
% D   = SPM MEEG object
% MNI = [nx3] MNI coordinates
% rad = radius to search around MNI coords
%
% dD  = 3D double dD(source x samples x conditions) where nsource is nMNI
% XYZ = the [MNI] location of the vertices found within radius
%
% Based partially on spm_eeg_inv_extract
% Notes: - trials will be averaged to account for trial number differences
%        - designed for multimodal data [MAGS+GRADS+EEG]
% AS2016


% Input options
try, XYZ = MNI;   catch,  return;     end
try, rad      ;   catch,  rad   = 5;  end


% Find sources near MNI coordinates
inv   = D.inv{end};
vert  = inv.mesh.tess_mni.vert(inv.inverse.Is, :);  % vertices
Ns    = size(XYZ, 1);                               % number of sources
svert = {};
for i = 1:Ns
    dist = sqrt(sum([vert(:,1) - XYZ(i,1), ...
                     vert(:,2) - XYZ(i,2), ...
                     vert(:,3) - XYZ(i,3)].^2, 2));
    if rad > 0
        svert{i} = find(dist < rad);
    else
        [junk,svert{i}] = min(dist);
        XYZ(i, :) = vert(svert{i}, :);
    end
end

% make sure we found proximal verices
Iy = find(cellfun('isempty', svert));
if ~isempty(Iy)
    disp(['No proximal vertices found for source(s) ' num2str(Iy)]);
    dD = [];
    return;
end

% call GetSS for these verts
js = spm_vec(svert);

for i = 1:length(svert)
    fprintf('Found %d vertices within radius for MNI set %d.\n These will be averged.\n',length(svert{i}),i);
    SS = GetSS(D,spm_vec(svert{i}));

    dD(i,:,:) = SS;

end

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
    fprintf('reconstructing trial %d\n',i);
    
    for j = 1:Nt
        for k = 1:length(U)
            Y       = D(Ic{k},It,c(j));
            UY{k,1} = U{k}*Y*scale(k);
        end
        Y = spm_cat(UY);
        
        if j > 1
            MY{i} = MY{i} + M*Y;
        else
            MY{i} = M*Y;
        end
    end
    MY{i} = MY{i} / j;
end

MY = squeeze(inner(MY));
SS = permute(MY,[3 2 1]);
SS =        (mean(SS,1));

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