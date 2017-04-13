function MNI = mesh2mni(D,o)
% convert mesh vertex indices to mni
% by looking it up


XYZ   = o;
inv   = D.inv{end};
vert  = inv.forward(end).mesh.vert;
nvert = inv.mesh.tess_mni.vert;

rad   = 1;   
s     = 500;
CL    = 'r';

Ns    = size(XYZ, 1); % n points to plot
svert = {};

% convert to MNI

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