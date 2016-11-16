function [x,y,z] = ind2mni(D,ox,oy,oz)

if ~isobject(D); D = spm_eeg_load(D); end


% Get mesh verts
vert  = D.inv{end}.forward(end).mesh.vert;
rad   = 1;   
XYZ   = [ox oy oz];

Ns    = size(XYZ, 1); % n points to plot
svert = {};

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

% closest vertices:
for i = 1:Ns
    x(i,:) = (vert(svert{i},1));
    y(i,:) = (vert(svert{i},2));
    z(i,:) = (vert(svert{i},3));
end