function [x,y,z] = ind2mni(D,ox,oy,oz)
% obsolete
%
%
%

if ~isobject(D); D = spm_eeg_load(D); end


% Get mesh verts
XYZ    = D.inv{end}.forward(end).mesh.vert;
rad    = 1;   
vert   = [ox oy oz];

Ns    = size(vert, 1); % n points to plot
svert = {};

for i = 1:Ns
    dist = sqrt(sum([vert(i,1) - XYZ(:,1), ...
                     vert(i,2) - XYZ(:,2), ...
                     vert(i,3) - XYZ(:,3)].^2, 2));                 
%     if rad > 1
%         for j = 1:rad
%             [junk,svert{i,j}] = min(dist);
%             dist(svert{i,j}) = NaN;
%         end
%     else
        [junk,svert{i}] = min(dist);
        xyz(i, :) = XYZ(svert{i}, :);
%     end
end

x = xyz(:,1);
y = xyz(:,2);
z = xyz(:,3);


% closest vertices:
% for i = 1:Ns
%     x(i,:) = (vert(svert{i},1));
%     y(i,:) = (vert(svert{i},2));
%     z(i,:) = (vert(svert{i},3));
% end