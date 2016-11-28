function [x,y,z] = get_xyz(n,f)

try n ; catch n = 1; end

if ischar(f); f = spm_eeg_load(f); end

rotate3d off;
datacursormode on;
waitforbuttonpress

% get the closest vertex
for i = 1:n
    vert        = f.inv{end}.forward(end).mesh.vert(f.inv{end}.inverse.Is,:);
    coord       = get(gca, 'CurrentPoint');
    dist        = sum((vert - repmat(coord(1, :), size(vert, 1), 1)).^2, 2);
    [junk, ind] = min(dist);
    coord       = vert(ind, :);
    
    % retain per pos
    x(i) = coord(1);
    y(i) = coord(2);
    z(i) = coord(3);
end

datacursormode off;
rotate3d on;

