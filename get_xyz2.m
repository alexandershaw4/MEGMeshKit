function [x,y,z] = get_xyz2(n)
% Capture x y z coordinates from crosshair click
%
%
%
%
% AS

try n ; catch n = 1; end


rotate3d off;
datacursormode on;
waitforbuttonpress

% get the closest vertex
for i = 1
    %vert        = f.inv{end}.forward(end).mesh.vert(f.inv{end}.inverse.Is,:);
    ax = get(gcf,'children');
    ax = ax(n+1);
    dataObjs = get(ax, 'Children');
    vert  = get(dataObjs(1),'Vertices');

    
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

if nargout<3
    x = [x y z];
end
