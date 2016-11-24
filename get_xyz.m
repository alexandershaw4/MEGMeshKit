function [x,y,z] = get_xyz(n)

try n ; catch n = 1; end

%state = uisuspend(gcf);
%set(gcf,'Pointer','fullcross')

for i=1:n
    
    rotate3d off;
    datacursormode on;
    waitforbuttonpress
    
    % get the closest vertex
    coord = get(gca, 'CurrentPoint');
    
    % retain per pos
    nx(i) = coord(1);
    ny(i) = coord(2);
    nz(i) = coord(3);
    
    datacursormode off;
    rotate3d on;

    % record click at peak:
    x(i,:) = c_info.Position(1);
    y(i,:) = c_info.Position(2);
    z(i,:) = c_info.Position(3);

end

%uirestore(state);
set(gcf,'Pointer','arrow')
refresh
