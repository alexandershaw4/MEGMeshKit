function [x,y,z] = get_xyz(n)

try n ; catch n = 1; end

%state = uisuspend(gcf);
%set(gcf,'Pointer','fullcross')

for i=1:n
    
   % if ~waitforbuttonpress
        
%         clickedPt = get(gca,'currentpoint');
%         %XYZ = XYZ(1,:);
%         %XYZ=~(XYZi(1,:)-XYZi(2,:)).*XYZi(1,:);
%         
%         VMtx = view(gca);
%         point2d = VMtx * [clickedPt(1,:) 1]';
%         XYZ = point2d(1:3)';
        
    dcm_obj = datacursormode(gcf);
    waitforbuttonpress

    set(dcm_obj,'DisplayStyle','datatip',...
    'SnapToDataVertex','off','Enable','on');
    drawnow
    c_info = getCursorInfo(dcm_obj);
    
    vert        = handles.vert(handles.Is, :);
    coord       = get(handles.sources_axes, 'CurrentPoint');
    dist        = sum((vert - repmat(coord(1, :), size(vert, 1), 1)).^2, 2);
    [junk, ind] = min(dist);
    coord       = vert(ind, :);
    
    
    trg = get(c_info.Target);
    %Position: [-0.0221 0.0749 0.0765]
    
    
    
    % record click at peak:
    x(i,:) = c_info.Position(1);
    y(i,:) = c_info.Position(2);
    z(i,:) = c_info.Position(3);

%         point = get(gca, 'CurrentPoint'); % mouse click position
%         camPos = get(gca, 'CameraPosition'); % camera position
%         camTgt = get(gca, 'CameraTarget'); % where the camera is pointing to
%         camDir = camPos - camTgt; % camera direction
%         camUpVect = get(gca, 'CameraUpVector'); % camera 'up' vector
%         zAxis = camDir/norm(camDir);
%         upAxis = camUpVect/norm(camUpVect);
%         xAxis = cross(upAxis, zAxis);
%         yAxis = cross(zAxis, xAxis);
%         
%         rot = [xAxis; yAxis; zAxis]; % view rotation
%         rotatedPointFront = rot * point' ;
        
        
    %end
%     x(i,:) = rotatedPointFront(1);
%     y(i,:) = rotatedPointFront(2);
%     z(i,:) = rotatedPointFront(3);
end

%uirestore(state);
set(gcf,'Pointer','arrow')
refresh
