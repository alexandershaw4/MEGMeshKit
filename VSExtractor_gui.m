function varargout = VSExtractor_gui(varargin)
% An interactive plotting tool for extracting timeseries from a specific
% source location in source localised MEEG SPM datasets.
%
% Works over group data to find the nearest vertices [in individual
% datasets] to the clicked location. Can take multiple clicks / locations
% and can export the individuals extracted data either to the matlab
% workspace or into a new n-channel LFP spm file.
%
% *currently can only do 1 virtual sensor location at a time
% see also Click_to_extract_sensors for script version
%
% AS2016

% Last Modified by GUIDE v2.5 17-Nov-2016 17:29:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @VSExtractor_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @VSExtractor_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before VSExtractor_gui is made visible.
function VSExtractor_gui_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for VSExtractor_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes VSExtractor_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = VSExtractor_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)


% calculate and plot
if isfield(handles,'MNI'); M = handles.MNI; else M = []; end
if isfield(handles,'M'); M = handles.M; else M = []; end
if isfield(handles,'t'); t = handles.t; else t = 1 ; end

if isfield(handles,'woi');  woi  = handles.woi;  else woi  = []; end
if isfield(handles,'foi');  foi  = handles.foi;  else foi  = []; end
if isfield(handles,'type'); type = handles.type; else type = 'evoked'; end
if isfield(handles,'CL');   CL   = handles.CL;   else CL   = 'm'; end
if isfield(handles,'trans');trans= handles.trans;else trans= []; end

f = handles.G;

plotmesh_fo_grp(f,M,t,woi,foi,type,CL,trans);
h=rotate3d;
set(h,'Enable','on');
drawnow

% Save the handles structure.
guidata(hObject,handles); 




function edit1_Callback(hObject, eventdata, handles)

% Get time values if specified
handles.woi = str2num( get(hObject,'String') );

% Save the handles structure.
guidata(hObject,handles); 


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)


% Get trial / condition selection if specified
try 
    if iscell( eval( get(hObject,'String') ) );
       handles.t = ( eval( get(hObject,'String') ) );
    end
catch
       handles.t = str2num( get(hObject,'String') );
end

% Save the handles structure.
guidata(hObject,handles); 


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)

% Load datasets
%-------------------------------
%try 
    [G,Pth] = uigetfile('*.mat','MultiSelect','on');
    
    if ~iscell(G); G = {G}; end
    
    for i = 1:length(G)
        GG{i} = [Pth G{i}];
    end
    
    GG = loadarrayspm(GG);
    
    for i = 1:length(GG)
        handles.G{i,:} = GG{i}; % append paths
    end
    
set(gca,'visible','off');

% Save the handles structure.
guidata(hObject,handles); 

%catch; return;
%end



function edit3_Callback(hObject, eventdata, handles)

% Get MNI coords for projection on surface
handles.MNI = str2num( get(hObject,'String') );

% Save the handles structure.
guidata(hObject,handles); 


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)


% Radius to search when extracting
handles.radius = str2num( get(hObject,'String') );

% Save the handles structure.
guidata(hObject,handles); 


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)


% Number positions to click
handles.npos = str2num(get(hObject,'String'));

% Save the handles structure.
guidata(hObject,handles); 


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)

% Do crosshair capture
%----------------------------------------------
try npos = handles.npos; catch npos = 1; end
for k = 1:npos
    [rx(k,:),ry(k,:),rz(k,:)] = get_xyz();  % 
end
G = handles.G;


% Do conversion to verts per person
%----------------------------------------------
for j = 1:length(G)
    for i = 1:npos
        
        [nx(j,i,:) ny(j,i,:) nz(j,i,:)]  = ind2mni(G{j},rx(i),ry(i),rz(i));      % find nearest vertex to selected point for each subject
        
    end
end

handles.capt.x = nx;
handles.capt.y = ny;
handles.capt.z = nz;

handles.orig.x = rx;
handles.orig.y = ry;
handles.orig.z = rz;


% Save the handles structure.
guidata(hObject,handles); 


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)


% extract this sensor
f  = handles.G;
rx = handles.orig.x;
ry = handles.orig.y;
rz = handles.orig.z;
Q  = @squeeze;

try r  = handles.radius; catch r = 4; end

for i = 1:length(f)
    
    [dD{i},C,T,XYZ{i},FT] = MNI2TS(f{i},[Q(rx) Q(ry) Q(rz)],r); % extract full dataset [samples by trials] for each identified vertex

end

assignin('base','D',dD);
assignin('base','Conds',C);
assignin('base','Time',T);
assignin('base','XYZ',XYZ);



function edit6_Callback(hObject, eventdata, handles)

% Get frequencies of interest
handles.foi = str2num( get(hObject,'String') );

% Save the handles structure.
guidata(hObject,handles); 


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)

% Project coordinates onto brain
%M(:,1) = handles.orig.x;
%M(:,2) = handles.orig.y;
%M(:,3) = handles.orig.z;
X = handles.capt.x; X = (X(1,:));
Y = handles.capt.y; Y = (Y(1,:));
Z = handles.capt.z; Z = (Z(1,:));

handles.M  = [X' Y' Z'];
handles.CL = 'b';

% Save the handles structure.
guidata(hObject,handles); 

pushbutton1_Callback(hObject, eventdata, handles)


function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

X = handles.capt.x; 
Y = handles.capt.y;
Z = handles.capt.z;

assignin('base','x',X);
assignin('base','y',Y);
assignin('base','z',Z);


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

f = handles.G;
Q = @squeeze;
n = handles.npos;
%M(:,:,1) = handles.orig.x;
%M(:,:,2) = handles.orig.y;
%M(:,:,3) = handles.orig.z;

M(:,:,1) = handles.capt.x;
M(:,:,2) = handles.capt.y;
M(:,:,3) = handles.captc.z;


for i = 1:n
    L{i} = ['Sensor_' num2str(i)];
end

try r = handles.r; catch r = 5; end

fprintf('Please wait...\n');
for i = 1:length(f)
    D                     = f{i};
    D.inv{end}.source.label = L';
    
    D.inv{end}.source.XYZ   = Q(M(i,:,:));
    D.inv{end}.source.rad   = r;
    D.inv{end}.source.fname = ['Extract_' D.fname];
    D.inv{end}.source.type  = 'trials';

    [Ds, D] = spm_eeg_inv_extract(D);
end
fprintf('Done\n');


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)

% type in MNI coords
opts.Interpreter =  'tex';
tMNI = inputdlg('Enter MNIs - e.g. [7 10 20]','MNIs',6,{'[0 47 12]'},opts);

for i = 1:size(tMNI{1},1)
    M(i,:) = str2num(tMNI{:}(i,:));
end

handles.M = M;

% Save the handles structure.
guidata(hObject,handles); 


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Clear! start over
try handles = rmfield(handles,'M');end
try handles = rmfield(handles,'G');end
try handles = rmfield(handles,'MNI');end
try handles = rmfield(handles,'radius'); end
try handles = rmfield(handles,'foi');end
try handles = rmfield(handles,'woi');end
try handles = rmfield(handles,'CL');end
try handles = rmfield(handles,'orig'); end

clc;
cla();
% Save the handles structure.
guidata(hObject,handles); 


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Set transparency slider
v = get(hObject,'Value');% * ( get(hObject,'Min') / get(hObject,'Max') );

handles.trans = v;

cla();
pushbutton1_Callback(hObject, eventdata, handles);

% Save the handles structure.
guidata(hObject,handles); 




% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla();
