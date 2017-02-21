function varargout = contrast_plot_gui(varargin)
% Plot window for contrasts run in contrast_gui
% See help VSExtractor
%
%
%
% AS2016

% Edit the above text to modify the response to help contrast_plot_gui

% Last Modified by GUIDE v2.5 30-Nov-2016 17:20:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @contrast_plot_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @contrast_plot_gui_OutputFcn, ...
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


% --- Executes just before contrast_plot_gui is made visible.
function contrast_plot_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to contrast_plot_gui (see VARARGIN)

% Choose default command line output for contrast_plot_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes contrast_plot_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);
drawnow
Redoplot(hObject,handles)

 


% --- Outputs from this function are returned to the command line.
function varargout = contrast_plot_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function Redoplot(hObject,handles)
% pick up data from contrast gui when called
%-----------------------------------------------------
global ConPlot
global thr
global trs

nc   = ConPlot.nc;
tmap = ConPlot.tmap;
trans = ConPlot.trans;
G     = ConPlot.G;
cla;
drawnow
    
set(0, 'CurrentFigure', handles.figure1)

for k = 1:nc
    % t ->
    subplot(2,nc,k);
    plotmesh_fo_tmap(G,tmap{k},thr,trs,trans);
    % t <-
    subplot(2,nc,k+nc);
    plotmesh_fo_tmap(G,tmap{k}*-1,thr,trs,trans);
end
ac = allchild(handles.figure1);
linksubplots(ac);
rotate3d on;

% Save the handles structure.
guidata(hObject,handles); 


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% slider 1 - alpha
global ConPlot
v = get(hObject,'Value');% * ( get(hObject,'Min') / get(hObject,'Max') );
ConPlot.trans = v;

cla();drawnow
Redoplot(hObject,handles)
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


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% slider 2: threshold for t(+/- n) = 0
v = get(hObject,'Value');% * ( get(hObject,'Min') / get(hObject,'Max') );

global thr 
thr = v;

cla();drawnow
Redoplot(hObject,handles)
% Save the handles structure.
guidata(hObject,handles); 


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% slider 3:
% blankify
level = get(hObject,'Value'); 

global trs
trs = level;

cla();drawnow
Redoplot(hObject,handles)
% Save the handles structure.
guidata(hObject,handles); 


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% slider 4:
global inflate
inflate = get(hObject,'Value'); 
cla();drawnow
Redoplot(hObject,handles)
% Save the handles structure.
guidata(hObject,handles); 


% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% save image
% opts.Interpreter =  'tex';
% fname = inputdlg('Enter savename','',6,{'file.png'},opts);
% 
% % new figure, same axes [picked preferntially by gui]
% handles.figure2=figure(1);
% handles.axes2 = axes;
% Redoplot(handles);
% 
% % get view options from gui
% pos = get(handles.axes1);
% set(handles.axes2,'View',pos.View);
% 
% print(handles.figure2,fname{:},'-dpng','-r600');
% 
% handles = rmfield(handles,'figure2');
% set(0, 'CurrentFigure', handles.figure1);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ch = allchild(handles.figure1);
cla(ch);drawnow;
% Save the handles structure.
guidata(hObject,handles); 

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rotate3d on;
drawnow;
% Save the handles structure.
guidata(hObject,handles); 
