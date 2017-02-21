function varargout = contrast_gui(varargin)
% t-Contrast sub-gui for VSExtractor_gui
%
%
%
%
% AS2016

% Edit the above text to modify the response to help contrast_gui

% Last Modified by GUIDE v2.5 30-Nov-2016 15:48:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @contrast_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @contrast_gui_OutputFcn, ...
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


% --- Executes just before contrast_gui is made visible.
function contrast_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to contrast_gui (see VARARGIN)

% Choose default command line output for contrast_gui
handles.output = hObject;

global guihand
set(handles.popupmenu1, 'String', guihand.G{1}.condlist);
set(handles.popupmenu2, 'String', guihand.G{1}.condlist);


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes contrast_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = contrast_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

% popup condition 1:
global guihand
%current_entries = cellstr(get(handles.popupmenu1, 'String'));
%current_entries{end+1} = Addition;
set(handles.popupmenu1, 'String', guihand.G{1}.condlist);
contents = cellstr(get(hObject,'String'));
Sel = contents{get(hObject,'Value')};
Sel = find(strcmp(Sel,guihand.G{1}.condlist));

handles.T1 = Sel;

% Save the handles structure.
guidata(hObject,handles); 


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

% string conition 1:
handles.T1 = get(hObject,'String');

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


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2

% popup condition 2:
global guihand
%current_entries = cellstr(get(handles.popupmenu1, 'String'));
%current_entries{end+1} = Addition;
set(handles.popupmenu2, 'String', guihand.G{1}.condlist);
contents = cellstr(get(hObject,'String'));
Sel = contents{get(hObject,'Value')};
Sel = find(strcmp(Sel,guihand.G{1}.condlist));

handles.T2 = Sel;

% Save the handles structure.
guidata(hObject,handles); 


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

% string condition 2:
handles.T2 = get(hObject,'String');

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



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double

% time of interest:
woi = (get(hObject,'String'));
woi = str2num(woi);
handles.woi = woi;

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
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double

% freqs of interest:
foi = (get(hObject,'String'));
foi = str2num(foi);
handles.foi = foi;

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


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% run the function and do plot
global guihand

G   = guihand.G;
T1  = handles.T1; 
T2  = handles.T2;
try toi = handles.woi; catch; toi = []; end
try foi = handles.foi; catch; foi = []; end

fprintf('Trying to compute t-contrasts...\n');

if iscell(T1) && length(T1) > 1;
    % we are running multiple contrasts
    if length(T1) ~= length(T2);
        error('Must specify contrasts pairs!');
    end
    
    for k = 1:length(T1)
        fprintf('Running contrast %d of %d\n',k,length(T1));
        tmap{k} = contrastmesh(G,T1{k},T2{k},toi,foi);
    end
    if isfield(guihand,'trans');trans= guihand.trans;else trans= []; end
    
    global ConPlot
    ConPlot.trans = trans;
    ConPlot.nc    = length(T1);
    ConPlot.G = G{1};
    ConPlot.tmap = tmap;
    
    contrast_plot_gui
    
%     % new figure but to inherit plot settings
%     ch = figure(1);
%     nc = length(T1);
%     if isfield(guihand,'trans');trans= guihand.trans;else trans= []; end
%     global thr
%     global trs
%     
%     for k = 1:nc
%         % t ->
%         subplot(2,nc,k), plotmesh_fo_tmap(G{1},tmap{k},thr,trs,trans);
%         % t <-
%         subplot(2,nc,k+nc),plotmesh_fo_tmap(G{1},tmap{k}*-1,thr,trs,trans);
%     end
%     ac = allchild(ch);
%     linksubplots(ac);
    
else
    % we are running only 1 contrast
    tmap = contrastmesh(G,T1,T2,toi,foi);
    
    if ~isempty(tmap)
        
        % switch to VSExtractor_gui figure axes
        %----------------------------------------
        set(0, 'CurrentFigure', guihand.figure1)
        axes(guihand.axes1);
        
        if isfield(guihand,'trans');trans= guihand.trans;else trans= []; end
        global thr
        global trs
        
        plotmesh_fo_tmap(G{1},tmap,thr,trs,trans)
        
    end
end

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.T1 = [];
handles.T2 = [];

% Save the handles structure.
guidata(hObject,handles); 



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double

% List t(1)s
T1 = get(hObject,'String');
handles.T1 = cellstr(T1);


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



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double

% List t(2)s
T2 = get(hObject,'String');
handles.T2 = cellstr(T2);

% Save the handles structure.
guidata(hObject,handles); 


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
