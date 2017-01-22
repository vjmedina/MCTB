function varargout = display_calibration(varargin)
% DISPLAY_CALIBRATION MATLAB code for display_calibration.fig
%      DISPLAY_CALIBRATION, by itself, creates a new DISPLAY_CALIBRATION or raises the existing
%      singleton*.
%
%      H = DISPLAY_CALIBRATION returns the handle to a new DISPLAY_CALIBRATION or the handle to
%      the existing singleton*.
%
%      DISPLAY_CALIBRATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DISPLAY_CALIBRATION.M with the given input arguments.
%
%      DISPLAY_CALIBRATION('Property','Value',...) creates a new DISPLAY_CALIBRATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before display_calibration_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to display_calibration_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help display_calibration

% Last Modified by GUIDE v2.5 17-Apr-2015 14:48:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @display_calibration_OpeningFcn, ...
                   'gui_OutputFcn',  @display_calibration_OutputFcn, ...
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


% --- Executes just before display_calibration is made visible.
function display_calibration_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to display_calibration (see VARARGIN)

% Choose default command line output for display_calibration
handles.output = hObject;

set(handles.rgb2xyz_matrix_table, 'Data', zeros(3,3));
set(handles.xyz2rgb_matrix_table, 'Data', zeros(3,3));

set(handles.disp_eq_txtfld,'String','y = (gain*x + offset)^gamma');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes display_calibration wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = display_calibration_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in red_xyz_paste_bttn.
function red_xyz_paste_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to red_xyz_paste_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = clipboard('paste');
C = textscan(str,'%f');
A = reshape(C{1,1},3, numel(C{1,1})/3)';

oldData = get(handles.red_xyz_table,'Data');

if (isnumeric(oldData))
    newData = [oldData; A];
else
    newData = A;
end

set(handles.red_xyz_table,'Data',newData);

% --- Executes on button press in red_xyz_clear_bttn.
function red_xyz_clear_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to red_xyz_clear_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.red_xyz_table,'Data',{'','',''});

% --- Executes on button press in green_xyz_paste_bttn.
function green_xyz_paste_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to green_xyz_paste_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = clipboard('paste');
C = textscan(str,'%f');
A = reshape(C{1,1},3, numel(C{1,1})/3)';

oldData = get(handles.green_xyz_table,'Data');

if (isnumeric(oldData))
    newData = [oldData; A];
else
    newData = A;
end

set(handles.green_xyz_table,'Data',newData);

% --- Executes on button press in green_xyz_clear_bttn.
function green_xyz_clear_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to green_xyz_clear_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.green_xyz_table,'Data',{'','',''});

% --- Executes on button press in blue_xyz_paste_bttn.
function blue_xyz_paste_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to blue_xyz_paste_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = clipboard('paste');
C = textscan(str,'%f');
A = reshape(C{1,1},3, numel(C{1,1})/3)';

oldData = get(handles.blue_xyz_table,'Data');

if (isnumeric(oldData))
    newData = [oldData; A];
else
    newData = A;
end

set(handles.blue_xyz_table,'Data',newData);

% --- Executes on button press in blue_xyz_clear_bttn.
function blue_xyz_clear_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to blue_xyz_clear_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.blue_xyz_table,'Data',{'','',''});

% --- Executes on button press in blue_rgb_clear_bttn.
function blue_rgb_clear_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to blue_rgb_clear_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.blue_rgb_table,'Data',{'','',''});

% --- Executes on button press in blue_rgb_paste_bttn.
function blue_rgb_paste_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to blue_rgb_paste_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = clipboard('paste');
C = textscan(str,'%f');
A = reshape(C{1,1},3, numel(C{1,1})/3)';

oldData = get(handles.blue_rgb_table,'Data');

if (isnumeric(oldData))
    newData = [oldData; A];
else
    newData = A;
end

set(handles.blue_rgb_table,'Data',newData);

% --- Executes on button press in green_rgb_clear_bttn.
function green_rgb_clear_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to green_rgb_clear_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.green_rgb_table,'Data',{'','',''});

% --- Executes on button press in green_rgb_paste_bttn.
function green_rgb_paste_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to green_rgb_paste_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = clipboard('paste');
C = textscan(str,'%f');
A = reshape(C{1,1},3, numel(C{1,1})/3)';

oldData = get(handles.green_rgb_table,'Data');

if (isnumeric(oldData))
    newData = [oldData; A];
else
    newData = A;
end

set(handles.green_rgb_table,'Data',newData);

% --- Executes on button press in red_rgb_clear_bttn.
function red_rgb_clear_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to red_rgb_clear_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.red_rgb_table,'Data',{'','',''});

% --- Executes on button press in red_rgb_paste_bttn.
function red_rgb_paste_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to red_rgb_paste_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = clipboard('paste');
C = textscan(str,'%f');
A = reshape(C{1,1},3, numel(C{1,1})/3)';

oldData = get(handles.red_rgb_table,'Data');

if (isnumeric(oldData))
    newData = [oldData; A];
else
    newData = A;
end

set(handles.red_rgb_table,'Data',newData);

% --- Executes on selection change in calib_model_menu.
function calib_model_menu_Callback(hObject, eventdata, handles)
% hObject    handle to calib_model_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns calib_model_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from calib_model_menu
% show_calib_params(handles);

model_index = get(hObject,'Value');

switch(model_index)
    case 1 %GOG
%         set(handles.offset2_txtfld,'Visible','Off');
%         set(handles.offset2_label,'Visible','Off');
        set(handles.disp_eq_txtfld,'String','y = (gain*x + offset)^gamma');
    case 2 %GOGO
%         set(handles.offset1_label,'String','Offset 1:');
        set(handles.disp_eq_txtfld,'String','y = (gain*x + offset1)^gamma + offset2');
     case 3 %GGO
%         set(handles.offset2_txtfld,'Visible','Off');
%         set(handles.offset2_label,'Visible','Off');
        set(handles.disp_eq_txtfld,'String','y = (gain*x)^gamma + offset');
    case 4 %Simple model
%         set(handles.offset1_txtfld,'Visible','Off');
%         set(handles.offset1_label,'Visible','Off');
%         set(handles.offset2_txtfld,'Visible','Off');
%         set(handles.offset2_label,'Visible','Off');
%         set(handles.gain_txtfld,'Visible','Off');
%         set(handles.gain_label,'Visible','Off');
        set(handles.disp_eq_txtfld,'String','y = x^gamma');
end
% 
% function show_calib_params(handles)
%     set(handles.gain_txtfld,'Visible','On');
%     set(handles.gain_label,'Visible','On');
%     set(handles.gamma_txtfld,'Visible','On');
%     set(handles.gamma_label,'Visible','On');
%     set(handles.offset1_txtfld,'Visible','On');
%     set(handles.offset1_label,'String','Offset');
%     set(handles.offset1_label,'Visible','On');
%     set(handles.offset2_txtfld,'Visible','On');
%     set(handles.offset2_label,'Visible','On');
%  


% --- Executes during object creation, after setting all properties.
function calib_model_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calib_model_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gain_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to gain_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gain_txtfld as text
%        str2double(get(hObject,'String')) returns contents of gain_txtfld as a double


% --- Executes during object creation, after setting all properties.
function gain_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gain_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gamma_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to gamma_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gamma_txtfld as text
%        str2double(get(hObject,'String')) returns contents of gamma_txtfld as a double


% --- Executes during object creation, after setting all properties.
function gamma_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function offset1_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to offset1_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of offset1_txtfld as text
%        str2double(get(hObject,'String')) returns contents of offset1_txtfld as a double


% --- Executes during object creation, after setting all properties.
function offset1_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to offset1_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function offset2_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to offset2_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of offset2_txtfld as text
%        str2double(get(hObject,'String')) returns contents of offset2_txtfld as a double


% --- Executes during object creation, after setting all properties.
function offset2_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to offset2_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in calibrate_bttn.
function calibrate_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to calibrate_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% red_ok = 0;
% green_ok = 0;
% blue_ok = 0;

model = get(handles.calib_model_menu,'Value');
Mrgb2xyz = zeros(3,3);
model_params_data = zeros(5,3);
model_info = get_calibration_model_info (model);

params = struct('a',0,'b',0,'c',0,'m',0);
calib_data = struct('model',model,'red_rgb_in',zeros(26,3),'red_xyz_out',zeros(26,3),'red_rgb_out',zeros(26,3),...
                    'red_rgb_estim',zeros(26,3),'red_params',params,'red_R2',0,...
                    'green_rgb_in',zeros(26,3),'green_xyz_out',zeros(26,3),'green_rgb_out',zeros(26,3),...
                    'green_rgb_estim',zeros(26,3),'green_params',params,'green_R2',0,...
                    'blue_rgb_in',zeros(26,3),'blue_xyz_out',zeros(26,3),'blue_rgb_out',zeros(26,3),...
                    'blue_rgb_estim',zeros(26,3),'blue_params',params,'blue_R2',0,...
                    'Mrgb2xyz',zeros(3,3),'Mxyz2rgb',zeros(3,3));

% Check for data consistency and construct the color transformation matrix.
channels = {'red', 'green', 'blue'};
primaries = [[255 0 0]; [0 255 0]; [0 0 255]];
for i=1:3
    switch(i)
        case 1
            xyz_out = get(handles.red_xyz_table,'Data');
            rgb_in = get(handles.red_rgb_table,'Data');
            calib_data.red_rgb_in = rgb_in;
            calib_data.red_xyz_out = xyz_out;
        case 2
            xyz_out = get(handles.green_xyz_table,'Data');
            rgb_in = get(handles.green_rgb_table,'Data');            
            calib_data.green_rgb_in = rgb_in;
            calib_data.green_xyz_out = xyz_out;
            
        case 3
            xyz_out = get(handles.blue_xyz_table,'Data');
            rgb_in = get(handles.blue_rgb_table,'Data');            
            calib_data.blue_rgb_in = rgb_in;
            calib_data.blue_xyz_out = xyz_out;
    end
    
    channel_name = channels{i};    
    if (~isa(rgb_in, 'double') || ~isa(xyz_out, 'double'))
        errordlg(sprintf('The measurementes for the %s channel contain errors or missing data',channel_name),'Missing data');
        return;
    else
        primary_position = findTarget(rgb_in',primaries(i,:));
        if (isempty(primary_position))
            errordlg(sprintf('The %s primary is missing from the list',channel_name),'Missing data');
            return;
        else
            Mrgb2xyz(:,i) = xyz_out(primary_position,:);
        end
    end
end

Mxyz2rgb = inv(Mrgb2xyz);

calib_data.Mrgb2xyz = Mrgb2xyz;
calib_data.Mxyz2rgb = Mxyz2rgb;

% Fill in the matrix tables
set(handles.rgb2xyz_matrix_table,'Data',Mrgb2xyz);
set(handles.xyz2rgb_matrix_table,'Data',Mxyz2rgb);

% Convert xyz_out to rgb_out using the matrix.
for j=1:size(calib_data.red_rgb_in,1), calib_data.red_rgb_out(j,:) = (Mxyz2rgb * calib_data.red_xyz_out(j,:)')'; end
for j=1:size(calib_data.green_rgb_in,1), calib_data.green_rgb_out(j,:) = (Mxyz2rgb * calib_data.green_xyz_out(j,:)')'; end
for j=1:size(calib_data.blue_rgb_in,1), calib_data.blue_rgb_out(j,:) = (Mxyz2rgb * calib_data.blue_xyz_out(j,:)')'; end

% Perform data fitting and estimation for each channel's data.

f = fittype(model_info.expression);

[calib_data.red_params, calib_data.red_R2] = estimate_display_model (calib_data.red_rgb_out(:,1), calib_data.red_rgb_in(:,1), model);
model_params_data(:,1)=[calib_data.red_params.a, calib_data.red_params.m, calib_data.red_params.b, calib_data.red_params.c, calib_data.red_R2];
calib_data.red_rgb_estim = eval_expression (f, calib_data.red_rgb_out(:,1), model_info, calib_data.red_params);

[calib_data.green_params, calib_data.green_R2] = estimate_display_model (calib_data.green_rgb_out(:,2), calib_data.green_rgb_in(:,2), model);
model_params_data(:,2)=[calib_data.green_params.a, calib_data.green_params.m, calib_data.green_params.b, calib_data.green_params.c, calib_data.green_R2];
calib_data.green_rgb_estim = eval_expression (f, calib_data.green_rgb_out(:,2), model_info, calib_data.green_params);

[calib_data.blue_params, calib_data.blue_R2] = estimate_display_model (calib_data.blue_rgb_out(:,3), calib_data.blue_rgb_in(:,3), model);
model_params_data(:,3)=[calib_data.blue_params.a, calib_data.blue_params.m, calib_data.blue_params.b, calib_data.blue_params.c, calib_data.blue_R2];
calib_data.blue_rgb_estim = eval_expression (f, calib_data.blue_rgb_out(:,3), model_info, calib_data.blue_params);

% Fill in the parameters table
set(handles.display_params_table,'Data',model_params_data);

handles.calib_data = calib_data;
guidata(hObject, handles);
    

% xyz_out = get(handles.green_xyz_table,'Data');
% rgb_in = get(handles.green_rgb_table,'Data');
% primary_position = findTarget(rgb_in',[0 255 0]);
% if (isempty(primary_position))
%     errordlg('The green primary is missing from the list','Missing data');
%     return;
% else
%     Mrgb2xyz(:,2) = xyz_out(primary_position,:);
% end
% 
% xyz_out = get(handles.blue_xyz_table,'Data');
% rgb_in = get(handles.blue_rgb_table,'Data');
% primary_position = findTarget(rgb_in',[0 0 255]);
% if (isempty(primary_position))
%     errordlg('The blue primary is missing from the list','Missing data');
%     return;
% else
%     Mrgb2xyz(:,3) = xyz_out(primary_position,:);
% end
% 
% Mxyz2rgb = inv(Mrgb2xyz);
% 
% set(handles.rgb2xyz_matrix_table,'Data',Mrgb2xyz);
% set(handles.xyz2rgb_matrix_table,'Data',Mxyz2rgb);

% rgb_in = get(handles.red_rgb_table,'Data');
% rgb_out = rgb_in;
% xyz_out = get(handles.red_xyz_table,'Data');
% if (isa(rgb_in, 'double') && isa(xyz_out, 'double'))    
%     for i=1:size(rgb_in,1), rgb_out(i,:) = (Mxyz2rgb * xyz_out(i,:)')'; end    
%     [params, R2] = estimate_display_model (rgb_out(:,1), rgb_in(:,1), model);
%     model_params_data(:,1)=[params.a, params.b, params.c, params.m, R2];
%     %red_ok = 1;
% else
%     errordlg('The measurementes for the red channel contain errors or missing data','Missing data');
%     return;    
% end
%     
% rgb_in = get(handles.green_rgb_table,'Data');
% xyz_out = get(handles.green_xyz_table,'Data');
% if (isa(rgb_in, 'double') && isa(xyz_out, 'double'))
%     for i=1:size(rgb_in,1), rgb_out(i,:) = (Mxyz2rgb * xyz_out(i,:)')'; end
%     [params, R2] = estimate_display_model (rgb_out(:,2), rgb_in(:,2), model);
%     model_params_data(:,2)=[params.a, params.b, params.c, params.m, R2];
%     %green_ok = 1;
% else
%     errordlg('The measurementes for the green channel contain errors or missing data','Missing data');
%     return;    
% end
% 
% rgb_in = get(handles.blue_rgb_table,'Data');
% xyz_out = get(handles.blue_xyz_table,'Data');
% if (isa(rgb_in, 'double') && isa(xyz_out, 'double'))
%     for i=1:size(rgb_in,1), rgb_out(i,:) = (Mxyz2rgb * xyz_out(i,:)')'; end
%     [params, R2] = estimate_display_model (rgb_out(:,3), rgb_in(:,3), model);
%     model_params_data(:,3)=[params.a, params.b, params.c, params.m, R2];
%     %blue_ok = 1;
% else
%     errordlg('The measurementes for the blue channel contain errors or missing data','Missing data');
%     return;    
% end

% if (red_ok==0 || green_ok == 0 || blue_ok==0)
%     warning_msg = 'The following color channels contain errors or missing data and were skipped from the calculation:';
%     if (red_ok==0), warning_msg = sprintf('%s\n RED',warning_msg); end
%     if (green_ok==0), warning_msg = sprintf('%s\n GREEN',warning_msg); end
%     if (blue_ok==0), warning_msg = sprintf('%s\n BLUE',warning_msg); end
%     warndlg(warning_msg,'Missing data');
% end

% set(handles.display_params_table,'Data',model_params_data);


% --- Executes on button press in show_graph_bttn.
function show_graph_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to show_graph_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Xr = handles.calib_data.red_rgb_in(:,1);
Xg = handles.calib_data.green_rgb_in(:,2);
Xb = handles.calib_data.blue_rgb_in(:,3);

Yr = handles.calib_data.red_rgb_estim;
Yg = handles.calib_data.green_rgb_estim;
Yb = handles.calib_data.blue_rgb_estim;

figure,
hold on,

scatter(Xr,Yr,'r','s','MarkerFaceColor','r');
scatter(Xg,Yg,'g','d','MarkerFaceColor','g');
scatter(Xb,Yb,'b','^','MarkerFaceColor','b');

xlabel('Real Input', 'FontSize', 12);
ylabel('Estimated Input', 'FontSize', 12);
set(gca,'fontSize',12);

grid on,
hold off,


% --- Executes on button press in display_params_copy_bttn.
function display_params_copy_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to display_params_copy_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
copied_values = get(handles.display_params_table,'Data');
num2clip(copied_values);
h = msgbox('The selected values were successfully copied to the clipboard','Operation Completed');


% --- Executes on button press in rgb2xyz_matrix_copy_bttn.
function rgb2xyz_matrix_copy_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to rgb2xyz_matrix_copy_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
copied_values = get(handles.rgb2xyz_matrix_table,'Data');
num2clip(copied_values);
h = msgbox('The selected values were successfully copied to the clipboard','Operation Completed');

% --- Executes on button press in xyz2rgb_matrix_copy_bttn.
function xyz2rgb_matrix_copy_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to xyz2rgb_matrix_copy_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
copied_values = get(handles.xyz2rgb_matrix_table,'Data');
num2clip(copied_values);
h = msgbox('The selected values were successfully copied to the clipboard','Operation Completed');
