function varargout = image_calibration_batch(varargin)
% IMAGE_CALIBRATION_BATCH MATLAB code for image_calibration_batch.fig
%      IMAGE_CALIBRATION_BATCH, by itself, creates a new IMAGE_CALIBRATION_BATCH or raises the existing
%      singleton*.
%
%      H = IMAGE_CALIBRATION_BATCH returns the handle to a new IMAGE_CALIBRATION_BATCH or the handle to
%      the existing singleton*.
%
%      IMAGE_CALIBRATION_BATCH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGE_CALIBRATION_BATCH.M with the given input arguments.
%
%      IMAGE_CALIBRATION_BATCH('Property','Value',...) creates a new IMAGE_CALIBRATION_BATCH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before image_calibration_batch_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to image_calibration_batch_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help image_calibration_batch

% Last Modified by GUIDE v2.5 05-May-2015 18:50:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @image_calibration_batch_OpeningFcn, ...
                   'gui_OutputFcn',  @image_calibration_batch_OutputFcn, ...
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


% --- Executes just before image_calibration_batch is made visible.
function image_calibration_batch_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to image_calibration_batch (see VARARGIN)

% Choose default command line output for image_calibration_batch
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Set the size of the calibration tables
set(handles.cam_matrix_table, 'Data', zeros(3,4));
set(handles.disp_matrix_table, 'Data', zeros(3,3));

%set(handles.offset2_txtfld,'Visible','Off');
%set(handles.offset2_label,'Visible','Off');
set(handles.disp_eq_txtfld,'String','y = (gain*x + offset)^gamma');

% Display the default original and calibrated images
% axes(handles.orig_image_axes);
% im = imread('no_image.png');
% imshow(im, []);
% axes(handles.calib_image_axes);
% im = imread('no_calibration.png');
% imshow(im, []);


% UIWAIT makes image_calibration_batch wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = image_calibration_batch_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in cammat_paste_bttn.
function cammat_paste_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to cammat_paste_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str = clipboard('paste');
C = textscan(str,'%f');
A = reshape(C{1,1},numel(C{1,1})/3,3)';
if isa(A,'double')    
    set(handles.cam_matrix_table,'Data',A);
    handles.cam_matrix = A;
    guidata(hObject, handles);
else
    errordlg('The clipboard does not contain numbers','Wrong type of data');
    return; 
end


% --- Executes on button press in dispmat_paste_bttn.
function dispmat_paste_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to dispmat_paste_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str = clipboard('paste');
C = textscan(str,'%f');
A = reshape(C{1,1},3, numel(C{1,1})/3)';
if isa(A,'double')
    set(handles.disp_matrix_table,'Data',A);
    handles.disp_matrix = A;
    guidata(hObject, handles);
else
    errordlg('The clipboard does not contain numbers','Wrong type of data');
    return; 
end


% --- Executes on button press in calibrate_bttn.
function calibrate_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to calibrate_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

image_list = cellstr(get(handles.orig_image_list,'String'));
num_images = size(image_list,1);

if (num_images > 0)
    if isfield(handles,'cam_matrix')
        if (isfield(handles,'display_params_red') && isfield(handles,'display_params_red') && isfield(handles,'display_params_red'))
            set(handles.calib_image_list,'String',[]);
            %[result, calib_model] = check_display_calibration_model(handles);
            model_index = get(handles.calib_model_menu,'Value');
            display_model = struct('model',model_index,'red_params',handles.display_params_red, 'green_params',handles.display_params_green,'blue_params',handles.display_params_blue);
            handles.display_model = display_model;
            
            calib_file_suffix = get(handles.calib_file_suffix,'String');
            if strcmp(calib_file_suffix,'')
                calib_file_suffix = 'calib';
            end
            
            for current_entry = 1:num_images
                imFilepath = image_list{current_entry};
                imData = imread(imFilepath);

                % Convert RGB values to XYZ using the camera transfer function.
                Ixyz_disp_out = cam_rgb2xyz(imData, handles.cam_matrix);            

                % Adapt XYZ values to the display's range, based on the maximum
                % luminance.
                use_display_lum = get(handles.use_display_lum_chck,'Value');
                if (use_display_lum)
                    max_display_lum = str2num(get(handles.max_display_lum_txtfld,'String'));
                    max_img_lum = max(max(Ixyz_disp_out(:,:,2)));

                    if (max_img_lum > max_display_lum)
                        factor = max_display_lum / max_img_lum;
                        Ixyz_disp_out = Ixyz_disp_out .* factor;
                    end
                end

                % Convert display XYZ values to RGB using the conversion matrix.
                Irgb_disp_out = disp_xyz2rgb (Ixyz_disp_out, handles.disp_matrix);
                Irgb_disp_out(find(Irgb_disp_out<0)) = 0;     

                % Clip all negative values to 0
                %Ixyz(find(Ixyz<0)) = 0;

                % Reverse the display's gamma using the display's fitting model.
                Irgb_disp_in = display_color_transform (Irgb_disp_out, handles.display_model);
                Irgb_disp_in(find(Irgb_disp_in<0)) = 0;
                Irgb_disp_in = uint8(Irgb_disp_in);
                
                
%                 % Convert RGB values to XYZ using the camera transfer function.
%                 Ixyz = cam_rgb2xyz(imData, handles.cam_matrix);
% 
%                 
%                 % Clip all negative values to 0
%                 Ixyz(find(Ixyz<0)) = 0;
% 
%                 % Convert the pixels values back to RGB using the display's
%                 % calibration parameters.
%                 Irgb_d = display_color_transform (Ixyz, handles.display_model);
%                 Irgb_d = uint8(Irgb_d);
                
                Irgb_d = Irgb_disp_in;

                calib_list = cellstr(get(handles.calib_image_list,'String'));
                if cellfun(@isempty,calib_list)
                    new_value = [imFilepath];
                else
                    new_value = [calib_list; imFilepath];
                end
                set(handles.calib_image_list,'String',new_value)
                
                [~, filename, ~] = fileparts(imFilepath);
                newname = sprintf('%s%s',filename,calib_file_suffix);
                calib_filepath = strrep(imFilepath,filename,newname);
                imwrite(Irgb_d,calib_filepath);
                
                status_msg = sprintf('[Image Tools]->[Image Calibration (Batch)] %s...calibrated',imFilepath);    
                display(sprintf('%s',status_msg));
            end
            
            set(handles.orig_image_list,'Value',size(new_value,1));
            guidata(hObject, handles);

        else
            errordlg('The displays calibration parameters have not been set','Missing display parameters');
            return;
        end
    else
        errordlg('The Camera matrix has not been set','Missing matrix');
        return;
    end
else
    errordlg('You must add some images first','Missing matrix');
    return;
end

function [result, calib_model] = check_display_calibration_model(handles)
    model_index = get(handles.calib_model_menu,'Value');
    params = struct('a',1,'b',0,'c',0,'m',1);
    calib_model = struct('model',model_index,'params',params);
    result = 0;
    
    switch (model_index)
        case 1
            if (isnumeric(str2double(get(handles.gain_txtfld,'String'))) ... 
             && isnumeric(str2double(get(handles.gamma_txtfld,'String'))) ... 
             && isnumeric(str2double(get(handles.offset1_txtfld,'String'))))
                
                calib_model.params.a = str2double(get(handles.gain_txtfld,'String'));
                calib_model.params.m = str2double(get(handles.gamma_txtfld,'String'));
                calib_model.params.b = str2double(get(handles.offset1_txtfld,'String'));
                result = 1;
            end
        case 2
            if (isnumeric(str2double(get(handles.gain_txtfld,'String'))) ... 
             && isnumeric(str2double(get(handles.gamma_txtfld,'String'))) ... 
             && isnumeric(str2double(get(handles.offset1_txtfld,'String'))) ... 
             && isnumeric(str2double(get(handles.offset2_txtfld,'String'))))
                
                calib_model.params.a = str2double(get(handles.gain_txtfld,'String'));
                calib_model.params.m = str2double(get(handles.gamma_txtfld,'String'));
                calib_model.params.b = str2double(get(handles.offset1_txtfld,'String'));
                calib_model.params.c = str2double(get(handles.offset2_txtfld,'String'));
                result = 1;
            end
        case 3
            if (isnumeric(str2double(get(handles.gain_txtfld,'String'))) ... 
             && isnumeric(str2double(get(handles.gamma_txtfld,'String'))) ... 
             && isnumeric(str2double(get(handles.offset1_txtfld,'String'))))
                
                calib_model.params.a = str2double(get(handles.gain_txtfld,'String'));
                calib_model.params.m = str2double(get(handles.gamma_txtfld,'String'));
                calib_model.params.b = str2double(get(handles.offset1_txtfld,'String'));
                result = 1;
            end
        case 4
            if (isnumeric(str2double(get(handles.gamma_txtfld,'String'))))                
                calib_model.params.m = str2double(get(handles.gamma_txtfld,'String'));
                result = 1;
            end
    end
   

% --- Executes on button press in open_image_bttn.
function open_image_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to open_image_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [filename, pathname] = uigetfile( ...
    {'*.bmp;*.jpg;*.jpeg;*.png;*.ppm;*.tif;*.tiff',...
     'All image files (*.bmp;*.jpg;*.jpeg;*.png;*.ppm;*.tiff)';
     '*.bmp',  'Bitmap images(*.bmp)'; ...
     '*.jpeg;*.jpg',  'JPEG images(*.jpeg, *.jpg)'; ...
     '*.png',  'Portable Network Graphics images(*.png)'; ...
     '*.tif;*.tiff',  'Tagged Image File Format images (*.tif, *.tiff)'}, ...
     'Pick an image');

    if ~isnumeric(filename)
        filepath = strcat(pathname,filename);
        %handles.filename = filename;
        %handles.pathname = pathname;
        %handles.filepath = filepath;
        
        initial_values = cellstr(get(handles.orig_image_list,'String'));
        if cellfun(@isempty,initial_values)
            new_value = [filepath];
        else
            new_value = [initial_values; filepath];
        end
        set(handles.orig_image_list,'String',new_value)
        set(handles.orig_image_list,'Value',size(new_value,1));
        guidata(hObject, handles);

        status_msg = sprintf('[Image Tools]->[Image Calibration (Batch)] %s...added',filepath);    
        display(sprintf('%s',status_msg));

        %warning('off','images:initSize:adjustingMag');
        
        %imData = imread(handles.filepath);
        %handles.imData = imData;
        %guidata(hObject, handles);
        
        %axes(handles.orig_image_axes);
        %imshow(handles.imData, []);
    end


% --- Executes on selection change in calib_model_menu.
function calib_model_menu_Callback(hObject, eventdata, handles)
% hObject    handle to calib_model_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns calib_model_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from calib_model_menu

%show_calib_params(handles);
model_index = get(hObject,'Value');
switch(model_index)
    case 1
%         set(handles.offset2_txtfld,'Visible','Off');
%         set(handles.offset2_label,'Visible','Off');
        set(handles.disp_eq_txtfld,'String','y = (gain*x + offset)^gamma');
    case 2
%         set(handles.offset1_label,'String','Offset 1:');
        set(handles.disp_eq_txtfld,'String','y = (gain*x + offset1)^gamma + offset2');
%     case 3
%         set(handles.offset2_txtfld,'Visible','Off');
%         set(handles.offset2_label,'Visible','Off');
        set(handles.disp_eq_txtfld,'String','y = (gain*x)^gamma + offset');
    case 4
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


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in disp_params_copy_bttn.
function disp_params_copy_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to disp_params_copy_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = clipboard('paste');
C = textscan(str,'%f');
A = reshape(C{1,1}',3,4)';
if isa(A,'double')    
    set(handles.display_params_table,'Data',A);
    handles.display_params_red = struct('a',A(1,1),'b',A(3,1),'c',A(4,1),'m',A(2,1));
    handles.display_params_green = struct('a',A(1,2),'b',A(3,2),'c',A(4,2),'m',A(2,2));
    handles.display_params_blue = struct('a',A(1,3),'b',A(3,3),'c',A(4,3),'m',A(2,3));
    guidata(hObject, handles);
else
    errordlg('The clipboard does not contain numbers','Wrong type of data');
    return; 
end


% --- Executes on selection change in orig_image_list.
function orig_image_list_Callback(hObject, eventdata, handles)
% hObject    handle to orig_image_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns orig_image_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from orig_image_list


% --- Executes during object creation, after setting all properties.
function orig_image_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to orig_image_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in orig_image_delete_bttn.
function orig_image_delete_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to orig_image_delete_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    selected_indices = get(handles.orig_image_list,'Value');
    initial_values = cellstr(get(handles.orig_image_list,'String'));
    list_size = numel(initial_values);
    new_list = [];
    
    if (numel(initial_values{1,1})~=0)
        for current_index = 1:list_size
            if isempty(find (selected_indices == current_index, 1))
                new_list = [new_list(:); initial_values(current_index)];
            end
        end
        
        set(handles.orig_image_list,'String',new_list);
        if (selected_indices(1)==list_size)
            set(handles.orig_image_list,'Value',numel(new_list));
        else
            set(handles.orig_image_list,'Value',selected_indices(1));
        end
    end


% --- Executes on selection change in calib_image_list.
function calib_image_list_Callback(hObject, eventdata, handles)
% hObject    handle to calib_image_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns calib_image_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from calib_image_list


% --- Executes during object creation, after setting all properties.
function calib_image_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calib_image_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function calib_file_suffix_Callback(hObject, eventdata, handles)
% hObject    handle to calib_file_suffix_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of calib_file_suffix_label as text
%        str2double(get(hObject,'String')) returns contents of calib_file_suffix_label as a double


% --- Executes during object creation, after setting all properties.
function calib_file_suffix_label_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calib_file_suffix_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function calib_file_suffix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calib_file_suffix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in dispmat_paste_bttn.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to dispmat_paste_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function max_display_lum_slider_Callback(hObject, eventdata, handles)
% hObject    handle to max_display_lum_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function max_display_lum_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_display_lum_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function max_display_lum_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to max_display_lum_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_display_lum_txtfld as text
%        str2double(get(hObject,'String')) returns contents of max_display_lum_txtfld as a double


% --- Executes during object creation, after setting all properties.
function max_display_lum_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_display_lum_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in display_lum_set_bttn.
function display_lum_set_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to display_lum_set_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in use_display_lum_chck.
function use_display_lum_chck_Callback(hObject, eventdata, handles)
% hObject    handle to use_display_lum_chck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_display_lum_chck

checked = get(hObject,'Value');
if (checked)
    set (handles.max_display_lum_slider,'Enable','On');
    set (handles.max_display_lum_txtfld,'Enable','On');
    set (handles.display_lum_set_bttn,'Enable','On');
else
    set (handles.max_display_lum_slider,'Enable','Off');
    set (handles.max_display_lum_txtfld,'Enable','Off');
    set (handles.display_lum_set_bttn,'Enable','Off');
end
    
