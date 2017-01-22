function varargout = colorimetric_tests(varargin)
% COLORIMETRIC_TESTS MATLAB code for colorimetric_tests.fig
%      COLORIMETRIC_TESTS, by itself, creates a new COLORIMETRIC_TESTS or raises the existing
%      singleton*.
%
%      H = COLORIMETRIC_TESTS returns the handle to a new COLORIMETRIC_TESTS or the handle to
%      the existing singleton*.
%
%      COLORIMETRIC_TESTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COLORIMETRIC_TESTS.M with the given input arguments.
%
%      COLORIMETRIC_TESTS('Property','Value',...) creates a new COLORIMETRIC_TESTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before colorimetric_tests_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to colorimetric_tests_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help colorimetric_tests

% Last Modified by GUIDE v2.5 21-May-2015 18:10:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @colorimetric_tests_OpeningFcn, ...
                   'gui_OutputFcn',  @colorimetric_tests_OutputFcn, ...
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


% --- Executes just before colorimetric_tests is made visible.
function colorimetric_tests_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to colorimetric_tests (see VARARGIN)

% Choose default command line output for colorimetric_tests
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
axes(handles.orig_image_axes);
im = imread('no_image.png');
imshow(im, []);
axes(handles.output1_image_axes);
im = imread('no_calibration.png');
imshow(im, []);
axes(handles.output2_image_axes);
imshow(im, []);
axes(handles.output3_image_axes);
imshow(im, []);


% UIWAIT makes colorimetric_tests wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = colorimetric_tests_OutputFcn(hObject, eventdata, handles) 
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

if isfield(handles,'imData')
    if isfield(handles,'cam_matrix') && isfield(handles,'disp_matrix')
        if (isfield(handles,'display_params_red') && isfield(handles,'display_params_red') && isfield(handles,'display_params_red'))
            white_index = get(handles.ref_white_list,'Value');
            switch(white_index)
                case 1
                    handles.illuminant = ColorUtils.ILLUM_CIE_D65;
                    handles.illuminant_LMS = ColorUtils.XYZ2LMS([ColorUtils.CIE_D65_xy10(1) ColorUtils.CIE_D65_xy10(2) 100], ColorUtils.CAT_VONKRIES);                    
                case 2
                    handles.illuminant = ColorUtils.ILLUM_CUSTOM;                    
                    l = str2num(get(handles.white_l_txtfld,'String'));
                    m = str2num(get(handles.white_m_txtfld,'String'));
                    s = str2num(get(handles.white_s_txtfld,'String'));
                    
                    if (isempty(l) || isempty(m) || isempty(s))
                        errordlg('Some of the reference white''s tristimulus values are missing');
                        return;
                    end
                    
                    handles.illuminant_LMS = [l m s];
            end
                
            %[result, calib_model] = check_display_calibration_model(handles);
            model_index = get(handles.calib_model_menu,'Value');
            display_model = struct('model',model_index,'red_params',handles.display_params_red, 'green_params',handles.display_params_green,'blue_params',handles.display_params_blue);
            handles.display_model = display_model;
            
            % Convert RGB values to LMS using the camera transfer function.
            ILMS_disp_out = cam_rgb2xyz(handles.imData, handles.cam_matrix);
            
            maxLMS_orig = max(max(ILMS_disp_out));
            set(handles.orig_maxXYZ_txtfld,'String',mat2str(maxLMS_orig(:),7));
            
            % Adapt LMS values to the display's range, based on the maximum
            % luminance.
            use_display_lum = get(handles.use_display_lum_chck,'Value');
            if (use_display_lum)
                max_display_lum = str2num(get(handles.max_display_lum_txtfld,'String'));
                max_img_lum = max(max(ILMS_disp_out(:,:,2)));
                
                if (max_img_lum > max_display_lum)
                    factor = max_display_lum / max_img_lum;
                    ILMS_disp_out = ILMS_disp_out .* factor;
                end
            end
            
            ILMS_disp_out(find(ILMS_disp_out<0)) = 0;
            L = ILMS_disp_out(:,:,1);
            M = ILMS_disp_out(:,:,2);
            S = ILMS_disp_out(:,:,3);
            
            model = get(get(handles.color_model_bttn_group,'SelectedObject'),'Tag');
            switch model
                case 'model2'
                    color_space=2;
                    % We ignore the sign because it only indicates which color 
                    % it's closer to in the Red-Green line. A negative
                    % value indicates that it's closer to red than a positive 
                    % one (closer to green) or the other way around.
                    %
                    % Red ------- 0 ++++++++ Green
                    %
                    ILMS_disp_out(:,:,1) = L + M;
                    ILMS_disp_out(:,:,2) = (L - (2*M)); 
                    ILMS_disp_out(:,:,3) = (L + M - S);                    
                    set(handles.output1_label,'String','Luminance: L + M');
                    set(handles.output2_label,'String','Red-Green opponency: L - 2M');
                    set(handles.output3_label,'String','Yellow-Blue opponency: L + M - S');
                case 'model3'
                    color_space=3;
                    set(handles.output1_label,'String','Hue (H)');
                    set(handles.output2_label,'String','Saturation (S)');
                    set(handles.output3_label,'String','Value (V)');
                otherwise
                    color_space=0;
                    set(handles.output1_label,'String','Long wavelengths (L)');
                    set(handles.output2_label,'String','Medium wavelengths (M)');
                    set(handles.output3_label,'String','Short wavelengths (S)');
            end
                   
            % Color balance the LMS using Von Kries' method to adapt the 
            % HVS response to the environment's illuminant.
            %M_white = eye(3,3)./repmat(handles.illuminant_LMS',1,3);
            %ILMS_disp_out = img_array_mult(ILMS_disp_out, M_white);
            
            % We do not have the LMS or XYZ tristimulus values of a perfect
            % white (illuminant's white point) under each illumination. We
            % only have those of the Macbeth color chart's white, but since
            % the chart and samples plates are made of different materials 
            % we often find much brighter pixels in plate images than in
            % the calibration color chart. Therefore, to avoid
            % overexposure, we normalize using the maximum values between
            % the Macbeth's white and the image's maximum values. 
            %
            
            %max_LMS_img = reshape(max(max(ILMS_disp_out)),1,3);
            for i=1:3
                %ILMS_disp_out(:,:,i) = ILMS_disp_out(:,:,i)./handles.illuminant_LMS(i);
                %%ILMS_disp_out(:,:,i) = ILMS_disp_out(:,:,i)./max_LMS_img(i);
            end
            
            
            %% Convert display LMS values to RGB using the conversion matrix.
            Irgb_disp_out = disp_xyz2rgb (ILMS_disp_out, handles.disp_matrix);
            Irgb_disp_out(find(Irgb_disp_out<0)) = 0;
            %Irgb_disp_in = ILMS_disp_out;
            %for i=1:3
                %Irgb_disp_in(:,:,1) = ILMS_disp_out(:,:,i) - min(min(ILMS_disp_out(:,:,i)));                
                %Irgb_disp_in(:,:,1) = Irgb_disp_in(:,:,i) ./ 0.2;
                %Irgb_disp_in(:,:,1) = Irgb_disp_in(:,:,i) .* 255;
            %end
            
            %Irgb_disp_in(:,:,2) = (ILMS_disp_out(:,:,2) + 0.05) .* 2550;
            %Irgb_disp_in(:,:,3) = (ILMS_disp_out(:,:,3) + 0.05) .* 2550;
            
            % Clip all negative values to 0
            %Ixyz(find(Ixyz<0)) = 0;

            %% Reverse the display's gamma using the display's fitting model.
            Irgb_disp_in = display_color_transform (Irgb_disp_out, handles.display_model);

            if color_space == 3
                %[sx, sy, sz] = size(Irgb_disp_in);
                %R = double(Irgb_disp_in(:,:,1));
                %G = double(Irgb_disp_in(:,:,2));
                %B = double(Irgb_disp_in(:,:,3));                
                
                %temp = ColorUtils.RGB2HSV(R(:), G(:), B(:));
                %Irgb_disp_in = reshape(temp,sx,sy,sz);
                Irgb_disp_in = rgb2hsv(Irgb_disp_in);
                
            end
            
            Irgb_disp_in(find(Irgb_disp_in<0)) = 0;
            Irgb_disp_in = uint8(Irgb_disp_in);
            
            handles.imOutput1Data = Irgb_disp_in(:,:,1);
            handles.imOutput2Data = Irgb_disp_in(:,:,2);
            handles.imOutput3Data = Irgb_disp_in(:,:,3);

            axes(handles.output1_image_axes);
            imshow(uint8(handles.imOutput1Data));
            axes(handles.output2_image_axes);
            imshow(uint8(handles.imOutput2Data));
            axes(handles.output3_image_axes);
            imshow(uint8(handles.imOutput3Data));

            guidata(hObject, handles);
            set (handles.save_output1_bttn,'Enable','On');
            set (handles.save_output2_bttn,'Enable','On');
            set (handles.save_output3_bttn,'Enable','On');
            
            maxXYZ = max(max(ILMS_disp_out));
            set(handles.calib_maxXYZ_txtfld,'String',mat2str(maxXYZ(:),7));
                        
            maxRGB = max(max(Irgb_disp_in));
            set(handles.calib_maxRGB_txtfld,'String',mat2str(maxRGB(:),7));
        else
            errordlg('The displays calibration parameters have not been set','Missing display parameters');
            return;
        end
    else
        errordlg('The camera and/or display matrix has not been set','Missing matrix');
        return;
    end
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
        handles.filename = filename;
        handles.pathname = pathname;
        handles.filepath = filepath;

        status_msg = sprintf('[Image Tools]->[Image Calibration] %s...loaded',filepath);    
        display(sprintf('%s',status_msg));

        warning('off','images:initSize:adjustingMag');
        
        imData = imread(handles.filepath);
        handles.imData = imData;
        guidata(hObject, handles);
        
        axes(handles.orig_image_axes);
        imshow(handles.imData, []);
        
        set (handles.save_output1_bttn,'Enable','Off');
            
        axes(handles.output1_image_axes);
        im = imread('no_calibration.png');
        imshow(im, []);
        
        maxRGB = max(max(imData));
        set(handles.orig_maxRGB_txtfld,'String',mat2str(maxRGB(:),7));
        
        set(handles.orig_maxXYZ_txtfld,'String','');
        set(handles.calib_maxXYZ_txtfld,'String','');
        set(handles.calib_maxRGB_txtfld,'String','');
        
%         maxXYZ = max(max(imData));
%         set(handles.orig_maxXYZ_txtfld,'String',mat2str(maxXYZ(:)));
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
A = reshape(C{1,1}',3,5)';
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


% --- Executes on button press in save_output1_bttn.
function save_output1_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to save_output1_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles.imData = imData;
% handles.filename = filename;
%         handles.pathname = pathname;
%         handles.filepath = filepath;
        
        [~, filename, ~] = fileparts(handles.filepath);
        newname = sprintf('%s%s',filename,'_output1');
        output_filepath = strrep(handles.filepath,filename,newname);
        %newname = sprintf('%s\%s_calib%s',pathstr,filename,ext);
        %newname = sprintf('%s%s',handles.filename,calib_file_suffix);
        [fn,pathname] = uiputfile('*.bmp; *.jpg; *.png','Save image as...',output_filepath);
        if (fn ~= 0)
            output_filepath = sprintf('%s%s',pathname,fn);
            imwrite(handles.imOutput1Data,output_filepath);
            %status_msg = sprintf('[Image Tools]->[Image Calibration (Batch)] %s...saved',calib_filepath);
            %display(sprintf('%s',status_msg));
        end


% --- Executes on button press in dispmat_paste_bttn.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to dispmat_paste_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on slider movement.
function max_display_lum_slider_Callback(hObject, eventdata, handles)
% hObject    handle to max_display_lum_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
max_display_lum = get (hObject,'value');
set(handles.max_display_lum_txtfld,'String',max_display_lum);

% --- Executes during object creation, after setting all properties.
function max_display_lum_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_display_lum_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


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
max_display_lum = str2num(get (handles.max_display_lum_txtfld,'String'));
set(handles.max_display_lum_slider,'Value',max_display_lum);


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
    


% --- Executes on button press in save_output2_bttn.
function save_output2_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to save_output2_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[~, filename, ~] = fileparts(handles.filepath);
        newname = sprintf('%s%s',filename,'_output1');
        output_filepath = strrep(handles.filepath,filename,newname);
        %newname = sprintf('%s\%s_calib%s',pathstr,filename,ext);
        %newname = sprintf('%s%s',handles.filename,calib_file_suffix);
        [fn,pathname] = uiputfile('*.bmp; *.jpg; *.png','Save image as...',output_filepath);
        if (fn ~= 0)
            output_filepath = sprintf('%s%s',pathname,fn);
            imwrite(handles.imOutput2Data,output_filepath);
            %status_msg = sprintf('[Image Tools]->[Image Calibration (Batch)] %s...saved',calib_filepath);
            %display(sprintf('%s',status_msg));
        end

% --- Executes on button press in save_output3_bttn.
function save_output3_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to save_output3_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

        [~, filename, ~] = fileparts(handles.filepath);
        newname = sprintf('%s%s',filename,'_output1');
        output_filepath = strrep(handles.filepath,filename,newname);
        %newname = sprintf('%s\%s_calib%s',pathstr,filename,ext);
        %newname = sprintf('%s%s',handles.filename,calib_file_suffix);
        [fn,pathname] = uiputfile('*.bmp; *.jpg; *.png','Save image as...',output_filepath);
        if (fn ~= 0)
            output_filepath = sprintf('%s%s',pathname,fn);
            imwrite(handles.imOutput3Data,output_filepath);
            %status_msg = sprintf('[Image Tools]->[Image Calibration (Batch)] %s...saved',calib_filepath);
            %display(sprintf('%s',status_msg));
        end


% --- Executes when selected object is changed in color_model_bttn_group.
function color_model_bttn_group_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in color_model_bttn_group 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
%calibrate_bttn_Callback(hObject, eventdata, handles);



function white_l_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to white_l_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of white_l_txtfld as text
%        str2double(get(hObject,'String')) returns contents of white_l_txtfld as a double


% --- Executes during object creation, after setting all properties.
function white_l_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to white_l_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function white_m_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to white_m_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of white_m_txtfld as text
%        str2double(get(hObject,'String')) returns contents of white_m_txtfld as a double


% --- Executes during object creation, after setting all properties.
function white_m_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to white_m_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function white_s_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to white_s_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of white_s_txtfld as text
%        str2double(get(hObject,'String')) returns contents of white_s_txtfld as a double


% --- Executes during object creation, after setting all properties.
function white_s_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to white_s_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
     


% --- Executes on selection change in ref_white_list.
function ref_white_list_Callback(hObject, eventdata, handles)
% hObject    handle to ref_white_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ref_white_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ref_white_list
index = get(hObject,'Value');
switch(index)
    case 1
        set(handles.white_l_txtfld,'Enable','Off');
        set(handles.white_m_txtfld,'Enable','Off');
        set(handles.white_s_txtfld,'Enable','Off');
    case 2
        set(handles.white_l_txtfld,'Enable','On');
        set(handles.white_m_txtfld,'Enable','On');
        set(handles.white_s_txtfld,'Enable','On');
end

% --- Executes during object creation, after setting all properties.
function ref_white_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ref_white_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
