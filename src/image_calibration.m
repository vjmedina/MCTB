function varargout = image_calibration(varargin)
% IMAGE_CALIBRATION MATLAB code for image_calibration.fig
%      IMAGE_CALIBRATION, by itself, creates a new IMAGE_CALIBRATION or raises the existing
%      singleton*.
%
%      H = IMAGE_CALIBRATION returns the handle to a new IMAGE_CALIBRATION or the handle to
%      the existing singleton*.
%
%      IMAGE_CALIBRATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGE_CALIBRATION.M with the given input arguments.
%
%      IMAGE_CALIBRATION('Property','Value',...) creates a new IMAGE_CALIBRATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before image_calibration_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to image_calibration_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help image_calibration

% Last Modified by GUIDE v2.5 08-Jan-2016 18:24:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @image_calibration_OpeningFcn, ...
                   'gui_OutputFcn',  @image_calibration_OutputFcn, ...
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


% --- Executes just before image_calibration is made visible.
function image_calibration_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to image_calibration (see VARARGIN)

% Choose default command line output for image_calibration
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
im = imread('images/no_image.png');
imshow(im, []);
axes(handles.computed_image_axes);
im = imread('images/no_calibration.png');
imshow(im, []);

%icon=imread('eyedropper_icon.png');
%set(handles.eyedropper_bttn,'CData',icon);

% UIWAIT makes image_calibration wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = image_calibration_OutputFcn(hObject, eventdata, handles) 
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


% 1. CHECK THAT AN INPUT IMAGE HAS BEEN ENTERED AND GIVE AN ERROR OTHERWISE.

    if ~isfield(handles,'imData')
       errordlg('Please select an input image fisrt','Input image not selected');
       return;
    end

% 2. COMPUTE THE XYZ VALUES CORRESPONDING TO THE INPUT IMAGE ACCORDING TO 
%    THE CAMERA'S CHARACTERIZATION MATRIX.

    %if ~isfield(handles,'ImxyzData_orig')
        % Convert RGB values to XYZ using the camera transfer function,
        % and make all the neccesary modifications to the image
        % according to the options selected in the GUI.
        [res, Ixyz_disp_out_orig, Ixyz_disp_out_resampled] = obtain_cam_xyz_image(handles.imData, hObject, handles);
        if (res==1)
            handles.ImxyzData_orig = Ixyz_disp_out_orig;
            handles.ImxyzData = Ixyz_disp_out_resampled;
        else
            if (res==-1)
                errordlg('The camera matrix has not been set','Missing camera matrix');
                return;
            end
        end
    %else
    %    Ixyz_disp_out_orig = handles.ImxyzData_orig;
    %    Ixyz_disp_out_resampled = handles.ImxyzData_orig;
    %end

% 3. COMPUTE THE DISPLAY'S INPUT RGB REQUIRED TO OBTAIN THE DESIRED XYY
%    OUTPUT, ACCORDING TO THE DISPLAY'S CHARACTERIZATION DATA.
 
    % Has the display characterization model been entered?
    if isfield(handles,'disp_matrix')
        % Have the parameters of the display fitting model been entered?
        if (isfield(handles,'display_params_red') && isfield(handles,'display_params_green') && isfield(handles,'display_params_blue'))
            % Read the display characterization model parameters.            
            model_index = get(handles.calib_model_menu,'Value');
            display_model = struct('model',model_index,'red_params',handles.display_params_red, 'green_params',handles.display_params_green,'blue_params',handles.display_params_blue);
            handles.display_model = display_model;
        else
            errordlg('The display fitting model parameters have not been set','Missing display parameters');
            return;
        end
    else
       errordlg('The display matrix has not been set','Missing display matrix');
       return;
    end


    % Convert display XYZ values to display RGB using the 
    % conversion matrix.
    Irgb_disp_out = disp_xyz2rgb (Ixyz_disp_out_resampled, handles.disp_matrix);
    Irgb_disp_out(find(Irgb_disp_out<0)) = 0;     

    % Reverse the display's gamma using the display's fitting model.
    Irgb_disp_in = display_color_transform (Irgb_disp_out, handles.display_model);
    Irgb_disp_in(find(Irgb_disp_in<0)) = 0;
    Irgb_disp_in = uint8(Irgb_disp_in);

    % Save the calibrated data as a global variable.
    handles.imComputedData = Irgb_disp_in;

    % Display the calibrated image in the corresponding axes.
    handles.operation_mode = 0; % 0 is for display RGB images, and 1 for scene XYZ.
    axes(handles.computed_image_axes);
    imPlot = imshow(uint8(handles.imComputedData), []);
    handles.imPlot = imPlot;

% 4.OBTAIN ALL THE LUMINANCE VALUES FROM THE ORIGINAL AND MODIFIED XYZ 
%   IMAGES AND DISPLAY THEIR HISTOGRAMS WITH "nbins" BINS AS A BAR PLOT IN
%   THE CORRESPONDING AXES.
% 
%   This is done for experimentation purposes.

    nbins = 100;
    im = Ixyz_disp_out_orig(:,:,2); im = im(:);
    axes(handles.xyz1_histogram_axes);            
    [nb,xb] = hist(im,nbins);
    bh = bar(xb,nb);
    set(bh,'facecolor',[1 0 0],'EdgeColor',[0.7 0 0],'LineWidth',0.25);
    set(gca,'FontSize',6);

    im = Ixyz_disp_out_resampled(:,:,2); im = im(:);
    axes(handles.xyz2_histogram_axes);
    [nb,xb] = hist(im,nbins);
    bh = bar(xb,nb);
    set(bh,'facecolor',[0 0 1],'EdgeColor',[0 0 0.7],'LineWidth',0.25);
    set(gca,'FontSize',6);
    drawnow,


    handles.vis_draw_space = 0;  % 0 is for RGB color space, and 1 for XYZ.
    handles.sel_channel = 4;
    
    % Commit changes in global variables so they can be available
    % to other functions and modules.
    guidata(hObject, handles);

% 5. ONCE THE OUTPUT IMAGE HAS BEEN COMPUTED AND DISPLAYED, WE CAN ENABLE
%    PERFORM SOME ACTIONS REQUIRED TO UPDATE THE GUI.
    
    gui_post_image_displayed (handles);
    

% After an image is displayed, we perform some action required to update 
% GUI, such as enabling the "Save" button or display some image statistics.
function [] = gui_post_image_displayed (handles)
    set(handles.channel_selection_grp,'Visible','on');
%     set(handles.channel1_chkbx,'Visible','on');
%     set(handles.channel2_chkbx,'Visible','on');
%     set(handles.channel3_chkbx,'Visible','on');
    
    if handles.operation_mode == 0 %Display RGB
        set(handles.color_space_grp,'Visible','on');
        set(handles.color_space_grp,'selectedobject',handles.color_space_rgb_radio);
        set(handles.channel_selection_grp,'selectedobject',handles.channel4_chkbx);
        set(handles.color_space_rgb_radio,'Enable','on');
        set(handles.channel1_chkbx,'String','R');
        set(handles.channel2_chkbx,'String','G');
        set(handles.channel3_chkbx,'String','B');
        set(handles.channel4_chkbx,'Visible','on');
    else
        if handles.operation_mode == 1 %Scene XYZ
            set(handles.color_space_grp,'Visible','on');
            set(handles.color_space_grp,'selectedobject',handles.color_space_xyz_radio);
            set(handles.channel_selection_grp,'selectedobject',handles.channel1_chkbx);
            set(handles.color_space_rgb_radio,'Enable','off');
            set(handles.channel1_chkbx,'String','X');
            set(handles.channel2_chkbx,'String','Y');
            set(handles.channel3_chkbx,'String','Z');
            set(handles.channel4_chkbx,'Visible','off');
        end
    end
    
    set(handles.save_image_bttn,'Enable','On');

    maxXYZ_calib = max(max(handles.ImxyzData));
    set(handles.calib_maxXYZ_txtfld,'String',mat2str(maxXYZ_calib(:),4));           

    maxRGB = max(max(handles.imComputedData));
    set(handles.calib_maxRGB_txtfld,'String',mat2str(maxRGB(:),4));

    meanRGB = mean(mean(handles.imComputedData));
    set(handles.calib_meanRGB_txtfld,'String',mat2str(meanRGB(:),4));

    a=reshape(handles.imComputedData(:,:,1)>=255,[],1);
    b=reshape(handles.imComputedData(:,:,2)>=255,[],1);
    c=reshape(handles.imComputedData(:,:,3)>=255,[],1);
    saturation=([sum(a) sum(b) sum(c)]./numel(a)).*100;
    set(handles.calib_satRGB_txtfld,'String',strcat(mat2str(saturation(:),3),'%'));

    a=reshape(handles.imComputedData(:,:,1)>=255,[],1);
    b=reshape(handles.imComputedData(:,:,2)>=255,[],1);
    c=reshape(handles.imComputedData(:,:,3)>=255,[],1);
    saturation=([sum(a) sum(b) sum(c)]./numel(a)).*100;
    set(handles.calib_satRGB_txtfld,'String',strcat(mat2str(saturation(:),3),'%'));


function [res, Ixyz_disp_out_orig, Ixyz_disp_out_resampled] = obtain_cam_xyz_image(imData, hObject, handles)
    % Has the camera characterization model been entered?
    if isfield(handles,'cam_matrix')
        % Convert RGB values to XYZ using the camera transfer function.
        Ixyz_disp_out_orig = cam_rgb2xyz(imData, handles.cam_matrix);            

        % Set the "max XYZ" info text for the input image.
        maxXYZ_orig = max(max(Ixyz_disp_out_orig));
        set(handles.orig_maxXYZ_txtfld,'String',mat2str(maxXYZ_orig(:),4));

        % Translate the XYZ values to modify the dynamic range
        % according to the given factor. If no factor is given, the
        % default value is 1, which means the image does not change.
        %
        % This is used for experimentation purposes, to test the effect
        % of different dynamic ranges.
        resampling_factor = str2num(get(handles.resampling_factor_fld,'String'));
        handles.resampling_factor = resampling_factor;
        Ixyz_disp_out_resampled = Ixyz_disp_out_orig .* resampling_factor;

        % Adapt XYZ values to the display's range, according to its 
        % maximum achievable luminance, so that all values in the image
        % can be represented by the display.
        use_display_lum = get(handles.use_display_lum_chck,'Value');
        if (use_display_lum)
            max_display_lum = str2num(get(handles.max_display_lum_txtfld,'String'));
            max_img_lum = max(max(Ixyz_disp_out_resampled(:,:,2)));

            if (max_img_lum > max_display_lum)
                factor = max_display_lum / max_img_lum;
                Ixyz_disp_out_resampled = Ixyz_disp_out_resampled .* factor;
            end
        end

        % Save both the original and modified image as global variables
        % for future use. 
        %handles.ImxyzData_orig = Ixyz_disp_out_orig;
        %handles.ImxyzData = Ixyz_disp_out_resampled;

        % Commit changes in global variables so they can be available
        % to other functions and modules.
        guidata(hObject, handles);
        
        res = 1;
    else
       Ixyz_disp_out_orig = [];
       Ixyz_disp_out_resampled = [];
       res = -1; % Return error -1: Missing camera matrix.
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
        
        set (handles.save_image_bttn,'Enable','Off');
            
        axes(handles.computed_image_axes);
        im = imread('no_calibration.png');
        imshow(im, []);
        
        maxRGB = uint16(max(max(imData)));
        set(handles.orig_maxRGB_txtfld,'String',mat2str(maxRGB(:),7));
        
        meanRGB = uint16(mean(mean(imData)));
        set(handles.orig_meanRGB_txtfld,'String',mat2str(meanRGB(:),7));
        
        a=reshape(imData(:,:,1)>=65535,[],1);
        b=reshape(imData(:,:,2)>=65535,[],1);
        c=reshape(imData(:,:,3)>=65535,[],1);
        saturation=([sum(a) sum(b) sum(c)]./numel(a)).*100;
        set(handles.orig_satRGB_txtfld,'String',strcat(mat2str(saturation(:),3),'%'));
        
        set(handles.orig_maxXYZ_txtfld,'String','');
        set(handles.calib_maxXYZ_txtfld,'String','');
        set(handles.calib_maxRGB_txtfld,'String','');
        
        set(handles.resampling_factor_fld,'String','1');
        
        set(handles.color_space_grp,'Visible','off');
        set(handles.channel_selection_grp,'Visible','off');
        
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

model_index = get(hObject,'Value');
switch(model_index)
    case 1
        set(handles.disp_eq_txtfld,'String','y = (gain*x + offset)^gamma');
    case 2
        set(handles.disp_eq_txtfld,'String','y = (gain*x + offset1)^gamma + offset2');
%     case 3 
        set(handles.disp_eq_txtfld,'String','y = (gain*x)^gamma + offset');
    case 4
        set(handles.disp_eq_txtfld,'String','y = x^gamma');
end

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
A = reshape(C{1,1},3,numel(C{1,1})/3)';
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


% --- Executes on button press in save_image_bttn.
function save_image_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to save_image_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles.imData = imData;
% handles.filename = filename;
%         handles.pathname = pathname;
%         handles.filepath = filepath;
        
        [~, filename, ext] = fileparts(handles.filepath);
        %newname = sprintf('%s\%s_calib%s',pathstr,filename,ext);
        %newname = sprintf('%s%s',handles.filename,calib_file_suffix);
        
        %selected_channel = get(get(channel_selection_grp,'SelectedObject'),'String');
        %newname = sprintf('%s%s',filename,selected_channel);
        newname = sprintf('%s%s',filename,'_post');
        image_filepath = strrep(handles.filepath,filename,newname);
        
        switch handles.operation_mode
            case 0 %Display RGB
                switch handles.vis_draw_space
                    case 0
                        [fn,pathname,filterindex] = uiputfile(...
                            {'*.bmp; *.jpg; *.png', 'Selected channel (*.bmp, *.jpg, *.png)';
                            '*.bmp; *.jpg; *.png', 'Display calibrated image (*.bmp, *.jpg, *.png)';
                            '*.mat', 'XYZ data file (*.mat)'},...
                            'Save calibrated image as...',image_filepath);
                    case 1
                        image_filepath = strrep(image_filepath,ext,'.mat');
                        [fn,pathname,filterindex] = uiputfile(...
                            {'*.mat', 'Selected channel (*.mat)';
                            '*.bmp; *.jpg; *.png', 'Display calibrated image (*.bmp, *.jpg, *.png)';
                            '*.mat', 'XYZ data file (*.mat)'},...
                            'Save calibrated image as...',image_filepath);
                end
            case 1 %Scene XYZ
                image_filepath = strrep(image_filepath,ext,'.mat');
                [fn,pathname,filterindex] = uiputfile(...
                    {'*.mat', 'Selected channel (*.mat)';
                    '*.mat', 'XYZ data file (*.mat)'},...
                    'Save calibrated image as...',image_filepath);
        end
        if (fn ~= 0)            
            image_filepath = sprintf('%s%s',pathname,fn);
            Id = get_displayed_image(handles);
            switch handles.operation_mode
                case 0
                    switch filterindex
                        case 1
                            switch handles.sel_channel
                                case 4
                                    imwrite(Id,image_filepath);
                                otherwise
                                    save(image_filepath,'Id');
                            end
                        case 2
                            imwrite(Id,image_filepath);
                        case 3
                            Id = handles.ImxyzData; 
                            save(image_filepath,'Id');
                    end
                case 1
                    if filterindex == 2
                        Id = handles.ImxyzData; 
                    end
                    save(image_filepath,'Id');
            end

            status_msg = sprintf('[Image Tools]->[Image Calibration (Batch)] %s...calibrated',image_filepath);
            display(sprintf('%s',status_msg));
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
    



function resampling_factor_fld_Callback(hObject, eventdata, handles)
% hObject    handle to resampling_factor_fld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of resampling_factor_fld as text
%        str2double(get(hObject,'String')) returns contents of resampling_factor_fld as a double


% --- Executes during object creation, after setting all properties.
function resampling_factor_fld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to resampling_factor_fld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in auto_luminance_bttn.
function auto_luminance_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to auto_luminance_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
% TODO: If only the XYZ image has been computed, and not the RGB, this
% function will give an error because it works over data from Irgb_disp_in,
% which is the non-computed RGB image. 
%

    if isfield(handles,'imComputedData')
		Irgb_disp_in = handles.imComputedData;
        
        %tic
        
        a=reshape(Irgb_disp_in(:,:,1)>=255,[],1);
        b=reshape(Irgb_disp_in(:,:,2)>=255,[],1);
        c=reshape(Irgb_disp_in(:,:,3)>=255,[],1);
        saturation=([sum(a) sum(b) sum(c)]./numel(a)).*100;        
        prev_sat = sum(saturation>=0.3);
        
		prev_max = max(max(Irgb_disp_in));
		prev_max = prev_max(:,:,2);
        
		if (prev_sat >= 1)
            incr_sign = -1;
            incr = 0.5;
        else
            incr_sign = 1; 
            incr = 2;
        end
        
        if ~isfield(handles,'imComputedData')
            handles.resampling_factor = str2num(get(handles.resampling_factor_fld,'String'));
        end
        factor = handles.resampling_factor + (incr * incr_sign);
        
		while(incr > 0.1)
			Ixyz_disp_out_resampled = handles.ImxyzData_orig .* factor;

			% Convert display XYZ values to RGB using the conversion matrix.
			Irgb_disp_out = disp_xyz2rgb (Ixyz_disp_out_resampled, handles.disp_matrix);
			Irgb_disp_out(find(Irgb_disp_out<0)) = 0;     

			% Reverse the display's gamma using the display's fitting model.
			Irgb_disp_in = display_color_transform (Irgb_disp_out, handles.display_model);
			Irgb_disp_in(find(Irgb_disp_in<0)) = 0;
			Irgb_disp_in = uint8(Irgb_disp_in);

            
            a=reshape(Irgb_disp_in(:,:,1)>=255,[],1);
            b=reshape(Irgb_disp_in(:,:,2)>=255,[],1);
            c=reshape(Irgb_disp_in(:,:,3)>=255,[],1);
            saturation=([sum(a) sum(b) sum(c)]./numel(a)).*100;        
            total_sat = sum(saturation>=0.3);
            
			maxRGB = max(max(max(Irgb_disp_in)));

			%if ((incr_sign > 0) && (maxRGB >= 255) && (prev_max < 255)) || ((incr_sign < 0) && (maxRGB < 255) && (prev_max >= 255))
            if ((incr_sign > 0) && (total_sat >= 1) && (prev_sat < 1)) || ((incr_sign < 0) && (total_sat < 1) && (prev_sat >= 1))
                factor = factor - (incr * incr_sign);
				incr = incr / 2;
            else
                prev_sat = total_sat;
                prev_max = maxRGB;
			end
			
			factor = factor + (incr * incr_sign);
		end

        set(handles.resampling_factor_fld,'String',mat2str(factor));
        
        %toc
    else 
        errordlg('You must calibrate first','Calibration needed');
        return;
    end


% --- Executes on button press in radiance_bttn.
function radiance_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to radiance_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    % CHECK THAT AN INPUT IMAGE HAS BEEN ENTERED AND GIVE AN ERROR OTHERWISE.
    if ~isfield(handles,'imData')
       errordlg('Please select an input image fisrt','Input image not selected');
       return;
    end
    
    % If the XYZ image has not been already computed via another
    % module, do it here.
    %if ~isfield(handles,'ImxyzData_orig')
        % Convert RGB values to XYZ using the camera transfer function,
        % and make all the neccesary modifications to the image
        % according to the options selected in the GUI.
        [res, Ixyz_disp_out_orig, Ixyz_disp_out_resampled] = obtain_cam_xyz_image(handles.imData, hObject, handles);
        if (res==-1)
            errordlg('The camera matrix has not been set','Missing camera matrix');
            return;
        else
            % Save the computed data as a global variables.
            handles.ImxyzData_orig = Ixyz_disp_out_orig;
            handles.ImxyzData = Ixyz_disp_out_resampled;
            handles.imComputedData = Ixyz_disp_out_resampled;
            handles.operation_mode = 1; % 0 is for display RGB images, and 1 for scene XYZ.
            handles.vis_draw_space = 1; % 0 is for RGB color space, and 1 for XYZ.
            handles.sel_channel = 1;
            
            % Commit changes in global variables so they can be available
            % to other functions and modules.
            guidata(hObject, handles);
        end
    %end
    
    % Display the calibrated image in the corresponding axes. 
    axes(handles.computed_image_axes);
    handles.sel_channel = 1;
    imPlot = imshow(uint8(handles.ImxyzData(:,:,handles.sel_channel)),[]);
    handles.imPlot = imPlot;
    guidata(hObject, handles);

    gui_post_image_displayed (handles);
    
function Ixyz_d = displayableXYZ (Ixyz)
    vmax = max(max(max(Ixyz)));
    vmin = min(min(min(Ixyz)));
    
    Ixyz_d = ((Ixyz-vmin)/(vmax-vmin)).*255;


% --- Executes when selected object is changed in channel_selection_grp.
function channel_selection_grp_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in channel_selection_grp 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

channel = get(hObject,'String');
if (strcmp(channel,'X') || strcmp(channel,'R')), handles.sel_channel = 1;
else
    if (strcmp(channel,'Y') || strcmp(channel,'G')), handles.sel_channel = 2;
    else
        if (strcmp(channel,'Z') || strcmp(channel,'B')), handles.sel_channel = 3;
        else
            if (strcmp(channel,'All')), handles.sel_channel = 4;
            end
        end
    end
end

handles.imPlot = update_displayed_image (handles);
guidata(hObject, handles);

% --- Executes when selected object is changed in color_space_grp.
function color_space_grp_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in color_space_grp 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

space = get(hObject,'String');
if (strcmp(space,'XYZ'))
    handles.vis_draw_space = 1;
    handles.sel_channel = 1;
    set(handles.channel1_chkbx,'String','X');
    set(handles.channel2_chkbx,'String','Y');
    set(handles.channel3_chkbx,'String','Z');
    set(handles.channel4_chkbx,'Visible','off');
    set(handles.channel_selection_grp,'selectedobject',handles.channel1_chkbx)
else
    if (strcmp(space,'RGB'))
        handles.vis_draw_space = 0;
        handles.sel_channel = 4;
        set(handles.channel1_chkbx,'String','R');
        set(handles.channel2_chkbx,'String','G');
        set(handles.channel3_chkbx,'String','B');
        set(handles.channel4_chkbx,'Visible','on');
        set(handles.channel_selection_grp,'selectedobject',handles.channel4_chkbx)
    end
end

handles.imPlot = update_displayed_image (handles);
guidata(hObject, handles);


function Id = get_displayed_image(handles)
    switch handles.operation_mode
        case 0 %Display RGB
            switch handles.vis_draw_space
                case 0
                    if (handles.sel_channel == 4)
                        Id = handles.imComputedData (:,:,:);
                    else
                        Id = handles.imComputedData (:,:,handles.sel_channel);
                    end
                case 1
                    Id = handles.ImxyzData (:,:,handles.sel_channel);
            end
        case 1 %Scene XYZ
            Id = handles.ImxyzData (:,:,handles.sel_channel);
    end


function imPlot = update_displayed_image (handles)
   Id = get_displayed_image(handles);
   axes(handles.computed_image_axes);
   imPlot = imshow(uint8(Id),[]);
