function varargout = camera_calibration(varargin)
% CAMERA_CALIBRATION MATLAB code for camera_calibration.fig
%      CAMERA_CALIBRATION, by itself, creates a new CAMERA_CALIBRATION or raises the existing
%      singleton*.
%
%      H = CAMERA_CALIBRATION returns the handle to a new CAMERA_CALIBRATION or the handle to
%      the existing singleton*.
%
%      CAMERA_CALIBRATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CAMERA_CALIBRATION.M with the given input arguments.
%
%      CAMERA_CALIBRATION('Property','Value',...) creates a new CAMERA_CALIBRATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before camera_calibration_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to camera_calibration_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help camera_calibration

% Last Modified by GUIDE v2.5 06-Jan-2016 12:40:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @camera_calibration_OpeningFcn, ...
                   'gui_OutputFcn',  @camera_calibration_OutputFcn, ...
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


% --- Executes just before camera_calibration is made visible.
function camera_calibration_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to camera_calibration (see VARARGIN)

set(handles.m_rgb2xyz_table, 'Data', zeros(3,4));
set(handles.m_xyz2rgb_table, 'Data', zeros(3,4));

set(handles.corr_rgb2xyz_table, 'Data', zeros(1,3));
set(handles.corr_xyz2rgb_table, 'Data', zeros(1,3));

% Choose default command line output for camera_calibration
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes camera_calibration wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = camera_calibration_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in rgb_sample_list.
function rgb_sample_list_Callback(hObject, eventdata, handles)
% hObject    handle to rgb_sample_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns rgb_sample_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from rgb_sample_list


% --- Executes during object creation, after setting all properties.
function rgb_sample_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rgb_sample_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in slctSmpl_btn.
function slctSmpl_btn_Callback(hObject, eventdata, handles)
% hObject    handle to slctSmpl_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if isfield(handles,'imSelection')
        handles = rmfield(handles,'imSelection');
    end
    
    imSelection = imrect(get(handles.imPlot,'Parent'));
    handles.imSelection = imSelection;
    guidata(hObject, handles);

% --- Executes on button press in rgb_delete_bttn.
function rgb_delete_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to rgb_delete_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    selected_indices = get(handles.rgb_sample_list,'Value');
    initial_values = cellstr(get(handles.rgb_sample_list,'String'));
    list_size = numel(initial_values);
    new_list = [];
    
    %if (list_size >= 1)
    if (numel(initial_values{1,1})~=0)
        for current_index = 1:list_size
            if isempty(find (selected_indices == current_index, 1))
                new_list = [new_list(:); initial_values(current_index)];
            end
        end
        
        set(handles.rgb_sample_list,'String',new_list);
        if (selected_indices(1)==list_size)
            set(handles.rgb_sample_list,'Value',numel(new_list));
        else
            set(handles.rgb_sample_list,'Value',selected_indices(1));
        end
    end
    
    
    

% --- Executes on button press in calibrate_btn.
function calibrate_btn_Callback(hObject, eventdata, handles)
% hObject    handle to calibrate_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    % Read the RGB values
    
    %rgb_values = cellstr(get(handles.rgb_sample_list,'String'));
    rgb_values = get(handles.rgb_sample_list,'String');
    list_size = size(rgb_values,1);
        
    % If there are no elements in the list, quit the function.
    if (list_size == 0 || size(rgb_values,2)== 0)
        errordlg('No RGB samples selected','Wrong number of samples');
        return; 
    end
    RGB = zeros(list_size,3);
    
    % Read the colors rgb values
    for i=1:list_size
        temp = textscan(char(rgb_values(i)),'%f %f %f');
        RGB(i,:)=[temp{1,1} temp{1,2} temp{1,3}];
    end
    
    % Read the XYZ values
     XYZ = get(handles.xyz_uitable,'Data');
     
    %if (size(XYZ,1) == 0)
    % If the table is not of type double, it means that it is empty or corrupted
    if (~isa(XYZ, 'double'))
        errordlg('Missing radiances','Wrong data');
        return; 
    end
    
    % CALIBRATE
    
    % Check that they have the same number of samples
    if (size(RGB,1) == size(XYZ,1))
        %RGB_exp = RGB;
        %XYZ_exp = XYZ;

        % Check if over and/or under exposure thresholds are given and
        % store their values.
        if (get(handles.overexp_threshold_chkbx,'Value')==1)
            overexp_th = str2double(get(handles.overexp_threshold_txtfld,'String'));
        else
            overexp_th = max(max(RGB));
        end

        if (get(handles.underexp_threshold_chkbx,'Value')==1)
            underexp_th = str2double(get(handles.underexp_threshold_txtfld,'String'));
        else
            underexp_th = min(min(min(RGB)),0);
        end

        % REMOVE over and/or under exposed samples from the array. 
        
        index = true(1, size(RGB, 1));
        
        % find samples with one or more channels greater/smaller than
        % each threshold.
        gt = mod(find(RGB > overexp_th)-1,size(RGB,1))+1;
        lt = mod(find(RGB < underexp_th)-1,size(RGB,1))+1;
        a = unique([gt;lt]);
        index(a) = false;
                
        % Remove those entries from both XYZ and RGB to maintain the same
        % samples in both arrays.
        
        RGB = RGB(index, :);
        XYZ = XYZ(index, :);
        
        
        
        RGB_aug = [RGB, ones(size(RGB,1),1)];
        XYZ_aug = [XYZ, ones(size(XYZ,1),1)];

        M_rgb2xyz = ((RGB_aug' * RGB_aug) \ (RGB_aug' * XYZ))'; %inv([RGB|1]'*[RGB|1]) * ([RGB|1]' * XYZ)
        M_xyz2rgb = ((XYZ_aug' * XYZ_aug) \ (XYZ_aug' * RGB))';
        
        %XYZ_calib = (RGB_exp * M_rgb2xyz(:,1:3)') + repmat(M_rgb2xyz(:,4)',size(RGB_exp,1),1);        
        %XYZ_corr = [corrcoef(XYZ(:,1), XYZ_calib(:,1)) corrcoef(XYZ(:,2), XYZ_calib(:,2)) corrcoef(XYZ(:,3), XYZ_calib(:,3))];
        %XYZ_corr = (XYZ_corr(1,[2 4 6]));
        XYZ_calib = (RGB_aug * M_rgb2xyz');
        XYZ_corr = diag(corr(XYZ, XYZ_calib));
        
        %RGB_calib = (XYZ * M_xyz2rgb(:,1:3)') + repmat(M_xyz2rgb(:,4)',size(XYZ,1),1);        
        %RGB_corr = [corrcoef(RGB(:,1), RGB_calib(:,1)) corrcoef(RGB(:,2), RGB_calib(:,2)) corrcoef(RGB(:,3), RGB_calib(:,3))];
        %RGB_corr = (RGB_corr(1,[2 4 6]));
        RGB_calib = (XYZ_aug * M_xyz2rgb');
        RGB_corr = diag(corr(RGB, RGB_calib));
        
        set(handles.corr_rgb2xyz_table,'Data',XYZ_corr');
        set(handles.corr_xyz2rgb_table,'Data',RGB_corr');
        
        set(handles.m_rgb2xyz_table,'Data',M_rgb2xyz);
        set(handles.m_xyz2rgb_table,'Data',M_xyz2rgb);
        
        % This section is disabled because there is currently no way to
        % know the color of each sample. We could suppose that the user
        % will always input data from a 24-sample Macbeth color chart, but 
        % then we would produce false data if that is not the case. 
        % Moreover, since we may remove some entries when a threshold is 
        % used, this would alter the position of each color sample, which 
        % would again result in false data being produced. 
        
%         if size(RGB,1)>=24
%             rg_ratio = mean(RGB(19:24,1)./RGB(19:24,2));
%             bg_ratio = mean(RGB(19:24,3)./RGB(19:24,2));
%         else
%             rg_ratio = 0;
%             bg_ratio = 0;
%         end
%         
%         set(handles.rg_ratio_txtfld,'String',rg_ratio);
%         set(handles.bg_ratio_txtfld,'String',bg_ratio);
%         
%         RGB_th = [2^16*rg_ratio, 2^16, 2^16*bg_ratio];
%         XYZ_th = [RGB_th, 1] * M_rgb2xyz';
%         
%         set(handles.theo_wr_txtfld,'String',sprintf('%.4f',RGB_th(1)));
%         set(handles.theo_wg_txtfld,'String',sprintf('%.4f',RGB_th(2)));
%         set(handles.theo_wb_txtfld,'String',sprintf('%.4f',RGB_th(3)));
%         
%         set(handles.theo_wx_txtfld,'String',sprintf('%.4f',XYZ_th(1)));
%         set(handles.theo_wy_txtfld,'String',sprintf('%.4f',XYZ_th(2)));
%         set(handles.theo_wz_txtfld,'String',sprintf('%.4f',XYZ_th(3)));
        
%         set(handles.real_wx1_txtfld,'String',sprintf('%.4f',XYZ(19,1)));
%         set(handles.real_wy1_txtfld,'String',sprintf('%.4f',XYZ(19,2)));
%         set(handles.real_wz1_txtfld,'String',sprintf('%.4f',XYZ(19,3)));
%         
%         set(handles.real_wx2_txtfld,'String',sprintf('%.4f',XYZ(19,1)./sum(XYZ(19,:))));
%         set(handles.real_wy2_txtfld,'String',sprintf('%.4f',XYZ(19,2)./sum(XYZ(19,:))));
%         set(handles.real_wz2_txtfld,'String',sprintf('%.4f',XYZ(19,3)./sum(XYZ(19,:))));
        
        set(handles.M1_copy_bttn,'Enable','on');
        set(handles.M2_copy_bttn,'Enable','on');
                        
        plotCalibrationGraph(handles, RGB, RGB_calib, 'RGB calibration', {'Rcalib','Gcalib','Bcalib'});
        plotCalibrationGraph(handles, XYZ, XYZ_calib, 'XYZ calibration', {'Xcalib','Ycalib','Zcalib'});
    else
        errordlg('The number of XYZ and RGB samples does not match','Wrong number of samples');
        return;
    end
    

% --- Executes on button press in browse_btn.
function browse_btn_Callback(hObject, eventdata, handles)
% hObject    handle to browse_btn (see GCBO)
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

        status_msg = sprintf('[Image Tools]->[Camera Calibration] %s...loaded',filepath);    
        display(sprintf('%s',status_msg));
        %status_msg = get(handles.results_txtfld,'string');
        %new_msg = sprintf('%s...loaded',filepath);    
        %status_msg = sprintf('%s\n%s',status_msg,new_msg);    
        %set(handles.results_txtfld,'string',status_msg);

        %imshow(filepath, 'Parent', handles.img_viewer);
        warning('off','images:initSize:adjustingMag');
        imFigure = figure('name',filename);
        imData = imread(handles.filepath);        
        imPlot = imshow(imData);
        handles.imFigure = imFigure;
        handles.imData = imData;
        handles.imPlot = imPlot;
        guidata(hObject, handles);
    end
    

% --- Executes on button press in strSmpl_btn.
function strSmpl_btn_Callback(hObject, eventdata, handles)
% hObject    handle to strSmpl_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if isfield(handles,'imSelection')
        pos = getPosition(handles.imSelection);    
        
        x = ceil(pos(1));
        y = ceil(pos(2));
        w = ceil(pos(3));
        h = ceil(pos(4));

        % The rectangle coordinates are given in screen coordinates, whereas
        % the image pixels must be addressed in image coordinates (inversed),
        % which is why we exchange the x and y coordinates.
        
        res = mean(mean(handles.imData( y:y+h, x:x+w, : ) ));
        handles.last_sampled_value = res;      
        
        initial_values = cellstr(get(handles.rgb_sample_list,'String'));
        if cellfun(@isempty,initial_values)
            new_value = [sprintf('%.2f  %.2f  %.2f',res(1),res(2),res(3))];
        else
            new_value = [initial_values; sprintf('%.2f  %.2f  %.2f',res(1),res(2),res(3))];
        end
        set(handles.rgb_sample_list,'String',new_value)
        set(handles.rgb_sample_list,'Value',size(new_value,1));
        guidata(hObject, handles);
        
    else
        uiwait(msgbox('Please select pixels first','Error','modal'));
    end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isfield(handles,'imFigure')    % Check if the window has been created.
    if ishandle(handles.imFigure)  % Check is the window remains open (exists)
        close(handles.imFigure);
    end
end
delete(hObject);

% --- Executes on button press in radiance_paste_bttn.
function radiance_paste_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to radiance_paste_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str = clipboard('paste');
C = textscan(str,'%f');
A = reshape(C{1,1},3, numel(C{1,1})/3)';

oldData = get(handles.xyz_uitable,'Data');

if (isnumeric(oldData))
    newData = [oldData; A];
else
    newData = A;
end

set(handles.xyz_uitable,'Data',newData);

    
function cs = storeColors(RGB)

    [sx, sy] = size(RGB);

    cs = ColorSet();
    
    for i = 1:sx        
        c = ColorData(ColorUtils.COLOR_SURFACE_EMISSIVE, i, ColorUtils.COLOR_NON_SPECTRAL);
        %c = c.setIlluminant (il, spd_start, spd_end, delta);
        c = c.setTVs(ColorUtils.COLOR_SPACE_RGB, RGB(i,:));
        cs = cs.addColor(c);        
    end

function [] = plotCalibrationGraph(handles, Data_exp, Data_calib, chart_title, chart_legend)
    if nargin < 4
        chart_legend='';
    end

    if nargin < 3
        chart_title='Camera calibration chart';        
    end    
    
    cam = get(handles.cam_model_txtfld,'String');
    focal = get(handles.cam_focal_txtfld,'String');
    iso = get(handles.cam_iso_txtfld,'String');
    aperture = get(handles.cam_apert_txtfld,'String');
    shutter = get(handles.cam_shutter_txtfld,'String');
    
    CamInfoStr = sprintf('%s (F%s - %s ISO - %s mm - %s sec.)',cam,aperture,iso,focal,shutter);
    
    figure,
    hold on,
    
    title (sprintf('%s\n%s',chart_title,CamInfoStr));
    
    scatter(Data_exp(:,1), Data_calib(:,1),'r','s','MarkerFaceColor','r');
    scatter(Data_exp(:,2), Data_calib(:,2),'g','d','MarkerFaceColor','g');
    scatter(Data_exp(:,3), Data_calib(:,3),'b','^','MarkerFaceColor','b');
    
    if (~isempty(chart_legend))
        legend(chart_legend,'Location','Southeast');
    end
    
    xlabel('Measured data');
    ylabel('Estimated data');
    
    range = max(max(max(Data_exp,Data_calib)));
    axis([0 range 0 range]);
    
    grid on,
    
    hold off,

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over calibrate_btn.
function calibrate_btn_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to calibrate_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in radiance_clear_bttn.
function radiance_clear_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to radiance_clear_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.xyz_uitable,'Data',{'','',''});


% --- Executes on button press in xyz_output_bttn.
function xyz_output_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to xyz_output_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listData = get(handles.xyz_uitable,'Data');
set(handles.output_txtfld,'String',listData);


% --- Executes on button press in rgb_copy_bttn.
function rgb_copy_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to rgb_copy_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%index_selected = get(handles.rgb_sample_list,'Value');
%selected_value = file_list{index_selected};
copied_values = get(handles.rgb_sample_list,'String');

sel_length =  numel(copied_values);
RGB = zeros(sel_length,3);

% Split each selected entry into the three RGB components. 
for i=1:sel_length
    temp = textscan(char(copied_values(i)),'%f %f %f');
    RGB(i,:)=[temp{1,1} temp{1,2} temp{1,3}];
end
num2clip(RGB);
h = msgbox('The selected values were successfully copied to the clipboard','Operation Completed');


% --- Executes on button press in rgb_paste_bttn.
function rgb_paste_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to rgb_paste_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Copy data from the clipboard.
data = paste();
if ~isa(data,'double')
    errordlg('The clipboard does not contain numbers','Wrong type of data');
    return; 
end

% Convert "data" from an array of doubles to a cell array of strings.
tempStr =arrayfun(@num2str, data, 'unif', 0);

% Merge all three columns (RGB) into the same colum. This is because the
% listbox requires a 1-column cell array of strings where each row is
% stored in a separate line in the list. 
tempStr = strcat(tempStr(:,1),{' '},tempStr(:,2),{' '},tempStr(:,3));

% Store the data in the listbox.
set(handles.rgb_sample_list,'String',tempStr);
set(handles.rgb_sample_list,'Value',1);



function cam_model_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to cam_model_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cam_model_txtfld as text
%        str2double(get(hObject,'String')) returns contents of cam_model_txtfld as a double


% --- Executes during object creation, after setting all properties.
function cam_model_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cam_model_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cam_focal_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to cam_focal_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cam_focal_txtfld as text
%        str2double(get(hObject,'String')) returns contents of cam_focal_txtfld as a double


% --- Executes during object creation, after setting all properties.
function cam_focal_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cam_focal_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cam_iso_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to cam_iso_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cam_iso_txtfld as text
%        str2double(get(hObject,'String')) returns contents of cam_iso_txtfld as a double


% --- Executes during object creation, after setting all properties.
function cam_iso_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cam_iso_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cam_apert_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to cam_apert_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cam_apert_txtfld as text
%        str2double(get(hObject,'String')) returns contents of cam_apert_txtfld as a double


% --- Executes during object creation, after setting all properties.
function cam_apert_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cam_apert_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cam_shutter_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to cam_shutter_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cam_shutter_txtfld as text
%        str2double(get(hObject,'String')) returns contents of cam_shutter_txtfld as a double


% --- Executes during object creation, after setting all properties.
function cam_shutter_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cam_shutter_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in M1_copy_bttn.
function M1_copy_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to M1_copy_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

copied_values = get(handles.m_rgb2xyz_table,'Data');
num2clip(copied_values);
h = msgbox('The selected values were successfully copied to the clipboard','Operation Completed');


% --- Executes on button press in M2_copy_bttn.
function M2_copy_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to M2_copy_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

copied_values = get(handles.m_xyz2rgb_table,'Data');
num2clip(copied_values);
h = msgbox('The selected values were successfully copied to the clipboard','Operation Completed');


% --- Executes on button press in M1_corr_copy_bttn.
function M1_corr_copy_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to M1_corr_copy_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

copied_values = get(handles.corr_rgb2xyz_table,'Data');
num2clip(copied_values);
h = msgbox('The selected values were successfully copied to the clipboard','Operation Completed');

% --- Executes on button press in M2_corr_copy_bttn.
function M2_corr_copy_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to M2_corr_copy_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

copied_values = get(handles.corr_xyz2rgb_table,'Data');
num2clip(copied_values);
h = msgbox('The selected values were successfully copied to the clipboard','Operation Completed');



function rg_ratio_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to rg_ratio_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rg_ratio_txtfld as text
%        str2double(get(hObject,'String')) returns contents of rg_ratio_txtfld as a double


% --- Executes during object creation, after setting all properties.
function rg_ratio_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rg_ratio_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bg_ratio_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to bg_ratio_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bg_ratio_txtfld as text
%        str2double(get(hObject,'String')) returns contents of bg_ratio_txtfld as a double


% --- Executes during object creation, after setting all properties.
function bg_ratio_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bg_ratio_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function theory_XYZ_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to theory_XYZ_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of theory_XYZ_txtfld as text
%        str2double(get(hObject,'String')) returns contents of theory_XYZ_txtfld as a double


% --- Executes during object creation, after setting all properties.
function theory_XYZ_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to theory_XYZ_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function theory_RGB_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to theory_RGB_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of theory_RGB_txtfld as text
%        str2double(get(hObject,'String')) returns contents of theory_RGB_txtfld as a double


% --- Executes during object creation, after setting all properties.
function theory_RGB_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to theory_RGB_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function real_wx1_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to real_wx1_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of real_wx1_txtfld as text
%        str2double(get(hObject,'String')) returns contents of real_wx1_txtfld as a double


% --- Executes during object creation, after setting all properties.
function real_wx1_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to real_wx1_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function real_wx2_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to real_wx2_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of real_wx2_txtfld as text
%        str2double(get(hObject,'String')) returns contents of real_wx2_txtfld as a double


% --- Executes during object creation, after setting all properties.
function real_wx2_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to real_wx2_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function real_wy1_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to real_wy1_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of real_wy1_txtfld as text
%        str2double(get(hObject,'String')) returns contents of real_wy1_txtfld as a double


% --- Executes during object creation, after setting all properties.
function real_wy1_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to real_wy1_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function real_wz1_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to real_wz1_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of real_wz1_txtfld as text
%        str2double(get(hObject,'String')) returns contents of real_wz1_txtfld as a double


% --- Executes during object creation, after setting all properties.
function real_wz1_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to real_wz1_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function real_wy2_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to real_wy2_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of real_wy2_txtfld as text
%        str2double(get(hObject,'String')) returns contents of real_wy2_txtfld as a double


% --- Executes during object creation, after setting all properties.
function real_wy2_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to real_wy2_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function real_wz2_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to real_wz2_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of real_wz2_txtfld as text
%        str2double(get(hObject,'String')) returns contents of real_wz2_txtfld as a double


% --- Executes during object creation, after setting all properties.
function real_wz2_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to real_wz2_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function theo_wr_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to theo_wr_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of theo_wr_txtfld as text
%        str2double(get(hObject,'String')) returns contents of theo_wr_txtfld as a double


% --- Executes during object creation, after setting all properties.
function theo_wr_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to theo_wr_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function theo_wx_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to theo_wx_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of theo_wx_txtfld as text
%        str2double(get(hObject,'String')) returns contents of theo_wx_txtfld as a double


% --- Executes during object creation, after setting all properties.
function theo_wx_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to theo_wx_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function theo_wg_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to theo_wg_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of theo_wg_txtfld as text
%        str2double(get(hObject,'String')) returns contents of theo_wg_txtfld as a double


% --- Executes during object creation, after setting all properties.
function theo_wg_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to theo_wg_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function theo_wb_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to theo_wb_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of theo_wb_txtfld as text
%        str2double(get(hObject,'String')) returns contents of theo_wb_txtfld as a double


% --- Executes during object creation, after setting all properties.
function theo_wb_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to theo_wb_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function theo_wy_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to theo_wy_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of theo_wy_txtfld as text
%        str2double(get(hObject,'String')) returns contents of theo_wy_txtfld as a double


% --- Executes during object creation, after setting all properties.
function theo_wy_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to theo_wy_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function theo_wz_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to theo_wz_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of theo_wz_txtfld as text
%        str2double(get(hObject,'String')) returns contents of theo_wz_txtfld as a double


% --- Executes during object creation, after setting all properties.
function theo_wz_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to theo_wz_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function overexp_threshold_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to overexp_threshold_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of overexp_threshold_txtfld as text
%        str2double(get(hObject,'String')) returns contents of overexp_threshold_txtfld as a double


% --- Executes during object creation, after setting all properties.
function overexp_threshold_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to overexp_threshold_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in overexp_threshold_chkbx.
function overexp_threshold_chkbx_Callback(hObject, eventdata, handles)
% hObject    handle to overexp_threshold_chkbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of overexp_threshold_chkbx
if get(hObject,'Value')==1
    set(handles.overexp_threshold_txtfld,'Enable','on');
else
    set(handles.overexp_threshold_txtfld,'Enable','off');
end



function underexp_threshold_txtfld_Callback(hObject, eventdata, handles)
% hObject    handle to underexp_threshold_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of underexp_threshold_txtfld as text
%        str2double(get(hObject,'String')) returns contents of underexp_threshold_txtfld as a double


% --- Executes during object creation, after setting all properties.
function underexp_threshold_txtfld_CreateFcn(hObject, eventdata, handles)
% hObject    handle to underexp_threshold_txtfld (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in underexp_threshold_chkbx.
function underexp_threshold_chkbx_Callback(hObject, eventdata, handles)
% hObject    handle to underexp_threshold_chkbx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of underexp_threshold_chkbx
if get(hObject,'Value')==1
    set(handles.underexp_threshold_txtfld,'Enable','on');
else
    set(handles.underexp_threshold_txtfld,'Enable','off');
end
