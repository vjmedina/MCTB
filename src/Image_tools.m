function varargout = Image_tools(varargin)
% IMAGE_TOOLS MATLAB code for Image_tools.fig
%      IMAGE_TOOLS, by itself, creates a new IMAGE_TOOLS or raises the existing
%      singleton*.
%
%      H = IMAGE_TOOLS returns the handle to a new IMAGE_TOOLS or the handle to
%      the existing singleton*.
%
%      IMAGE_TOOLS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGE_TOOLS.M with the given input arguments.
%
%      IMAGE_TOOLS('Property','Value',...) creates a new IMAGE_TOOLS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Image_tools_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Image_tools_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Image_tools

% Last Modified by GUIDE v2.5 13-May-2015 17:01:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Image_tools_OpeningFcn, ...
                   'gui_OutputFcn',  @Image_tools_OutputFcn, ...
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

addpath('../src', '../images','../temp','./lib','./CSToolBox');

% --- Executes just before Image_tools is made visible.
function Image_tools_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Image_tools (see VARARGIN)

% Choose default command line output for Image_tools
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Image_tools wait for user response (see UIRESUME)
% uiwait(handles.figure1);


set(handles.crtSelection_bttn,'Enable','off');
set(handles.dltSelection_bttn,'Enable','off');
set(handles.average_bttn,'Enable','off');
            
% --- Outputs from this function are returned to the command line.
function varargout = Image_tools_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in average_bttn.
function average_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to average_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    if isfield(handles,'imSelection')
        pos = getPosition(handles.imSelection);    
        %pos = [xmin ymin width height]
        x = ceil(pos(1));
        y = ceil(pos(2));
        w = ceil(pos(3));
        h = ceil(pos(4));
        
        
        % The rectangle coordinates are given in screen coordinates, whereas
        % the image pixels must be addressed in image coordinates (inversed),
        % which is why we exchange the x and y coordinates.
        
        res = mean(mean(handles.imData( y:y+h, x:x+w, : ) ));
        if (ndims(handles.imData)==3)            
            resultMsg = sprintf('%.5f \t %.5f \t %.5f',res(1),res(2),res(3));
        else
            resultMsg = sprintf('%.5f',res);
        end     
        
        %resultMsg = sprintf('Pixel selection average: \n\n Red: %.5f \n Green: %.5f \n Blue: %.5f\n',res(1),res(2),res(3));
        %resultMsg = sprintf('Pixel selection average: \n\n %.5f \n %.5f \n %.5f\n',res(1),res(2),res(3));

        disp(resultMsg);
        %uiwait(msgbox(resultMsg,'Average','modal'));
    else
        uiwait(msgbox('Please select pixels first','Error','modal'));
    end

% --- Executes on button press in browse_bttn.
function browse_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to browse_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [filename, pathname, filterindex] = uigetfile( ...
    {'*.bmp;*.jpg;*.jpeg;*.png;*.ppm;*.tif;*.tiff',...
     'All image files (*.bmp;*.jpg;*.jpeg;*.png;*.ppm;*.tiff)';
     '*.bmp',  'Bitmap images(*.bmp)'; ...
     '*.jpeg;*.jpg',  'JPEG images(*.jpeg, *.jpg)'; ...
     '*.png',  'Portable Network Graphics images(*.png)'; ...
     '*.tif;*.tiff',  'Tagged Image File Format images (*.tif, *.tiff)'; ...
     '*.mat',  'Matlab array(*.mat)'}, ...
     'Pick an image');

    if ~isnumeric(filename)
        filepath = strcat(pathname,filename);
        handles.filename = filename;
        handles.pathname = pathname;
        handles.filepath = filepath;

        status_msg = sprintf('[Image Tools] %s...loaded',filepath);    
        display(sprintf('%s',status_msg));
        %status_msg = get(handles.results_txtfld,'string');
        %new_msg = sprintf('%s...loaded',filepath);    
        %status_msg = sprintf('%s\n%s',status_msg,new_msg);    
        %set(handles.results_txtfld,'string',status_msg);

        %imshow(filepath, 'Parent', handles.img_viewer);
        imFigure = figure('name',filename);
        
        if (filterindex == 6)
            found = 0;
            data = load(handles.filepath);
            for f = fieldnames(data)'
               var = data.(f{1});
               dims = ndims(var);
               if (dims == 2 || dims == 3)
                    found = found + 1;
                    if (size(var,3) > 3)
                        imData = var(:,:,1:3);
                    else
                        imData = var;
                    end
               end
            end
            if (found == 0)
                errordlg('The file does not contain any accepted array.','Missing array');
                return;
            end
            if (found > 1)
                warndlg('There are multiple possible arrays in the file. Only one was read.','Multiple arrays');
            end
            
            imPlot = imshow(uint8(imData),[]);
        else
            imData = imread(handles.filepath);
            imPlot = imshow(imData);
        end
        
        set(handles.crtSelection_bttn,'Enable','on');
        set(handles.dltSelection_bttn,'Enable','on');
        set(handles.average_bttn,'Enable','on');
        
        handles.imFigure = imFigure;
        handles.imData = imData;
        handles.imPlot = imPlot;
        guidata(hObject, handles);
    end


% --- Executes on button press in camera_calibration_btn.
function camera_calibration_btn_Callback(hObject, eventdata, handles)
% hObject    handle to camera_calibration_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    %uiwait(camera_calibration);
    camera_calibration;
    
    
% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in crtSelection_bttn.
function crtSelection_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to crtSelection_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if ~isfield(handles,'imSelection')
        imSelection = imrect(get(handles.imPlot,'Parent'));
        handles.imSelection = imSelection;
        guidata(hObject, handles);
    end


% --- Executes on button press in dltSelection_bttn.
function dltSelection_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to dltSelection_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if isfield(handles,'imSelection')
        delete(handles.imSelection);
        handles = rmfield(handles, 'imSelection');
        guidata(hObject, handles);
    end
    
    
    
    


% --- Executes on button press in display_calibration_bttn.
function display_calibration_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to display_calibration_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    %uiwait(display_calibration);
    handles.win_display_calibration = display_calibration;
    guidata(hObject, handles);

% --- Executes on button press in image_calibration_bttn.
function image_calibration_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to image_calibration_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    %uiwait(image_calibration);
    handles.win_image_calibration = image_calibration;
    guidata(hObject, handles);
    


% --- Executes on button press in image_calib_batch_bttn.
function image_calib_batch_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to image_calib_batch_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    %uiwait(image_calibration_batch);
    handles.win_image_calibration_batch = image_calibration_batch;
    guidata(hObject, handles);


% --- Executes on button press in color_bttn.
function color_bttn_Callback(hObject, eventdata, handles)
% hObject    handle to color_bttn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.win_colorimetric_tests = colorimetric_tests;
guidata(hObject, handles);
