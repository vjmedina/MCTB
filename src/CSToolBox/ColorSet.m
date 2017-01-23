classdef ColorSet
%# COLORSET - This class defines an interface for working with color sets.
%#            It consists of a list of colors, defined as type ColorData,
%#            and a series of common operation with color sets.
%
%Author: Victor Medina Heierle
%Last update: 23-Jan-2017
%
%  
    properties
        Colors % List colors associated with the color set.
        num_colors
    end
    properties (SetAccess = private)
        calibrated = 0; % Indicates whether or not the colorset contains calibration information.
    end
    methods (Static)
        
        function XYZ = getPrimaryXYZ(color)
            range = 400:0.1:700;
            spd = zeros(size(range));
            mask = [ColorUtils.CIE_RED_WL ColorUtils.CIE_GREEN_WL ColorUtils.CIE_BLUE_WL];
            
            mask = mask .* color;
            for i=1:3, 
                if (mask(i)>0)
                    index = int8((mask(i)-400)/0.1)+1;
                    spd(index)=1; 
                end
            end
            
            CMF = ColorUtils.loadCMF();
            XYZ = ColorUtils.spect2XYZ([range; spd], CMF, 400, 700, 0.1);            
        end
        
        function locus = compute_cie_locus_points(wl_start, wl_end, step, colorSpace)            
            if nargin < 3, step=5; end
            wl_range = ((wl_end-wl_start)/step)+1;
                        
            % Load the color matching functions	
            %CMF = loadColorMatchingFunctions(varargin{:});
            CMF = ColorUtils.loadCMF();
                        
            [sx,sy] = size(CMF);
            
            if ~((sx == wl_range && sy == 3) || (sx >= wl_range && sy == 4)), error('Wrong number of rows and/or columns for the color matching functions.'); end
            
            if sy == 4
                if ( (CMF(1,1) > wl_start) || (CMF(sx,1) < wl_end) ), error('Wrong color matching functions range.' ); end
                if ( (CMF(1,1) < wl_start) || (CMF(sx,1) > wl_end) )
                    temp = zeros(wl_range,3);
                    for i=1:3                        
                        t = resample_array ([CMF(:,1) CMF(:,i+1)], wl_start:step:wl_end);
                        temp(:,i) = t(:,2);
                    end
                    CMF = temp;                
                end                
            end
            
            mono_colors = zeros(wl_range,wl_range);
            locus = zeros(wl_range,2);
            for i=1:wl_range
                mono_colors(i,i) = 1;
                XYZ = ColorUtils.spect2XYZ([(wl_start:step:wl_end); mono_colors(i,:)]', [(wl_start:step:wl_end)' CMF], wl_start, wl_end, step);
                switch(colorSpace)
                    case ColorUtils.COLOR_SPACE_XYZ
                        ccord = ColorUtils.XYZ2xy(XYZ);
                    case ColorUtils.COLOR_SPACE_LUV
                        ccord = ColorUtils.XYZ2uv(XYZ);
                end
                locus(i,1) = ccord(1);
                locus(i,2) = ccord(2);
            end
            
            dx=(locus(61,1)-locus(1,1))/60;
            dy=(locus(61,2)-locus(1,2))/60;

            x=[locus(61,1), locus(1,1)];
            y=[locus(61,2), locus(1,2)];

            purple_boundary = [x(1)-dx:-dx:x(2);y(1)-dy:-dy:y(2)]';

            % Join the initial and final points of the spectrum locus to draw the purple boundary
            for i=1:length(purple_boundary)
                locus(61+i,:) = [purple_boundary(i,1),purple_boundary(i,2)];
            end                        
        end
        function plotLocus(colorSpace)
            
            wl_start = ColorData.MIN_WAVELENGTH;
            wl_end = ColorData.MAX_WAVELENGTH;
            wl_step = ColorData.DEFAULT_SAMPLING_STEP;
            
            locus = ColorSet.compute_cie_locus_points(wl_start,wl_end, wl_step, colorSpace);
            
            for i=1:size(locus,1)
                %fill(locus(:,1), locus(:,2),bgcolor);
                plot(locus(:,1), locus(:,2));
            end
        end    
    end
    methods
        function obj = ColorSet(data, spd_start, spd_end, spd_delta)
            %
            % Create a new Color Set.
            %
            % If no arguments are given, the resulting object will be
            % initialized to an empty Color space object. Otherwise, we can
            % pass an existing data matrix to fill in the data from. 
            %
            % Note:
            %   If we don't want to use spectral colors, but instead just
            %   give the stimuli's tristimulus values, we must create an
            %   empty ColorSet object and add the color stimuli one by one
            %   with the addColor function.
            %
            %
            % INPUT:            
            %        data --> A MxN data matrix containing the spd of each
            %                 color element to be added to the set. The matrix must
            %                 have the following structure:
            %
            %                 | r_11 r_12 r_13...r_1n |
            %                 | r_21 r_22 r_23...r_2n |
            %                 | r_31 r_32 r_33...r_3n |
            %                 | .... .... ....   .... | 
            %                 | r_m1 r_m2 r_m3...r_mn |
            %
            %                 where M is the number of color samples and N is the
            %                 number of measured wavelengths.
            %
            %   spd_start --> The initial wavelength of the spd band.
            %     spd_end --> The final wavelength of the spd band.
            %   spd_delta --> The step size between sampled wavelengths.    
            %
            % OUTPUT:
            %       obj --> The resulting color set object.
            %
            %
            
            obj.num_colors = 0;
            obj.Colors = ColorData.empty(2,0);
            
            if (nargin ~= 0)
                obj = obj.loadColors(data, spd_start, spd_end, spd_delta);
            end            
        end
        
        function obj = loadColors(obj, data, spd_start, spd_end, spd_delta)
            % 
            % Read color samples from a new data set variable.
            %
            % The result is the same as creating a new Color set object
            % with values from a data variable. The parameters are the 
            % same as those for the ColorSet constructor.
            %
            
            [sx, sy] = size(data);
            
            for i=1:sx
                c = ColorData(ColorUtils.COLOR_SURFACE_EMISSIVE, i, ColorUtils.COLOR_SPECTRAL);
                c = c.setSpectrum (data(i,:), spd_start, spd_end, spd_delta);
                obj = obj.addColor(c);
            end
        end
        
        function obj = setSampling(obj, newstart, newend, newdelta)            
            for i=1:obj.num_colors
                obj.Colors(i) = obj.Colors(i).setSampling(newstart, newend, newdelta);
            end            
        end
        
        function obj = addColor(obj, C) 
            if (~isa(C,'ColorData')), error('The parameter must be an object of type ColorData.'); end            
            obj.Colors(obj.num_colors+1) = C;
            obj.num_colors = obj.num_colors + 1;
        end
        
        function h = plotColor(obj, idx, size, colorSpace)
            if (nargin < 2)
                colorSpace = ColorUtils.COLOR_SPACE_XYZ;
            end
            if (nargin < 3)
                size=50;
            end
            
            h = obj.Colors(idx).plot(colorSpace, 'CircleWidth',size);
        end             
        function h = plotAllColors(obj, colorSpace, varargin)
            if (nargin < 2)
                colorSpace = ColorUtils.COLOR_SPACE_XYZ;
            end
            h = zeros(obj.num_colors, 1);
            for i=1:obj.num_colors
                h(i) = obj.Colors(i).plot(colorSpace, varargin{:});
            end
        end 
        function h = plotChromaticityDiagram(obj, colorSpace, graphTitle)
            if (nargin < 3)
                colorSpace_name = ColorUtils.getColorSpaceName(colorSpace);
                graphTitle = sprintf('%s(%s) Diagram',char(colorSpace_name(1)),char(colorSpace_name(2)));
            end
            
            figure;            
            
            switch(colorSpace)
                case ColorUtils.COLOR_SPACE_XYZ
                    bgImg = imread('CIE1931xy.png');
                    x1=0.76; y1=0.88;
                    x2=0.9; y2=0.9;
                case ColorUtils.COLOR_SPACE_LUV
                    bgImg = imread('CIE1976_UCS2.png');
                    x1=0.62; y1=0.59;
                    x2=0.62; y2=0.60;
                otherwise
                    warning('Chromaticity diagram not implemented. Showing the xy chromaticity diagram instead.');
                    colorSpace = ColorUtils.COLOR_SPACE_XYZ;
                    colorSpace_name = ColorUtils.getColorSpaceName(colorSpace);
                    graphTitle = sprintf('%s(%s) Diagram',colorSpace_name(1),colorSpace_name(2));
                    bgImg = imread('CIE1931xy.png');
                    x1=0.76; y1=0.88;
                    x2=0.9; y2=0.9;
            end
            
            imagesc([0 x1], [0 y1], flipdim(bgImg,1));
            title(graphTitle);            
            hold on,
            set(gca,'ydir','normal'); 
            ColorSet.plotLocus(colorSpace);
            P_locus = ColorUtils.PlanckianLocus(colorSpace);
            plot(P_locus(:,1), P_locus(:,2), 'w', 'LineWidth', 2);
            h = obj.plotAllColors(colorSpace, 'format', 80, 'fill', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
            axis([0,x2,0,y2]),
            set(gca,'XTick',0:0.1125:x2),
            set(gca,'YTick',0:0.1125:y2),
            grid on,
            hold off,
        end        
        function l = getColorList(obj, colorSpace, normalize)
            if (nargin < 3)
                normalize = 0;
            end
            % Normalize all values so that the maximum is 1 and the minimum
            % is 0.
            if (normalize == 1)
                minV = obj.getMinValue(colorSpace);
                maxV = obj.getMaxValue(colorSpace) - minV;                
            end
            
            l = zeros (obj.num_colors, 4);
            for i=1:obj.num_colors
                l(i,4) = obj.Colors(i).getId();
                l(i,1:3) = obj.Colors(i).getTVs(colorSpace);
            end
            if (normalize == 1)
                l(:,1:3) = l(:,1:3) - minV;
                l(:,1:3) = l(:,1:3) ./ maxV;
            end            
        end
        function obj = sortById(obj)
            temp = obj;
            for i=1:obj.num_colors
                 id = temp.Colors(i).getId();
                 obj.Colors(i) = temp.Colors(id);
            end            
        end
        function maxL = getMaxLuminance(obj)
            maxL = 0;
            for i=1:obj.num_colors
                xyz = obj.Colors(i).getTVs(ColorUtils.COLOR_SPACE_XYZ);
                if (xyz(2) > maxL)
                    maxL = xyz(2);
                end
            end
        end
        function maxV = getMaxValue(obj, colorSpace)
            maxV = 0;
            for i=1:obj.num_colors
                xyz = obj.Colors(i).getTVs(colorSpace);
                if (max(xyz) > maxV)
                    maxV = max(xyz);
                end
            end
        end 
        function minV = getMinValue(obj, colorSpace)
            minV = 9999999999;
            for i=1:obj.num_colors
                xyz = obj.Colors(i).getTVs(colorSpace);
                if (min(xyz) < minV)
                    minV = min(xyz);
                end
            end
        end 
                
        function obj = normalizeLuminance(obj)
            % This function is meant to normalize the tristimulus of all the
            % colors in the set to a maximum luminance of 100. 
            %
            % NOTE: This method has not been tested so it might not behave as
            % expected.            
            %
            maxY = obj.getMaxLuminance();
            factor = maxY/100;
            for i=1:obj.num_colors                
                obj.Colors(i) = obj.Colors(i).scaleXYZ (factor);
            end
        end
        
        function [M, V] = compute_calib_matrix (obj)
            RGB_exp = obj.getColorList(ColorUtils.COLOR_SPACE_RGB);
            RGB_exp = RGB_exp(:,1:3);
            LMS_exp = zeros(size(RGB_exp));
            for i=1:size(RGB_exp,1)
                LMS_exp(i,:) = ColorUtils.XYZ2LMS(ColorUtils.RGB2XYZ(RGB_exp(i,:)));
            end
            [sx,sy]=size(RGB_exp);
            RGB_ext = [RGB_exp ones(sx,1)]; % extend rgb with a column of ones.
            M_V = (inv(RGB_ext' * RGB_ext)*(RGB_ext' * LMS_exp))'; 
            
            % M_V should be a 3x4 matrix where the last column is V.
            M = M_V(:,1:3);
            V = M_V(:,4);
        end
        
        function [obj, XYZ_calib, corr] = calibrate (obj)
            % Convert RGB to XYZ
            RGB_exp = obj.getColorList(ColorUtils.COLOR_SPACE_RGB);
            RGB_exp = RGB_exp(:,1:3);
            XYZ_exp = zeros(size(RGB_exp));            
            for i=1:size(RGB_exp,1)
                XYZ_exp(i,:) = ColorUtils.XYZ2LMS(ColorUtils.RGB2XYZ(RGB_exp(i,:)));
            end
            
            % Calibrate
            [M, V] = obj.compute_calib_matrix();
            XYZ_calib = M * RGB_exp + V;
            
            % Compute the correlation between the measured and calibrated
            % values.
            corr = zeros(1,3);
            for i=1:3
                tmp_corr = corrcoef(XYZ_exp(:,i),XYZ_calib(:,i));
                corr(i) = tmp_corr(1,2);
            end
            
            obj.calibrated = 1;
        end
        
        function h = plot_calib_graph (obj)
            % Convert RGB to XYZ
            RGB_exp = obj.getColorList(ColorUtils.COLOR_SPACE_RGB);
            RGB_exp = RGB_exp(:,1:3);
            XYZ_exp = zeros(size(RGB_exp));
            for i=1:size(RGB_exp,1)
                XYZ_exp(i,:) = ColorUtils.XYZ2LMS(ColorUtils.RGB2XYZ(RGB_exp(i,:)));
            end
            
            figure;
            
            %h = 
            
            
        end
    end
end
