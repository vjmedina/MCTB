
%# This class stores a color stimulus and all its associated properties 
%# into one single object, reducing the amount of variables in the 
%# workspace and grouping all the color information for a single stimulus. 
%# At the same time, it allows to perform different typical operations with 
%# color data, such as color space conversions, chromatic adaptations or 
%# color analysis.

classdef ColorData
    properties (Access = private)
        % Current set of Space-independent tristimulus values in use (in the XYZ colorspace)
        TV_XYZ        
        % To avoid unneccessary conversion, we use variable TV_input to
        % keep the original input tristimulus untouched, so we can return
        % them when the object is queried for them, instead of converting 
        % to XYZ and then back to the input space with the consequent lose
        % of precision. Of course, this only applies to non-spectral 
        % colors.
        TV_input        
        % Input Color Space. To be used with variable TV_input. 
        CS_input        
        % Flag to indicate if the tristimulus values have already been calculated.
        TV_set=0;        
        % Space-independent color coordinates
        coords
        % Current color space index
        ColorSpace
        % Original Color Matching Functions (so we can upsample without re-reading the file)
        CMF
        % Downsampled Color Matching Functions (needed by several methods)
        CMF_sampled
        % Illuminant used for the color of non-emissive surfaces.
        illum
        % XYZ luminance normalization for the reference white. Typically
        % this will be 1 (for values in the range [0.0, 1.0]) or 100 (for
        % values in the range [0, 100]).
        whiteLum
        % Normalization: Indicates whether the tristimulus values  should
        % be normalized (=1) or not (=0). If enabled, the value of whiteLum
        % is used to normalize. If no normalization is done, the stimulus
        % will keep its original range. 
        norm
    end
    properties (Constant)
        MAX_WAVELENGTH = 780;
        MIN_WAVELENGTH = 380;
        DEFAULT_SAMPLING_STEP = 5;
    end
    properties
        % Color spectrum (arranged column-wise)
        spd
        % Start value of the wavelength band spectrum
        spd_start
        % End value of the wavelength band spectrum
        spd_end
        % Step size (delta) of the wavelength band spectrum
        spd_delta
        % Identifier used to refer to the color (optional)
        id
        % Number that represents whether the color belongs to a reflective(0), transmissive (1) or emissive (2) surface.
        colorType
        % Flag to indicate if the color stimulus is created via spectral
        % or tristimulus values.
        spectral = ColorUtils.COLOR_SPECTRAL;
        % Standard observer used to obtain the Color Matching Functions (2
        % Degrees by default).
        observer = ColorUtils.STD_OBS_CIE_1931_2;
    end
    properties (Dependent)
        csShortName
        csLongName
    end
    properties (Constant)        
        DEFAULT_ILLUMINANT = ColorUtils.CIE_D65_xy2;
    end
    methods(Access = public)
        %% ColorData -- Class constructor  
        function obj = ColorData (mode, id, spectral)
            % 
            % Creates a ColorData object.
            %
            % INPUT:
            %   mode --> Type of surface the color belongs to. Allowed
            %            values are: 
            %                 ColorUtils.COLOR_SURFACE_REFLECTIVE
            %                 ColorUtils.COLOR_SURFACE_TRANSMISSIVE
            %                 ColorUtils.COLOR_SURFACE_EMISSIVE (Default)
            %   id  --> Id to refer to the color (optional)
            %   spectral--> indicates if the color stimulus contains
            %               spectral data or the tristimulus values are 
            %               given instead. (optional: it can be set later
            %               on with function "setSpectral"). Allowed values
            %               are:
            %                 ColorUtils.COLOR_SPECTRAL
            %                 ColorUtils.COLOR_NON_SPECTRAL
            %
            % OUTPUT:
            %   obj --> The newly created object.
            %
            
            if (nargin < 1)
                mode = ColorUtils.COLOR_SURFACE_EMISSIVE;
            end

            if (mode ~= ColorUtils.COLOR_SURFACE_REFLECTIVE) && (mode ~= ColorUtils.COLOR_SURFACE_TRANSMISSIVE) && (mode ~= ColorUtils.COLOR_SURFACE_EMISSIVE)
                error('Invalid mode: %d (Use only ColorUtils.COLOR_SURFACE_REFLECTIVE, ColorUtils.COLOR_SURFACE_TRANSMISSIVE or ColorUtils.COLOR_SURFACE_EMISSIVE).', mode);
            end
            
            obj.colorType = mode;
            
            if (nargin > 1)
                obj.id = id; 
            end
            
            if (nargin > 2)
                obj.spectral = spectral;
                if (spectral ~= ColorUtils.COLOR_SPECTRAL) && (spectral ~= ColorUtils.COLOR_NON_SPECTRAL)
                    error('Invalid argument: spectral can only be ColorUtils.COLOR_SPECTRAL or ColorUtils.COLOR_NON_SPECTRAL');
                end
            end
            
            % By default normalization is enabled to a luminance of 100.
            obj = obj.setNormalize(100);
            obj = obj.setWhiteLuminance(100);
            
            % By default, we load the CMF for 2-degree standard observer.
            obj.setObserver(ColorUtils.STD_OBS_CIE_1931_2);
            obj.CMF = ColorUtils.loadCMF(obj.getObserver());
        end
        
        %% setNormalize -- Enables/Disabled XYZ Tristimulus normalization
        function obj = setNormalize(obj, value)
            obj.norm = value;
        end
        
        %% getNormalize -- Checks if tristimulus normalization is enabled.
        function value = getNormalize(obj)
            value = obj.norm;
        end
        
        %% setSpectral -- Indicates wheter or not the color stimulus is spectral.
        function obj = setSpectral (obj, spectral)
            % 
            % Sets the type of surface the color belongs to.
            %
            % INPUT:
            %      obj --> The initial object before the changes.
            % spectral --> Type of surface the color belongs to. Allowed
            %            values are: 
            %                 COLOR_SPECTRAL
            %                 COLOR_NON_SPECTRAL
            %
            % OUTPUT:
            %   obj --> The modified object.
            %
            
            obj.spectral = spectral;            
        end
        
        %% setColorType -- Sets the type of surface the color belongs to.
        function obj = setColorType (obj, mode)
            % 
            % Sets the type of surface the color belongs to.
            %
            % INPUT:
            %    obj --> The initial object before the changes.
            %   mode --> Type of surface the color belongs to. Allowed
            %            values are: 
            %                 COLOR_SURFACE_REFLECTIVE
            %                 COLOR_SURFACE_TRANSMISSIVE
            %                 COLOR_SURFACE_EMISSIVE
            %
            % OUTPUT:
            %   obj --> The modified object.
            %
            
            obj.colorType = mode;            
        end
        
        %% setIlluminant -- Sets the spd of the illuminant.
        function obj = setIlluminant(obj, illum, spd_start, spd_end, spd_delta)
            % 
            % Stores the illuminant's SPD into the data object.
            %
            % INPUT:
            %         obj --> The initial object before the changes.
            %       illum --> The Nx1 or 1xN vector with the illuminant spd.
            %   spd_start --> Start value of the wavelength band spectrum.
            %     spd_end --> End value of the wavelength band spectrum.
            %   spd_delta --> Step size (delta) of the wavelength band spectrum.
            %
            % OUTPUT:
            %         obj --> The modified object.
            %
            
            if (isempty(obj.spd))
                obj.spd_start = spd_start;
                obj.spd_end = spd_end;
                obj.spd_delta = spd_delta;

                obj.illum = illum;
            else
                if (obj.spd_start == spd_start && obj.spd_end == spd_end && obj.spd_delta == spd_delta)
                    obj.illum = illum;
                else
                    if (obj.spd_start <= spd_start && obj.spd_end >= spd_end)
                       wl_sampling_old = spd_start:spd_delta:spd_end;
                       wl_sampling_new = obj.getSampling();
                       %wl_sampling_new = obj.spd_start:obj.spd_delta:obj.spd_end;

                       obj.illum = ColorUtils.resample_array([wl_sampling_old'; illum], wl_sampling_new);
                    else
                       error('The sampling interval is shorter than the object''s spd. There are some missing wavelengths.');
                    end
                end
            end
        end
        
        %% getId -- Returns the id of the color object.
        function id = getId(obj)
            % 
            % Sets the type of surface the color belongs to.
            %
            % INPUT:
            %         obj --> The object to read the information from.
            %
            % OUTPUT:
            %          id --> The color's id.
            %
            
            id = obj.id;
        end
        
        %% getSampling -- Returns the wavelength sampling interval
        function sampling = getSampling(obj)
            if (obj.spectral == ColorUtils.COLOR_NON_SPECTRAL)
                sampling = ColorData.MIN_WAVELENGTH:ColorData.DEFAULT_SAMPLING_STEP:ColorData.MAX_WAVELENGTH;                
            else
                sampling = (obj.spd_start:obj.spd_delta:obj.spd_end);
            end
        end

        %% getIlluminant -- Returns the illuminant's spd.
        function il = getIlluminant(obj)
            % 
            % Reads the spd of the color's illuminant.
            % Reflective and Transmissive surfaces require an illuminant to
            % be seen, otherwise they look black. This is not the case for
            % Emissive surfaces, which can be considered as illuminants in
            % themselves. 
            % To avoid errors when no illuminant is provided, we will use 
            % a black illuminant to simulate the effect on reflective and 
            % transmissive surfaces. For emissive surfaces, we will use the
            % stimulus' spd as the illuminant's spd in the case of spectral
            % stimuli, and an empty illuminant otherwise.
            %
            % INPUT:
            %         obj --> The object to read the information from.
            %
            % OUTPUT:
            %          il --> The spd of the color's illuminant.
            %
                
            if (isempty(obj.illum))
                if (obj.colorType == ColorUtils.COLOR_SURFACE_EMISSIVE)
                    if (obj.spectral == ColorUtils.COLOR_SPECTRAL)
                        % If the color is spectral, the spd already
                        % contains the luminous power of the illumination
                        % and therefore we can convert to lumens later on. 
                        sampling = obj.getSampling();
                        il = ones(numel(sampling),2);
                    end
                else
                    il = zeros(size(obj.spd));
                    warning('The object does not have any associated illuminant. Zero spectrum illuminant returned.');
                end
            else
                il = obj.illum;
            end
        end
        
        %% getIlluminantXYZ -- Returns the XYZ tristimulus values of the color's illuminant.
        function illumXYZ = getIlluminantXYZ(obj)
            % 
            % Returns the XYZ tristimulus values of the color's illuminant.
            % This is indeed equivalent to returning the XYZ tristimulus of
            % a perfect white diffuser placed in front of the illuminant.
            %
            % INPUT:
            %      obj --> The object to read the information from.
            %
            % OUTPUT:
            %   illumXYZ --> The XYZ tristimulus values of the color's illuminant.
            %
            
            il = obj.getIlluminant();
            
            %sampling_range = (obj.spd_start:obj.spd_delta:obj.spd_end)';
            sampling_range = obj.getSampling()';
            P = [sampling_range il];
            I = ones(size(P));
            I(:,1) = P(:,1);
            illumXYZ = ColorUtils.spect2XYZ(P, obj.getCMF(), obj.spd_start, obj.spd_end, obj.spd_delta, I);
            
            % Normalize luminance.
            if (obj.getNormalize() == 1)
                illumXYZ = illumXYZ .* (obj.getWhiteLuminance()/illumXYZ(2));
            end
        end        
        
        %% setSpectrum -- Sets the color signal's spd.
        function obj = setSpectrum (obj, spd, spd_start, spd_end, delta)
            % 
            % Sets the color signal's spd.
            %
            % INPUT:
            %         obj --> The object to read the information from.
            %         spd --> The color signal's spd. It must be a Nx1
            %                 column vector where N is the number of
            %                 wavelengths.
            %   spd_start --> The initial wavelength of the spd band.
            %     spd_end --> The final wavelength of the spd band.
            %       delta --> The step size between sampled wavelengths.            
            %
            % OUTPUT:
            %         obj --> The modified object.
            %
            
            if (obj.spectral == ColorUtils.COLOR_SPECTRAL)
                
                [sx, sy] = size(spd);
                if (sy > sx), spd = spd'; end

                obj.spd = spd;
                obj.spd_start = spd_start;
                obj.spd_end = spd_end;
                obj.spd_delta = delta;            

                obj = obj.setSampling(spd_start, spd_end, delta);
                
                obj = obj.computeTV();
            else
                warning('Cannot assign a spectrum to a non-spectral color. METHOD CALL IGNORED.');
            end
        end 
        
        %% setTVs -- Sets the tristimulus values directly, without using an spd.
        function obj = setTVs(obj, colorSpace, ABC, M_W)
            % 
            % Sets the tristimulus values directly, without using an spd..
            %
            % INPUT:
            %         obj --> The object to read the information from.
            %  colorSpace --> The color space in which the values are
            %                 expressed.
            %         ABC --> The object to read the information from.
            %         M_W --> Dual argument. Its value depends on the colorSpace:
            %                  COLOR_SPACE_RGB: The conversion matrix for the color space's primaries.
            %                  COLOR_SPACE_LAB: The XYZ tristimulus values of the reference white.
            %                  COLOR_SPACE_LUV: The XYZ tristimulus values of the reference white.
            %
            % OUTPUT:
            %    illumXYZ --> The XYZ tristimulus values of the color's illuminant.
            %
            
            if (obj.spectral == ColorUtils.COLOR_NON_SPECTRAL)
                if colorSpace < 1 || colorSpace > numel(ColorUtils.spacenames_long)
                    error('The chosen color space does not exist. ''ColorData.listColorSpaces'' for more info.');                        
                end

                obj.CS_input = colorSpace;
                obj.TV_input = ABC;
                
                switch colorSpace
                    case ColorUtils.COLOR_SPACE_XYZ                        
                        obj = obj.setXYZ(ABC);
                        
                    case ColorUtils.COLOR_SPACE_xyY
                        obj = obj.setXYZ(ColorUtils.xyY2XYZ(ABC));
                        
                    case ColorUtils.COLOR_SPACE_LAB
                        if (nargin < 4)
                            error('A reference white is required for the Lab color space.');
                        end
                        obj = obj.setXYZ(ColorUtils.Lab2XYZ(ABC, M_W));
                        
                    case ColorUtils.COLOR_SPACE_RGB
                        sampling = obj.getSampling();
                        M = ColorProfile.compute_rgb2xyz_matrix_illum ([sampling' obj.getIlluminant()]);
                        if (nargin < 4)
                            %obj = obj.setXYZ(ColorUtils.RGB2XYZ(ABC));
                            obj = obj.setXYZ(ColorUtils.RGB2XYZ(ABC,M));
                        else
                            obj = obj.setXYZ(ColorUtils.RGB2XYZ(ABC, M_W));
                            %obj = obj.setXYZ(ColorUtils.RGB2XYZ(ABC, M));
                        end
                        
                    case ColorUtils.COLOR_SPACE_sRGB
                        sampling = obj.getSampling();
                        M = ColorProfile.compute_rgb2xyz_matrix_illum ([sampling' obj.getIlluminant()]);
                        obj = obj.setXYZ(ColorUtils.sRGB2XYZ(ABC, M));
                        
                    case ColorUtils.COLOR_SPACE_LUV
                        obj = obj.setXYZ(ColorUtils.Luv2XYZ(ABC, M_W));
                        
                    otherwise
                        error('Conversion method not yet implemented for the chosen color space.');                    
                end
            else
                warning('Cannot assign tristimulus values directly to a spectral color. You must first add a spectrum with setSpectrum() and then obtain the corresponding tristimulus values with getTVs().');
            end
            
        end
        
        %%
        function obj = setXYZ(obj, XYZ)            
            obj.TV_XYZ = XYZ;
            obj.TV_set = 1;
        end
        
        %%
        function obj = setSampling(obj, newstart, newend, newdelta)
            A = [(obj.spd_start:obj.spd_delta:obj.spd_end)' obj.spd];
            temp = ColorUtils.resample_array (A, newstart:newdelta:newend);
            
            obj.spd_start = newstart;
            obj.spd_end = newend;
            obj.spd_delta = newdelta;
            obj.spd = temp(:,2);
            
            %obj.CMF_sampled = zeros(size(obj.spd,1),size(obj.CMF_sampled,2));
            %for i=2:4                
            %    temp = ColorUtils.resample_array ([obj.CMF(:,1) obj.CMF(:,i)], newstart:newdelta:newend);
            %    obj.CMF_sampled(:,[1,i]) = temp;
            %end
            
        end
       
         %%
        function obj = setLuminance(obj, lum)
        % Set the color stimulus' luminance to a specific value. As opposed
        % to setWhiteLuminance, this function sets the stimulus' luminance
        % to the exact value specified as parameter. 
        %
        % INPUT:
        %   obj --> The ColorData object to be modified.
        %   lum --> The exact new luminance value of the estimulus.
        %
        % OUTPUT:
        %   obj --> The modified object.
        %   
            if (obj.TV_set==1)  
                TV = obj.getTVs(ColorUtils.COLOR_SPACE_XYZ);
                k = (lum/TV(2));

                obj = obj.setTVs(ColorUtils.COLOR_SPACE_XYZ, (TV .* k));
                obj.whiteLum = obj.whiteLum * k;
            else
                warning('The tristimulus values of the color have not been computed yet.');
            end
        end
        
        %%
        function obj = setWhiteLuminance(obj, lum)
        % Set the color stimulus' luminance according to its reference white. This
        % operation consists in modifying all the XYZ tristimulus values so
        % that the luminance of the reference white is exactly "lum".
        %
        % INPUT:
        %   obj --> The ColorDate object to be modified.
        %   lum --> The new luminance value of the reference white.
        %
        % OUTPUT:
        %   obj --> The modified object.
        %
                        
            if (obj.TV_set==1)
                TV = obj.getTVs(ColorUtils.COLOR_SPACE_XYZ);
                %cmf = obj.getCMF();
                %I = obj.getIlluminant();

                % Compute the normalizing constant.
                %k = lum/(sum( cmf(:,3) .* I(:) ));
                
                obj = obj.setTVs(ColorUtils.COLOR_SPACE_XYZ, (TV .* (lum/obj.whiteLum)));
            end
            obj.whiteLum = lum;
        end
        
        %%
        function lum = getLuminance (obj)
            TV = obj.getTVs(ColorUtils.COLOR_SPACE_XYZ);
            lum = TV(2);
        end
        
        %%
        function k = getWhiteLuminance(obj)
            %W_XYZ = obj.getIlluminantXYZ();
            %k = W_XYZ(2);
            k = obj.whiteLum;
        end

        %%
        function obs = getObserver(obj)
            obs = obj.observer;
        end
        
        %%
        function obj = setObserver(obj, obs)
            obj.observer = obs;
        end
        
        %%
        function TV = getTVs(obj, colorSpace, clipping, M_W)
            if nargin < 3
                clipping = 1;
            end
            
            if colorSpace < 1 || colorSpace > numel(ColorUtils.spacenames_long)
                error('The chosen color space does not exist. ''ColorData.listColorSpaces'' for more info.');                    
            end
            
            if (obj.spectral == ColorUtils.COLOR_NON_SPECTRAL) && (obj.CS_input == colorSpace)
                TV = obj.TV_input;
            else            
                switch colorSpace
                    case ColorUtils.COLOR_SPACE_XYZ
                        TV = obj.toXYZ();
                    case ColorUtils.COLOR_SPACE_xyY
                        TV = obj.toxyY();
                    case ColorUtils.COLOR_SPACE_LAB
                        TV = obj.toLab(M_W);
                    case ColorUtils.COLOR_SPACE_RGB
                        sampling = obj.getSampling();
                        ill = obj.getIlluminant();
                        if ~isempty(ill)
                            M = ColorProfile.compute_rgb2xyz_matrix_illum ([sampling' ill]);
                            M = inv(M);
                        else
                            if (nargin < 4)
                                warning('Missing XYZ2RGB conversion matrix. A standard one will be used');
                                M = ColorUtils.CIE_XYZ2RGB_MATRIX;
                            else
                                M = M_W;
                            end
                        end
                        TV = obj.toRGB(M);
                        
                        if clipping == 1
                            TV( TV < 0 ) = 0; % Set all negative values to 0.
                        end 
                    case ColorUtils.COLOR_SPACE_sRGB
                        sampling = obj.getSampling();
                        M = ColorProfile.compute_rgb2xyz_matrix_illum ([sampling' obj.getIlluminant()]);
                        M = inv(M);
                        TV = obj.tosRGB(M);
                    case ColorUtils.COLOR_SPACE_LUV
                        TV = obj.toLuv(M_W);
                    case ColorUtils.COLOR_SPACE_LMS
                        TV = obj.toLMS(M_W);
                    otherwise
                        error('Conversion method not yet implemented for the chosen color space.');                    
                end
            end
                         
        end
        
        %%
        function obj = scaleXYZ (obj, factor)
            obj.TV_XYZ = obj.TV_XYZ ./ factor;
        end
        
        %%
        function ccord = colorCoord(obj, colorSpace)
            if nargin < 2
                colorSpace = ColorUtils.COLOR_SPACE_XYZ;
            end

            if (obj.TV_set == 0)
                obj = obj.computeTV();
            end
            
            ccord(1) = obj.TV_XYZ(1)/sum(obj.TV_XYZ);
            ccord(2) = obj.TV_XYZ(2)/sum(obj.TV_XYZ);
            
            switch(colorSpace)
                case ColorUtils.COLOR_SPACE_LUV
                    ccord = ColorUtils.xy2uv(ccord);
                otherwise
                    if (colorSpace ~= ColorUtils.COLOR_SPACE_XYZ)
                        warning('The calculation of chromaticity coordinates for this color space is not implemented. xy coordinates returned instead.');
                    end
            end
        end

                
        %% getCMF -- Returns the resampled standard Color matching functions.
        function cmfs = getCMF(obj, deg)
            if (nargin < 2)                
                deg = obj.getObserver();
            end
            
            switch (deg)
                case ColorUtils.STD_OBS_CIE_1931_2
                    cmfs = ColorUtils.CIE_1931_XYZ_CMF_1nm';
                case ColorUtils.STD_OBS_CIE_1964_10
                    cmfs = ColorUtils.CIE_1964_XYZ10_CMF_1nm';
            end
            
            sampling = obj.getSampling();
            cmfs = ColorUtils.resample_array(cmfs,sampling);
        end
        
        %% plot -- Plot a point at the color object's coordinates.
        function h = plot(obj, colorSpace, varargin)
            if (nargin < 3)
                cwidth = 80;          
                formatArgs = {cwidth};
            else
                labelpoint = false;
                nVarargs = length(varargin);
                i=1;
                while i <= nVarargs
                    param = varargin{i};
                    if ischar(param)
                        switch(param)
%                             case 'CircleWidth'
%                                 cwidth = varargin{i+1};
                            case 'CCM'
                                % The Color Convertion matrix is only used for RGB.
                                M = varargin{i+1};
                                tv = obj.getTVs(colorSpace, 1, M);
%                             case 'White'
%                                 Wr = varargin{i+1};
%                                 tv = obj.getTVs(colorSpace, 1, Wr);
                            case 'label'
                                label = varargin{i+1};
                                labelpoint=true;
                            case 'format'
                                formatArgs = {varargin{i+1:nVarargs}};
                                break;
                        end
                    end
                    i=i+1;
                end
            end

            cc = obj.colorCoord(colorSpace);
            h = scatter(cc(1),cc(2),formatArgs{:});
            
            if (labelpoint)
                text(cc(1)-0.03, cc(2)+0.02, cellstr(label),'FontSize',15);
            end
        end
        
        %% isEqualTo -- Checks if two colors are equal.
        function eq = isEqualTo(obj, C, JND, ill_xy)
            if nargin < 3
                JND = 2.3; % By default, we consider two color as dfiferent if differ by more than 3% of its maximun value.
            end
            
            if nargin < 4
                ill_xy = ColorUtils.CIE_D65_xy2;
            end
           
            ill_XYZ = ColorUtils.xyY2XYZ ([ill_xy obj.getWhiteLuminance()]);            
            eq = 0;
            
            lab1 = obj.getTVs(ColorUtils.COLOR_SPACE_LAB, 0, ill_XYZ);
            lab2 = C.getTVs(ColorUtils.COLOR_SPACE_LAB, 0, ill_XYZ);
            
            dE00 = ColorUtils.difference_CIEDE2000 (lab1, lab2);            

            if ( dE00 < JND )
                eq = 1;                    
            end
        end
        
        %% colorDiffLuv -- CIELuv distance between two colors.
        function dE = colorDiffLuv(obj, C, ill_xy)
            if nargin < 3                
                ill_xyY = ColorUtils.XYZ2xyY(obj.getIlluminantXYZ());
                ill_xy = ill_xyY(1:2);
            end
           
            ill_XYZ = ColorUtils.xyY2XYZ ([ill_xy obj.getWhiteLuminance()]);
                        
            luv1 = obj.getTVs(ColorUtils.COLOR_SPACE_LUV, 0, ill_XYZ);
            luv2 = C.getTVs(ColorUtils.COLOR_SPACE_LUV, 0, ill_XYZ);
            
            dE = ColorUtils.Euclidean_distance (luv1, luv2);
        end
        
        %% colorDiff76 -- CIE76 distance between two colors.
        function dE = colorDiff76(obj, C, ill_xy)           
            if nargin < 3
                ill_xy = ColorUtils.CIE_D65_xy2;
            end
           
            ill_XYZ = ColorUtils.xyY2XYZ ([ill_xy obj.getWhiteLuminance()]);            
            %eq = 0;
            
            lab1 = obj.getTVs(ColorUtils.COLOR_SPACE_LAB, 0, ill_XYZ);
            lab2 = C.getTVs(ColorUtils.COLOR_SPACE_LAB, 0, ill_XYZ);
            
            dE = ColorUtils.Euclidean_distance (lab1, lab2);
        end
        
        %% colorDiff94 -- CIE94 distance between two colors.
        function dE = colorDiff94(obj, C, ill_xy)
            %if nargin < 3
            %    JND = 2.3; % By default, we consider two color as dfiferent if differ by more than 3% of its maximun value.
            %end
            
            if nargin < 3
                ill_xy = ColorUtils.CIE_D65_xy2;
            end
           
            ill_XYZ = ColorUtils.xyY2XYZ ([ill_xy obj.getWhiteLuminance()]);            
            %eq = 0;
            
            lab1 = obj.getTVs(ColorUtils.COLOR_SPACE_LAB, 0, ill_XYZ);
            lab2 = C.getTVs(ColorUtils.COLOR_SPACE_LAB, 0, ill_XYZ);
            
            dE = ColorUtils.difference_CIE94 (lab1, lab2);

            %if ( dE00 < JND )
            %    eq = 1;                    
            %end            
        end
        
        %% colorDiff00 -- CIE2000 distance between two colors.
        function dE = colorDiff00(obj, C, ill_xy)
            %if nargin < 3
            %    JND = 2.3; % By default, we consider two color as dfiferent if differ by more than 3% of its maximun value.
            %end
            
            if nargin < 3
                ill_xy = ColorUtils.CIE_D65_xy2;
            end
           
            ill_XYZ = ColorUtils.xyY2XYZ ([ill_xy obj.getWhiteLuminance()]);            
            %eq = 0;
            
            lab1 = obj.getTVs(ColorUtils.COLOR_SPACE_LAB, 0, ill_XYZ);
            lab2 = C.getTVs(ColorUtils.COLOR_SPACE_LAB, 0, ill_XYZ);
            
            dE = ColorUtils.difference_CIEDE2000 (lab1, lab2);

            %if ( dE00 < JND )
            %    eq = 1;                    
            %end            
        end
    end
    methods(Static)
        %% loadIlluminant -- Reads the object's illuminant into a variable.
        function Ill = loadIlluminant (illum, wl_start, wl_end, wl_step, CMFfile)
            if (nargin < 5)
                CMFfile='.\Illuminants.xlsx';
                if (exist(CMFfile, 'file') ~= 2)
                    error('Illuminants.xlsx was not found in the current folder');
                end
            end
            
            switch (illum)
                case ColorUtils.ILLUM_D55
                    if (wl_step ~= 5)
                        error('Wrong step size. Use only 5 nm.');
                    end
                    Ill = xlsread(CMFfile,'Illuminant D55 5nm');
                case ColorUtils.ILLUM_D65
                    switch wl_step
                        case 1
                            Ill = xlsread(CMFfile,'Illuminant D65 1nm');
                        case 5
                            Ill = xlsread(CMFfile,'Illuminant D65 5nm');
                        case 10
                            Ill = xlsread(CMFfile,'Illuminant D65 10nm');
                        otherwise
                            error('Wrong step size. Use only 1, 5 or 10 nm.')
                    end
                    
                case ColorUtils.ILLUM_A
                    if (wl_step ~= 5)
                        error('Wrong step size. Use only 5 nm.');
                    end
                    Ill = xlsread(CMFfile,'Illuminant A 5nm');
                case ColorUtils.ILLUM_F7
                    if (wl_step ~= 5)
                        error('Wrong step size. Use only 5 nm.');
                    end
                    Ill = xlsread(CMFfile,'Illuminant F7 5nm');                        
                otherwise
            end
            
            start_row = find(Ill(:,1)==wl_start);
            end_row = find(Ill(:,1)==wl_end);
            
            if (isempty(start_row) || isempty(end_row))
                error('Wrong wavelength range.');
            end
            
            Ill = Ill(start_row:end_row,:);            
        end
        
        %% listColorSpaces -- Lists all available color spaces in the class.
        function listColorSpaces
            str = sprintf('\nColor Spaces:\Variable\t\tLong name');
            str = sprintf('%s\n-------------------------------------------',str);                
            str = sprintf('%sCOLOR_SPACE_XYZ\t\t\t\t%s',str,char(ColorUtils.spacenames_long(ColorUtils.COLOR_SPACE_XYZ)));
            str = sprintf('%sCOLOR_SPACE_xyY\t\t\t\t%s',str,char(ColorUtils.spacenames_long(ColorUtils.COLOR_SPACE_xyY)));
            str = sprintf('%sCOLOR_SPACE_LAB\t\t\t\t%s',str,char(ColorUtils.spacenames_long(ColorUtils.COLOR_SPACE_LAB)));
            str = sprintf('%sCOLOR_SPACE_LUV\t\t\t\t%s',str,char(ColorUtils.spacenames_long(ColorUtils.COLOR_SPACE_LUV)));
            str = sprintf('%sCOLOR_SPACE_RGB\t\t\t\t%s',str,char(ColorUtils.spacenames_long(ColorUtils.COLOR_SPACE_RGB)));
            str = sprintf('%sCOLOR_SPACE_sRGB\t\t\t\t%s',str,char(ColorUtils.spacenames_long(ColorUtils.COLOR_SPACE_sRGB)));
            str = sprintf('%sCOLOR_SPACE_YIQ\t\t\t\t%s',str,char(ColorUtils.spacenames_long(ColorUtils.COLOR_SPACE_YIQ)));
            str = sprintf('%sCOLOR_SPACE_YUV\t\t\t\t%s',str,char(ColorUtils.spacenames_long(ColorUtils.COLOR_SPACE_YUV)));
            str = sprintf('%sCOLOR_SPACE_YCbCr\t\t\t\t%s',str,char(ColorUtils.spacenames_long(ColorUtils.COLOR_SPACE_YCbCr)));
            str = sprintf('%sCOLOR_SPACE_HSV\t\t\t\t%s',str,char(ColorUtils.spacenames_long(ColorUtils.COLOR_SPACE_HSV)));
            str = sprintf('%sCOLOR_SPACE_HSL\t\t\t\t%s',str,char(ColorUtils.spacenames_long(ColorUtils.COLOR_SPACE_HSL)));
            str = sprintf('%sCOLOR_SPACE_CMYK\t\t\t\t%s',str,char(ColorUtils.spacenames_long(ColorUtils.COLOR_SPACE_CMYK)));
            str = sprintf('%sCOLOR_SPACE_LMS\t\t\t\t%s',str,char(ColorUtils.spacenames_long(ColorUtils.COLOR_SPACE_LMS)));

            str = sprintf('%s\n',str);
            disp(str);
        end
    end
    methods(Access = private)        
        
        %%
        function obj = computeTV(obj)
            sampling_range = obj.getSampling()';
            
            if (obj.colorType == ColorUtils.COLOR_SURFACE_EMISSIVE)
                il = ones(size(obj.spd));
            else
                il = obj.getIlluminant();
            end
            cmf = obj.getCMF();
            
            if (obj.getNormalize()==1)                
                obj.TV_XYZ = ColorUtils.spect2XYZ([sampling_range obj.spd], cmf, obj.spd_start, obj.spd_end, obj.spd_delta, [sampling_range il], obj.getWhiteLuminance());
            else
                obj.TV_XYZ = ColorUtils.spect2XYZ([sampling_range obj.spd], cmf, obj.spd_start, obj.spd_end, obj.spd_delta, [sampling_range il]);
            end
            obj.TV_set = 1;
        end
        
        %%
        function XYZ = toXYZ(obj)
            if (obj.TV_set == 1)
                XYZ = obj.TV_XYZ;
            else                
                error('The object does not contain any color information yet.');
            end
        end
        
        %%
        function xyY = toxyY(obj)            
            if (obj.TV_set == 0)
                obj = obj.computeTV();
            end
            
            xyY(1) = obj.TV_XYZ(1)/sum(obj.TV_XYZ);
            xyY(2) = obj.TV_XYZ(2)/sum(obj.TV_XYZ);
            xyY(3) = obj.TV_XYZ(2);
        end
        
        %%
        function RGB = toRGB(obj, M)
            if (obj.TV_set == 0)
                obj = obj.computeTV();
            end
                        
            if (nargin < 2)                
                RGB = ColorUtils.XYZ2RGB (obj.TV_XYZ);
            else
                RGB = ColorUtils.XYZ2RGB (obj.TV_XYZ, M);
            end
        end
        
        %%
        function Lab = toLab(obj,W_XYZ)
            %W_XYZ = W.getTVs(ColorUtils.COLOR_SPACE_XYZ,0);
            XYZ = obj.getTVs(ColorUtils.COLOR_SPACE_XYZ,0);
            Lab = ColorUtils.XYZ2Lab(XYZ,W_XYZ);
        end
        
        %%
        function Luv = toLuv(obj,W_XYZ)
            %W_XYZ = W.getTVs(ColorUtils.COLOR_SPACE_XYZ,0);
            XYZ = obj.getTVs(ColorUtils.COLOR_SPACE_XYZ,0);
            Luv = ColorUtils.XYZ2Luv(XYZ,W_XYZ);
        end
        
        %%
        function LMS = toLMS(obj,W_XYZ)
            %W_XYZ = W.getTVs(ColorUtils.COLOR_SPACE_XYZ,0);
            XYZ = obj.getTVs(ColorUtils.COLOR_SPACE_XYZ,0);
            LMS = ColorUtils.XYZ2LMS(XYZ,W_XYZ);
        end
        
        %%
        function sRGB = tosRGB(obj, M)
            XYZ = obj.getTVs(ColorUtils.COLOR_SPACE_XYZ,0);
            sRGB = ColorUtils.XYZ2sRGB(XYZ, M);
        end 
       
    end
end
        
        
        