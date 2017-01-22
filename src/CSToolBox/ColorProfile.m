classdef ColorProfile < ColorSet
    properties
        PRIM  % ColorData Array with 3 elements. 
        W     % ColorData
    end
    properties (Constant)        
    end
    methods(Static)
                
        function il_prf = get_illuminant_profile(illum)
        %
        % Given the name of a standard CIE illuminant, or any illuminant's 
        % spectrum, it computes the XYZ tristimulus of its reference white 
        % and 3 primaries.
        %
        % INPUT:
        %   illum --> The variable that identifies the standard CIE illuminant
        %             These variables are defined in the ColorSet subclass 
        %             and can be:
        %                 ILLUM_CIE_A, ILLUM_CIE_B, ILLUM_CIE_C, 
        %                 ILLUM_CIE_D50, ILLUM_CIE_D55, ILLUM_CIE_D65,
        %                 ILLUM_CIE_D75, ILLUM_CIE_E, or ILLUM_CIE_F7.
        %             It can also be a 2xN or Nx2 matrix containing the spd
        %             of the illuminant, where N is the number of
        %             wavelengths, and each pair contains the wavelength
        %             value and the energy corresponding to that
        %             wavelength.
        %
        % OUTPUT: 
        %     il_prf --> A 4x3 matrix where the first row contains the
        %                 XYZ values of the reference white, and the
        %                 remaining 3 rows contain each of the three
        %                 primaries. The structure is as follows:
        %
        %                | Xw Yw Zw | 
        %                | Xr Yr Zr | 
        %                | Xg Yg Zg | 
        %                | Xb Yb Zb | 
        %
        % EXAMPLE:
        %   il_prf = get_illuminant_profile(ColorUtils.ILLUM_CIE_D50);
        %
            
            Ws_xy = ColorUtils.CIE_E_xy2;
            
            % This variable indicates whether or not the input illuminant
            % represents a standard illuminant (via an illuminant
            % identification constant) or another illuminant's spd.
            std_illum = 1;
            
            if (numel(illum) == 1)
                switch (illum)
                    case ColorUtils.ILLUM_CIE_A                    
                        Wd_xy = ColorUtils.CIE_A_xy2;
                    case ColorUtils.ILLUM_CIE_B
                        Wd_xy = ColorUtils.CIE_B_xy2;
                    case ColorUtils.ILLUM_CIE_C
                        Wd_xy = ColorUtils.CIE_C_xy2;
                    case ColorUtils.ILLUM_CIE_D50
                        Wd_xy = ColorUtils.CIE_D50_xy2;
                    case ColorUtils.ILLUM_CIE_D55
                        Wd_xy = ColorUtils.CIE_D55_xy2;
                    case ColorUtils.ILLUM_CIE_D65
                        Wd_xy = ColorUtils.CIE_D65_xy2;
                    case ColorUtils.ILLUM_CIE_D75
                        Wd_xy = ColorUtils.CIE_D75_xy2;
                    case ColorUtils.ILLUM_CIE_E
                        Wd_xy = ColorUtils.CIE_E_xy2;
                    case ColorUtils.ILLUM_CIE_F7
                        Wd_xy = ColorUtils.CIE_F7_xy2;
                    otherwise
                        error('Wrong illuminant identifier.');
                end
            else
                std_illum = 0;
                [sx, sy] = size(illum);
                if (sx==2 && sy>2) || (sy==2 && sx>2)
                    % Load the color matching functions	
                    CMF = ColorUtils.loadCMF();

                    if (sx==2)
                        illum = illum';
                        [sx, sy] = size(illum);
                    end

                    b_start = illum(1,1);
                    b_end = illum(sx,1);
                    lambda = illum(2,1) - illum(1,1);
                    W_XYZd = ColorUtils.spect2XYZ(illum, CMF, b_start, b_end, lambda);
                else
                    error('Wrong illuminant format. Size must be Nx2 or 2xN (N=number of wavelengths)');
                end
            end
                
                        
            % Build the profile with respect to the known coordinates under
            % illuminant E.
            W_XYZs = ColorUtils.xyY2XYZ([Ws_xy 1]);
            prim_R_XYZs = ColorUtils.RGB2XYZ (ColorUtils.CIE_RED_RGB_E);
            prim_G_XYZs = ColorUtils.RGB2XYZ (ColorUtils.CIE_GREEN_RGB_E);
            prim_B_XYZs = ColorUtils.RGB2XYZ (ColorUtils.CIE_BLUE_RGB_E);
            
            % Apply a color adaptation transform from E to the destination
            % illuminant (std_illum).
            if (std_illum == 1)
                W_XYZd = ColorUtils.xyY2XYZ([Wd_xy 1]);
            end
            
            prim_R_XYZd = ColorUtils.CAT(prim_R_XYZs, W_XYZs, W_XYZd, ColorUtils.CAT_BRADFORD);
            prim_G_XYZd = ColorUtils.CAT(prim_G_XYZs, W_XYZs, W_XYZd, ColorUtils.CAT_BRADFORD);
            prim_B_XYZd = ColorUtils.CAT(prim_B_XYZs, W_XYZs, W_XYZd, ColorUtils.CAT_BRADFORD);

            il_prf = [W_XYZd; prim_R_XYZd'; prim_G_XYZd'; prim_B_XYZd'];

        end
        
        function M = compute_rgb2xyz_matrix_illum (illum)
        %
        % Given the name of a standard CIE illuminant, or any illuminant's 
        % spectrum, it computes the conversion matrix to pass from RGB to 
        % XYZ (the inverse operation is also possible by just computing the 
        % inverse matrix).
        %
        % INPUT:
        %   std_illum --> The variable that identifies the standard
        %                 CIE illuminant. These variables are defined
        %                 in the ColorSet subclass and can be:
        %                 ILLUM_CIE_A, ILLUM_CIE_B, ILLUM_CIE_C, 
        %                 ILLUM_CIE_D50, ILLUM_CIE_D55, ILLUM_CIE_D65,
        %                 ILLUM_CIE_D75, ILLUM_CIE_E, or ILLUM_CIE_F7.
        %             It can also be a 2xN or Nx2 matrix containing the spd
        %             of the illuminant, where N is the number of wavelenghts.
        %
        % OUTPUT:
        %   M --> The conversion matrix for the given profile, such that 
        %        [XYZ] = M * [RGB]
        %        (NOTE: to pass from XYZ to RGB we can use inv(M), so that 
        %        [RGB] = inv(M) * [XYZ]). 
        %
        % EXAMPLE:
        %
        %   M = compute_rgb2xyz_matrix_illum (ILLUM_CIE_D50);
        %
        
            % Get the tristimulus values for the profile
            prf = ColorProfile.get_illuminant_profile(illum);
            
            % Convert the tristimulus values into chromaticity coordinates,
            % as required by function compute_rgb2xyz_matrix            
            prim_xyY = zeros(3,3);
            prim_xyY (1,:) = ColorUtils.XYZ2xyY(prf(2,:));
            prim_xyY (2,:) = ColorUtils.XYZ2xyY(prf(3,:));
            prim_xyY (3,:) = ColorUtils.XYZ2xyY(prf(4,:));
            
            % Compute the conversion matrix
            M = ColorProfile.compute_rgb2xyz_matrix (prim_xyY(:,1:2), prf(1,:));
        end
        
        function M = compute_rgb2xyz_matrix (prim_xy, W_XYZ)
        %
        % Computes the matrix M to convert a given set of RGB Tristimulus Values
        % into the XYZ color space. The matrix is computed from the xyY values
        % for each of the three primaries (R, G and B) in the given color space
        % and the xyY values of the Reference white in the same space. 
        %
        % Source: http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
        %
        % INPUT:
        %   prim_xy --> A 3x2 matrix containing the xy coordinates of the primaries 
        %               (each row corresponds to a different primary).
        %
        %                | R_x R_y | 
        %                | G_x G_y |
        %                | B_x B_y |
        %   
        %   W_XYZ --> A 1x3 row vector containing the XYZ tristimulus of the Reference
        %             White ([W_X W_Y W_Z]).
        %
        %
        % OUTPUT: 
        %   M --> The conversion matrix for the given primaries, such that 
        %        [XYZ] = M * [RGB]
        %        (NOTE: to pass from XYZ to RGB we can use inv(M), so that 
        %        [RGB] = inv(M) * [XYZ]). 
        %
        % EXAMPLE:
        %
        %   M = compute_rgb2xyz_matrix (xy, W_XYZ);
        %
            
            Xr = prim_xy(1,1)/prim_xy(1,2);
            Yr = 1;
            Zr = (1-prim_xy(1,1)-prim_xy(1,2))/prim_xy(1,2);
            
            Xg = prim_xy(2,1)/prim_xy(2,2);
            Yg = 1;
            Zg = (1-prim_xy(2,1)-prim_xy(2,2))/prim_xy(2,2);
            
            Xb = prim_xy(3,1)/prim_xy(3,2);
            Yb = 1;
            Zb = (1-prim_xy(3,1)-prim_xy(3,2))/prim_xy(3,2);
            
            S = ([[Xr Xg Xb]; [Yr Yg Yb]; [Zr Zg Zb]]) \ W_XYZ';
            
            Mr = S(1).*[Xr Yr Zr];
            Mg = S(2).*[Xg Yg Zg];
            Mb = S(3).*[Xb Yb Zb];
            
            M = [Mr' Mg' Mb'];
            
            % Normalize the matrix to the white luminance.
            M = M ./ W_XYZ(2);
        end
        
        function plotCIEDiagram(color)
            if (nargin == 0)
                color = 'white';
            end
            locus = ColorProfile.compute_cie_locus_points(380, 780, 5);
            fill(locus(:,1), locus(:,2), color);
        end
    end
    methods(Static)
        function locus = compute_cie_locus_points(wl_start, wl_end, step)            
            if nargin < 2, step=5; end
            wl_range = ((wl_end-wl_start)/step)+1;
                        
            % Load the color matching functions	
            CMF = ColorUtils.loadCMF();
                        
            [sx,sy] = size(CMF);
            
            if ~((sx == wl_range && sy == 3) || (sx >= wl_range && sy == 4)), error('Wrong number of rows and/or columns for the color matching functions.'); end
            
            if sy == 4
                if ( (CMF(1,1) > wl_start) || (CMF(sx,1) < wl_end) ), error('Wrong color matching functions range.' ); end
                if ( (CMF(1,1) < wl_start) || (CMF(sx,1) > wl_end) )
                    temp = zeros(wl_range,3);
                    for i=1:3                        
                        t = ColorUtils.resample_array ([CMF(:,1) CMF(:,i+1)], wl_start:step:wl_end);
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
                locus(i,1) = XYZ(1)./sum(XYZ);
                locus(i,2) = XYZ(2)./sum(XYZ);
            end
            
            dx=(locus(wl_range,1)-locus(1,1))/(wl_range-1);
            dy=(locus(wl_range,2)-locus(1,2))/(wl_range-1);

            x=[locus(wl_range,1), locus(1,1)];
            y=[locus(wl_range,2), locus(1,2)];

            purple_boundary = [x(1)-dx:-dx:x(2);y(1)-dy:-dy:y(2)]';

            % Join the initial and final points of the spectrum locus to draw the purple boundary
            for i=1:length(purple_boundary)
                locus(wl_range+i,:) = [purple_boundary(i,1),purple_boundary(i,2)];
            end                        
        end
    end
    methods
        function obj = ColorProfile(PRIM, W)
            if (~isa(W,'ColorData')), error('White point must be an object of type ColorData.'); end
            if (numel(PRIM)~=3 || ~isa(PRIM,'ColorData')), error('Primaries must be expressed as a 3-element array of type ColorData.'); end            
            
            obj.PRIM = PRIM;
            obj.W = W;
            obj.num_colors = 0;
            obj.Colors = ColorData.empty(2,0);
        end        
        function CM = computeConvMatrix(obj)
            prim_xy = [obj.PRIM(1).colorCoord(); obj.PRIM(2).colorCoord(); obj.PRIM(3).colorCoord()];
            W_XYZ = obj.W.getTVs(ColorUtils.COLOR_SPACE_XYZ);
            CM = inv(ColorProfile.compute_rgb2xyz_matrix(prim_xy,W_XYZ)); 
        end
        
        function CM = XYZ2RGBConvMatrix(obj, method)            
            CM = inv(RGB2XYZConvMatrix(varargin));
        end
        
        function CM = RGB2XYZConvMatrix(obj, method)
            if (nargin < 2), method = 0; end
            
            switch (method)
                case 0 % Direct method
                    CM = [obj.PRIM(1).getTVs(ColorUtils.COLOR_SPACE_XYZ); obj.PRIM(2).getTVs(ColorUtils.COLOR_SPACE_XYZ); obj.PRIM(3).getTVs(ColorUtils.COLOR_SPACE_XYZ)]';
                case 1 % Luminance adaptation method
                    prim_xy = [obj.PRIM(1).colorCoord(); obj.PRIM(2).colorCoord(); obj.PRIM(3).colorCoord()];
                    W_XYZ = obj.W.getTVs(ColorUtils.COLOR_SPACE_XYZ);
                    CM = ColorProfile.compute_rgb2xyz_matrix(prim_xy,W_XYZ); 
            end
        end
        function obj = setSampling(obj, newstart, newend, newdelta)
            obj = obj.setSampling@ColorSet(newstart, newend, newdelta);
            for i=1:3
                obj.PRIM(i) = obj.PRIM(i).setSampling(newstart, newend, newdelta);
            end            
            obj.W = obj.W.setSampling(newstart, newend, newdelta);
        end        
        function h = plotGamut(obj, lineColor, lineWidth)
            if (nargin < 2)
                lineColor = 'blue';
            end
            if (nargin < 3)
                lineWidth = 2;
            end
            xy = zeros(3,2);
            for i=1:length(obj.PRIM)
                xy(i,:) = obj.PRIM(i).colorCoord(ColorUtils.COLOR_SPACE_XYZ);
            end
            h = plot([xy(:,1); xy(1,1)],[xy(:,2); xy(1,2)], 'color', lineColor, 'LineWidth', lineWidth);
        end
        function h = plotPrimaries(obj, colorSpace)
            if (nargin < 2)
                colorSpace = ColorUtils.COLOR_SPACE_XYZ;
            end
            h = zeros(length(obj.PRIM), 1);
            
            switch (colorSpace)
                case ColorUtils.COLOR_SPACE_XYZ
                    for i=1:length(obj.PRIM)                
                        %h(i) = obj.PRIM(i).plot(colorSpace,'format','CircleWidth',80,'fill',1);                        
                        switch(i)
                           case 1
                            MarkerColor='r';
                           case 2
                            MarkerColor='g';
                           case 3
                            MarkerColor='b';
                        end
                        h(i) = obj.PRIM(i).plot(colorSpace,'format',80,'MarkerFaceColor',MarkerColor,'MarkerEdgeColor',MarkerColor);
                    end
                case ColorUtils.COLOR_SPACE_RGB
                    M = obj.computeConvMatrix();                    
                    for i=1:length(obj.PRIM)
                        switch(i)
                           case 1
                            MarkerColor='r';
                           case 2
                            MarkerColor='g';
                           case 3
                            MarkerColor='b';
                        end
                        h(i) = obj.PRIM(i).plot(colorSpace,'CCM',M,'format',80,'MarkerFaceColor',MarkerColor,'MarkerEdgeColor',MarkerColor);
                    end                    
            end       
        end
        function h = plotAllColors(obj, size, colorSpace)
            if (nargin < 3)
                colorSpace = ColorUtils.COLOR_SPACE_XYZ;
            end
            
            if (nargin < 2)
                size=50;
            end
            
            h = zeros(obj.num_colors, 1);
            
            switch (colorSpace)
                case ColorUtils.COLOR_SPACE_XYZ
                    for i=1:obj.num_colors
                        h(i) = obj.Colors(i).plot(colorSpace,'format',size);
                    end
                case ColorUtils.COLOR_SPACE_RGB
                    M = obj.computeConvMatrix();
                    for i=1:obj.num_colors
                        h(i) = obj.Colors(i).plot(colorSpace,'CCM',M,'format',size);
                    end
            end                
        end        
        function h = plotWhite(obj, label, MarkerColor, size, colorSpace)
            add_label = true;
            
            if (nargin < 5)
                colorSpace = ColorUtils.COLOR_SPACE_XYZ;
            end
            
            if (nargin < 4)
                size=50;
            end
            
            if (nargin < 3)
                MarkerColor = 'white';
            end
            
            if (nargin < 2)
                label = '';
                add_label=false;
            end
            
            M = obj.computeConvMatrix();
            
            if (add_label)
                h = obj.W.plot(colorSpace,'CCM',M,'label',label,'format',size,'MarkerFaceColor',MarkerColor,'MarkerEdgeColor',MarkerColor);
            else
                h = obj.W.plot(colorSpace,'CCM',M,'format',size,'MarkerFaceColor',MarkerColor,'MarkerEdgeColor',MarkerColor);
            end
        end
    end
end













