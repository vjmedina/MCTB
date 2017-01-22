function Irgb = disp_xyz2rgb (Ixyz_disp, Mdisp)  
    
    [sx, sy, sz] = size(Ixyz_disp);
    
    % Extract each color channel.    
    X = double(Ixyz_disp(:,:,1));
    Y = double(Ixyz_disp(:,:,2));
    Z = double(Ixyz_disp(:,:,3));
    
    % Arrange the pixels in columns. Each column is one pixels and each row
    % corresponds to a different color channel.     
    Ixyz_disp2 = [X(:), Y(:), Z(:)].';
    
    % Convert from Source RGB to Source XYZ, using the source
    % characterization matrices.
    Irgb = (Mdisp(:,1:3) * Ixyz_disp2);
    
    % Reshape the matrix into the original image dimensions.    
    %Irgb_d = reshape(Irgb_d',sx,sy,sz);
    Irgb = reshape(Irgb',sx,sy,sz);
    
end