function Ixyz = cam_rgb2xyz (Irgb_cam, Mcam)  
    
    [sx, sy, sz] = size(Irgb_cam);
    
    % Extract each color channel.    
    R = double(Irgb_cam(:,:,1));
    G = double(Irgb_cam(:,:,2));
    B = double(Irgb_cam(:,:,3));
    
    % Arrange the pixels in columns. Each column is one pixels and each row
    % corresponds to a different color channel.     
    Irgb_cam2 = [R(:), G(:), B(:)].';
    
    % Convert from Source RGB to Source XYZ, using the source
    % characterization matrices.
    Ixyz = (Mcam(:,1:3) * Irgb_cam2) + repmat(Mcam(:,4),1,size(Irgb_cam2,2));
    
    % Reshape the matrix into the original image dimensions.    
    %Irgb_d = reshape(Irgb_d',sx,sy,sz);
    Ixyz = reshape(Ixyz',sx,sy,sz);
    
end