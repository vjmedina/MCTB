function I2 = img_array_mult (I, M)
    
    [sx, sy, sz] = size(I);
    
    % Extract each color channel.    
    R = double(I(:,:,1));
    G = double(I(:,:,2));
    B = double(I(:,:,3));
    
    % Arrange the pixels in columns. Each column is one pixels and each row
    % corresponds to a different color channel.     
    I2 = [R(:), G(:), B(:)].';
    
    % Multiply the new array by the matrix.
    I2 = M * I2;
    
    % Reshape the image into the original dimensions.
    I2 = reshape(I2',sx,sy,sz);
    
end