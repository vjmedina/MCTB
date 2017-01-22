function Iout = display_color_transform (Iin, display_model_data)

[sx, sy, sz] = size(Iin);
Xr = Iin(:,:,1);
Xg = Iin(:,:,2);
Xb = Iin(:,:,3);

model_info = get_calibration_model_info (display_model_data.model);
f = fittype(model_info.expression);

Yr = eval_expression (f, Xr(:), model_info, display_model_data.red_params);
Yg = eval_expression (f, Xg(:), model_info, display_model_data.green_params);
Yb = eval_expression (f, Xb(:), model_info, display_model_data.blue_params);

Yr = reshape(Yr, sx, sy);
Yg = reshape(Yg, sx, sy);
Yb = reshape(Yb, sx, sy);

Iout=Iin;
Iout(:,:,1) = Yr;
Iout(:,:,2) = Yg;
Iout(:,:,3) = Yb;

%     switch(model_data.model)
%         case 1
%             Iout = (params.a .* Iin + params.b).^params.gamma;
%         case 2
%             Iout = (params.a .* Iin + params.b).^params.gamma + params.c;
%         case 3
%             Iout = (params.a .* Iin).^params.gamma + params.b;
%         case 4
%             Iout = Iin.^params.gamma;
%         otherwise
%             Iout = Iin;
%     end            
end