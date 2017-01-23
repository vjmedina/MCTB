function [params, R2] = estimate_display_model (X, Y, model)
%ESTIMATE_DISPLAY_MODEL - Given a series of corresponding input(Y) and output(X) data, 
%                         and a specific characterization model, this function estimates
%                         the neccesary parameters for the model, as well as the corresponding
%                         goodness of fit parameter, R2.
%
% [PARAMS, R2] = ESTIMATE_DISPLAY_MODEL(X, Y, MODEL)
%
%Author: Victor Medina Heierle
%Last update: 23-Jan-2017
%
%   
    warning('off', 'curvefit:fit:noStartPoint');
    
    if (nargin < 1)
        model = 3;
    end
    
    % Clip negative values to zero.
    X = max(X,zeros(size(X)));
    Y = max(Y,zeros(size(Y)));
    
    % Read the model information
    model_info = get_calibration_model_info (model);
        
    % Declare the parameter structure where we will store the results. 
    params = struct('a',1,'b',0,'c',0,'m',1);
    
    % Define the fitting options for each model
    s = fitoptions('Method','NonlinearLeastSquares', 'Lower', model_info.LowerBounds, 'Upper', model_info.UpperBounds);
    f = fittype(model_info.expression,'coefficients',model_info.coefficients,'independent',model_info.independent,'options', s);

    % Perform the Non-Linear Least-Squares data fitting according to the
    % chosen model.
    
    R2 = 0;
    num_loops = 0;
    while (R2 < 0.9 && num_loops < 5)
        if (size(X,1)==1)
            [f, gof] = fit(X',Y',f,s);
        else
            [f, gof] = fit(X,Y,f,s);
        end    
        % Store the results.
        R2 = gof.rsquare;
        num_loops = num_loops + 1;
    end
           
    coeff_names = coeffnames(f);
    coeff_values = coeffvalues(f);
    for i=1:numel(coeff_names)
        if (isfield(params, coeff_names{i}))
            params.(coeff_names{i}) = coeff_values(i);
        end
    end
        
%     if (model ~= 4 && model ~= 8) % Simple models
%         params.a = f.a;
%         params.b = f.b;
%     end
%     params.m = f.m;
%     if (model == 2 || model == 6) % GOGO model
%         params.c = f.c;
%     end
%     
end