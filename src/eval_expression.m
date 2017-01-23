function Y = eval_expression (f, X, model_info, params)
%EVAL_EXPRESSION - Evaluates a function call by passing it the corresponding
%                  arguments, extracted from a display characterization model structure.
%                  This structure is specific for the display characterization model defined
%                  for the MCTB toolbox.
%
% Y = EVAL_EXPRESSION(F, X, MODEL_INFO, PARAMS)
%
%Author: Victor Medina Heierle
%Last update: 23-Jan-2017
%
%    
    func_str = 'feval(f';
    
    for i=1:numel(model_info.coefficients)   
        if isfield(params, model_info.coefficients{i})
            func_str = sprintf('%s, %f',func_str,params.(model_info.coefficients{i}));
        end
    end
    
    func_str = sprintf('%s, X)',func_str);
    
    Y = eval(func_str);

end