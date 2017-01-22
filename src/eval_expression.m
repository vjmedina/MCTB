function Y = eval_expression (f, X, model_info, params)
    
    func_str = 'feval(f';
    
    for i=1:numel(model_info.coefficients)   
        if isfield(params, model_info.coefficients{i})
            func_str = sprintf('%s, %f',func_str,params.(model_info.coefficients{i}));
        end
    end
    
    func_str = sprintf('%s, X)',func_str);
    
    Y = eval(func_str);

end