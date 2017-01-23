function [info] = get_calibration_model_info (model)
%GET_CALIBRATION_MODEL_INFO - Given a display characterization model id number,
%                             as defined in the MCTB toolbox, this function creates 
%                             and returns a data structures containing all the neccesary
%                             data relative to the given model.
%
% INFO = GET_CALIBRATION_MODEL_INFO(MODEL)
%
%Author: Victor Medina Heierle
%Last update: 23-Jan-2017
%
%   
    % Define the equations for each model
    % Forward models
    GOG = '(a*x+b)^m';
    GOGO = '(a*x+b)^m + c';
    GGO = '(a*x)^m +b';
    simple = 'x^m';
    
    %Backward models
    GOG2 = 'a*((x^m)-b)';
    GOGO2 = 'a*((max(0,(x-c))^m)-b)';
    GGO2 = '(a*(max(0,(x-b))^m))';
    simple2 = 'x^m';
    
    % Define the fitting options for each model
    switch(model)
        case 1 % backward GOG model
            eq = GOG2;
            lower = [0 0 0.01];
            upper = [Inf 1 1];
            coeff = {'a','b','m'};
            ind_var = 'x';
        case 2 % backward GOGO model
            eq = GOGO2;
            lower = [0 0 0 0];
            upper = [Inf 1 1 1];
            coeff = {'a','b','c','m'};
            ind_var = 'x';
        case 3 % backward GGO model
            eq = GGO2;
            lower = [0 0 0];
            upper = [Inf 1 1];
            coeff = {'a','b','m'};
            ind_var = 'x';
        case 4 % backward simple model
            eq = simple2;
            lower = 0.01;
            upper = 1;
            coeff = {'m'};
            ind_var = 'x';
        case 5 % forward GOG model
            eq = GOG;
            lower = [0 0 1];
            upper = [2 2 5];
            coeff = {'a','b','m'};
            ind_var = 'x';
        case 6 % forward GOGO model
            eq = GOGO;
            lower = [0 0 -10 1];
            upper = [2 5 10 100];
            coeff = {'a','b','c','m'};
            ind_var = 'x';
        case 7 % forward GGO model
            eq = GGO;
            lower = [0 0 1];
            upper = [2 100 5];
            coeff = {'a','b','m'};
            ind_var = 'x';
        case 8 % forward Simple model
            eq = simple;
            lower = 1;
            upper = 100;
            coeff = {'m'};
            ind_var = 'x';
        otherwise % backward GOG model
            eq = GOG2;
            lower = [0 0 0.01];
            upper = [Inf 5 1];
            coeff = {'a','b','m'};
            ind_var = 'x';
    end
    
    % Declare the parameter structure where we will store the results. 
    %info = struct('expression',eq,'LowerBounds',lower,'UpperBounds',upper,'coefficients',coeff,'independent',ind_var);
    info = struct('expression',eq,'LowerBounds',lower,'UpperBounds',upper,'coefficients','','independent',ind_var);
    info.coefficients=coeff;
    
end
    
    
    