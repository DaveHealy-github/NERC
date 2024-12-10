classdef TriConstitutiveCoefficients

% DESCRIPTION: 
%   Base class for use in calculating poroelastic coefficients belonging to the 
%   triple-porosity model.  

% PARAMETERS:
%   G     - Grid structure
%   model - Takes the mechanical model

% RETURNS:
%   class instance
    
    properties

        alpha_m;
        alpha_m2;
        alpha_f;
        invM_m;
        invM_m2;
        invM_f;
        invM12;
        invM13;
        invM23;
       
    end
    
    methods
        function coefficients = TriConstitutiveCoefficients(model)
            % Check if we have the necessary fields required to calculate
            % the poroelastic coefficients
                                
            assert(isfield(model.rock, 'poro'),...
                    'rock struct requires a fracture porosity field - "poro"')
            assert(isfield(model.rock_matrix, 'poro'),...
                    'rock_matrix struct requires a matrix porosity field - "poro"')
            assert(isfield(model.rock_matrix2, 'poro'),...
                    'rock_matrix2 struct requires a matrix2 porosity field - "poro"')
            assert(isfield(model.rock, 'vol_fraction'),...
                    'rock struct requires a fracture volume fraction field - "vol_fraction"')
            assert(isfield(model.rock_matrix, 'vol_fraction'),...
                    'rock_matrix struct requires a matrix volume fraction field - "vol_fraction"')
            assert(isfield(model.rock_matrix2, 'vol_fraction'),...
                    'rock_matrix2 struct requires a matrix2 volume fraction field - "vol_fraction"')
            
              
            % calculate effective coefficients
            [coefficients.alpha_m, coefficients.alpha_m2, coefficients.alpha_f, ... 
                coefficients.invM_m, coefficients.invM_m2, coefficients.invM_f, coefficients.invM12, coefficients.invM13, coefficients.invM23]...
                = caculateCoefficients(coefficients, model);
        end
       
       
        function [alpha_m, alpha_m2, alpha_f, invM_m, invM_m2, invM_f, invM12, invM13, invM23] = caculateCoefficients(coefficients, model)
            error('Base class function not meant for direct use.');
        end
    end
end

