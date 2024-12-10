classdef TriIsotropicStiffCoefficientModels < TriConstitutiveCoefficients
 
% PARAMETERS:
%   G     - Grid structure
%   model - Takes the mechanical model

% RETURNS:
%   class instance
    
    properties

    end
    
    methods
        function coefficients = TriIsotropicStiffCoefficientModels(model)
            coefficients = coefficients@TriConstitutiveCoefficients(model);                                 
        end
        
        function [C_m, C_m2, C_f] = caculateIntrinsicStiffnesses(coefficients, model)
             
            C_m  = Enu2C(model.mech.E_m, model.mech.nu_m, model.G);
            C_m2 = Enu2C(model.mech.E_m2, model.mech.nu_m2, model.G);
            C_f  = Enu2C(model.mech.E_f, model.mech.nu_f, model.G);
            
            d = model.G.griddim;
            if d == 3 
                identity = bsxfun(@times, ones(model.G.cells.num, 1), [1, 1, 1, 0, 0, 0]);
            else 
                identity = bsxfun(@times, ones(model.G.cells.num, 1), [1, 1, 0]);
            end
            
            K_m  = (1/(d^2))*doubledot(identity, doubledot(identity, C_m, model.G), model.G);
            C_m  = K_m;
            K_m2 = (1/(d^2))*doubledot(identity, doubledot(identity, C_m2, model.G), model.G);
            C_m2 = K_m2;
            K_f  = (1/(d^2))*doubledot(identity, doubledot(identity, C_f, model.G), model.G);
            C_f  = K_f;
        end 
        
        function [alpha_m, alpha_m2, alpha_f, invM_m, invM_m2, invM_f, invM12, invM13, invM23] = caculateCoefficients(coefficients, model)

            [alpha_m, alpha_m2, alpha_f, invM_m, invM_m2, invM_f, invM12, invM13, invM23] = ...
                TriCalcIsotropicCoefficientModels(coefficients, model);
        end
    end
end

