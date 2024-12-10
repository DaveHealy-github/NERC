classdef TriTotalStress < StateFunction
    % Total stress imposed on model
    % Final equation is equivalent to equation 24 in Adamus et al. (2023))

    properties
    end
    
    methods
        function gp = TriTotalStress(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'pressure', 'pressure_matrix','pressure_matrix2'}, 'state');
            gp = gp.dependsOn({'EffectiveStress'});
        end
        function v = evaluateOnDomain(prop, model, state)
            [p, pm, pm2] = model.getProps(state, 'pressure', 'pressure_matrix','pressure_matrix2');
            effective_stress = prop.getEvaluatedDependencies(state, 'EffectiveStress');
            cM = model.mechModel.constitutive_coefficients_object;
            if model.G.griddim == 2
                pI = bsxfun(@times, p, [1, 1, 0]);
                pmI = bsxfun(@times, pm, [1, 1, 0]);
                pm2I = bsxfun(@times, pm2, [1, 1, 0]);
            else
                pI = bsxfun(@times, p, [1, 1, 1, 0, 0, 0]);
                pmI = bsxfun(@times, pm, [1, 1, 1, 0, 0, 0]);
                pm2I = bsxfun(@times, pm2, [1, 1, 1, 0, 0, 0]);
            end
            v = effective_stress - cM.alpha_m2.*pm2I - cM.alpha_m.*pmI - cM.alpha_f.*pI;  
        end
    end
end