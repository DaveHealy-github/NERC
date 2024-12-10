classdef TriPoroelasticFracturePoro < StateFunction
    % Effective porosity of a fracture set 
    % Final equation describes the fluid content change 
    % that is equivalent to equation 25 in Adamus et al. (2023))

    properties
    end
    
    methods
        function gp = TriPoroelasticFracturePoro(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'pressure', 'pressure_matrix', 'pressure_matrix2', 'xd'}, 'state');
        end
        function poro_f = evaluateOnDomain(prop, model, state)
            % Get effective fracture porosity given changes in strain and 
            % macropore, micorpore, and fracture pressures. 
            cM = model.mechModel.constitutive_coefficients_object;
            ref_poro_f = model.rock.poro;
            [p, pm, pm2, xd] = model.getProps(state, 'pressure', 'pressure_matrix', 'pressure_matrix2', 'xd');
            mechTerm = model.computeStrainTerms(xd);
            poro_f = ref_poro_f + (mechTerm.fracture + cM.invM13.*pm + cM.invM23.*pm2 + cM.invM_f.*p); % mechTerm.fracture = alpha_f:strain
        end
    end
end