classdef TriPoroelasticMatrix2Poro < StateFunction
    % Effective porosity of a micropore set 
    % Final equation describes the fluid content change 
    % that is equivalent to equation 25 in Adamus et al. (2023))

    properties
    end
    
    methods
        function gp = TriPoroelasticMatrix2Poro(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'pressure', 'pressure_matrix', 'pressure_matrix2', 'xd'}, 'state');
        end
        function poro_m2 = evaluateOnDomain(prop, model, state)
            % Get effective fracture porosity given changes in strain and 
            % macropore, micorpore, and fracture pressures. 
            cM = model.mechModel.constitutive_coefficients_object;
            ref_poro_m2 = model.rock_matrix2.poro;
            [p, pm, pm2, xd] = model.getProps(state, 'pressure', 'pressure_matrix', 'pressure_matrix2', 'xd');
            mechTerm = model.computeStrainTerms(xd);
            poro_m2 = ref_poro_m2 + (mechTerm.matrix2 + cM.invM12.*pm + cM.invM_m2.*pm2 + cM.invM23.*p); % mechTerm.matrix2 = alpha_m2:strain
        end
    end
end