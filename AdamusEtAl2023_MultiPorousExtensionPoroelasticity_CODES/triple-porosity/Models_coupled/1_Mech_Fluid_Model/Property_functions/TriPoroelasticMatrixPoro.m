classdef TriPoroelasticMatrixPoro < StateFunction
    % Effective porosity of a macropore set 
    % Final equation describes the fluid content change 
    % that is equivalent to equation 25 in Adamus et al. (2023))

    properties
    end
    
    methods
        function gp = TriPoroelasticMatrixPoro(model, varargin)
            gp@StateFunction(model, varargin{:});
            gp = gp.dependsOn({'pressure', 'pressure_matrix', 'pressure_matrix2', 'xd'}, 'state');
        end
        function m_poro = evaluateOnDomain(prop, model, state)
            % Get effective fracture porosity given changes in strain and 
            % macropore, micorpore, and fracture pressures. 
            cM = model.mechModel.constitutive_coefficients_object;
            ref_poro_m = model.rock_matrix.poro;
            [p, pm, pm2, xd] = model.getProps(state, 'pressure', 'pressure_matrix', 'pressure_matrix2', 'xd');
            mechTerm = model.computeStrainTerms(xd);
            m_poro = ref_poro_m + (mechTerm.matrix + cM.invM_m.*pm + cM.invM12.*pm2 + cM.invM13.*p); % mechTerm.matrix = alpha_m:strain
        end
    end
end