classdef TriPoroelasticFracturePerm < StateFunction
    % Permeability arising from poroelastic effects
    properties
    end
    
    methods
        function gp = TriPoroelasticFracturePerm(model, varargin)
            gp@StateFunction(model, varargin{:});
            if isfield(model.rock, 'nonLinearPerm')
                gp = gp.dependsOn({'TriPoroelasticFracturePoro'});
            end
        end
        function perm = evaluateOnDomain(prop, model, state)
            r = model.rock;
            if isfield(r, 'nonLinearPerm')
                poro_f = prop.getEvaluatedDependencies(state, 'TriPoroelasticFracturePoro');
                perm = prop.evaluateFunctionOnDomainWithArguments(r.nonLinearPerm, poro_f);
            end
        end
    end
end