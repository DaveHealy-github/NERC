classdef TC_PoroelasticPropertyFunctions < StateFunctionGrouping
    % Grouping for describing rock-fluid properties  
    
    properties
        PoroelasticMatrixPoro 
        PoroelasticMatrix2Poro 
        PoroelasticFracturePoro 
        Strain 
        EffectiveStress 
        TotalStress 
        PoroelasticFracturePerm  
    end

    methods
        function props = TC_PoroelasticPropertyFunctions(model)
            props@StateFunctionGrouping();
            props.PoroelasticMatrixPoro = TriPoroelasticMatrixPoro(model);
            props.PoroelasticMatrix2Poro = TriPoroelasticMatrix2Poro(model);
            props.PoroelasticFracturePoro = TriPoroelasticFracturePoro(model);
            props.Strain = TriStrain(model);
            props.EffectiveStress = TriEffectiveStress(model);
            props.TotalStress = TriTotalStress(model);
            props.PoroelasticFracturePerm = TriPoroelasticFracturePerm(model);
            
            % Define storage field in state
            props.structName = 'TC_PoroelasticProps';
        end
    end
end