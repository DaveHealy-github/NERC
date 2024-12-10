classdef TriContinuumReservoirModel < TriPorosityReservoirModel

% DESCRIPTION:
%   Base flow class for all models that show triple-permeability behaviour.

% REQUIRED PARAMETERS:
%   G     - Simulation grid.
%   rock  - Rock cell for fracture / macropores / micropores
%   fluid - Fluid cell for fracture / macropores / micropores

% RETURNS:
%   Class instance.


    properties    
        operators_matrix  % required to consider flow in macropore phase
        operators_matrix2 % required to consider flow in micropore phase
    end

    methods
        % --------------------------------------------------------------------%
        function model = TriContinuumReservoirModel(G, rock, fluid, varargin)
            
            model = model@TriPorosityReservoirModel(G, rock, fluid);
            model.operators_matrix = setupOperatorsTPFA(G, model.rock_matrix, 'deck', model.inputdata);
            model.operators_matrix2= setupOperatorsTPFA(G, model.rock_matrix2, 'deck', model.inputdata);
        end
    end
end
