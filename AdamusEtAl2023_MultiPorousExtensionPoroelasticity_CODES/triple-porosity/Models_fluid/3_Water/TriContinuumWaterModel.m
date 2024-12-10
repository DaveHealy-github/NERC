classdef TriContinuumWaterModel < TriContinuumReservoirModel

% DESCRIPTION: 
%   Single phase water model that will be required by the fully coupled 
%   and fixed stress triple-porosity water models.

% PARAMETERS:
%   G     - Simulation grid.
%   rock  - Rock cell for fracture / macropore / micropore phases
%   fluid - Fluid cell for fracture / macropore / micropore phases

% RETURNS:
%   class instance

    properties

    end
    
    methods
        function model = TriContinuumWaterModel(G, rock, fluid, varargin)
            model = model@TriContinuumReservoirModel(G, rock, fluid);
            % Only enable water
            model.oil = false;
            model.gas = false;
            model.water = true;
            model.useCNVConvergence = false;   
            model = merge_options(model, varargin{:});
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            % we can just keep the equations for the single-phase instance, as they will not be used anyway. 
            [problem, state] = equationsWater(state0, state, model,...
                                              dt, drivingForces,...
                                              varargin{:});           
        end
        
        function [frac_index, mat_index, mat2_index] = findDCWells(model, wellSol)
        % By definition the triple-porosity model assumes flow from each
        % porosity cluster. This function helps to find the wells in
        % wellSol that correspond to macropore, micropore, and fracture continua. 
        % It enforces that when we add in a well, we add a well for each
        % continua. Rows of the wellSol struct that correspond to fracture
        % and macropore and micropore wells are denoted frac_index, mat_index, and mat2_index 
            if ~ isempty(wellSol) 
                frac_index = find(contains({wellSol.name}, 'frac'));
                mat_index = find(contains({wellSol.name}, 'mat'));
                 mat2_index = find(contains({wellSol.name}, 'mat2'));
                assert(~ isempty(frac_index) && ~ isempty(mat_index) && ~ isempty(mat2_index), ...
                    ['well names are required to include fracture, macropore, and micropore labels e.g. '...
                     'Prod_frac, Prod_mat, Prod_mat2. If a single well connection is required '...
                     'for example to a fracture, create a dummy "matrix" and "matrix2" wells with a '...
                     'dummy structure containing a zero permeability field entry'])
            else
                % dummy indices
                frac_index = 1; mat_index = 2; mat2_index = 3;
            end
        end
        
    end
end

