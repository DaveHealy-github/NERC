classdef TriContMechFluidModel < TriContinuumReservoirModel

% DESCRIPTION: 
%   Base class model to set up fully coupled triple-permeability mechanical-fluid simulations. 

% PARAMETERS:
%   G            - Grid structure
%   rock         - Rock cell for fracture / macropore / micropore phase
%   fluid        - Fluid cell for fracture / macropore / micropore phase
%   mech_problem - Structure that contains the mechanical parameters of the system

% RETURNS:
%   class instance

    properties
        % Mechanical model
        mechModel;
        % List of primary variable names for the mechanical part
        MechPrimaryVars;
        % List of all the variable names for the mechanical part
        mechfds;
        % Fluid model
        fluidModel;
        % List of primary variable names for the fluid part
        FluidPrimaryVars;
        % List of all the variable names for the fluid part
        fluidfds;
        % Grouping for flow properties
        TC_PoroelasticPropertyFunctions; 
    end

    methods
        function model = TriContMechFluidModel(G, rock, fluid, mech_problem, varargin)
            model = model@TriContinuumReservoirModel(G, rock, fluid, varargin{:});
            % Process the grid for mechanical computation
            if ~ismember('createAugmentedGrid', model.G.type)
                model.G = createAugmentedGrid(model.G);
            end

            % Different fluid models may be used. 
            model.fluidModel = setupFluidModel(model); 
            model.fluidModel.extraStateOutput = 'True';
            model.fluidfds = model.fluidModel.getAllVarsNames(); % get variable names for fluid problem

            % Setup mech model
            rock_fracture = rock{1};
            rock_matrix = rock{2};
            rock_matrix2 = rock{3};
            model.mechModel = TriContMechMechanicModel(model.G, rock_fracture, rock_matrix, rock_matrix2, mech_problem); % Mechanical problem for FC sim
            model.mechfds = model.mechModel.getAllVarsNames(); % get variable names for mechanical problem
            model.TC_PoroelasticPropertyFunctions = [];
        end

        function fluidModel = setupFluidModel(model)
            error('Base class function not meant for direct use.');
        end

        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, varargin)

            error('Base class function not meant for direct use.');
        end

        function [fn, index] = getVariableField(model, name, varargin)
            if ismember(name, model.fluidfds)
                [fn, index] = model.fluidModel.getVariableField(name);
            elseif ismember(name, model.mechfds)
                [fn, index] = model.mechModel.getVariableField(name);
            else
                [fn, index] = getVariableField@TriContinuumReservoirModel(model, name, varargin{:});
            end
        end

        function mechTerm = computeStrainTerms(model, xd0, xd)
            error('Base class function not meant for direct use.');
        end

        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            fluidModel = model.fluidModel; 
            mechModel  = model.mechModel;
            [state, fluidReport] = fluidModel.updateState(state, problem, dx, []);
            [state, mechReport]  = mechModel.updateState(state, problem, dx, []);
            report = [];
        end

        function model = validateModel(model, varargin)
            % let the fluid model deal with the FacilityModel
            model.fluidModel.FacilityModel = model.FacilityModel;
            model.fluidModel = model.fluidModel.validateModel(); % setup wells 
            if isempty(model.TC_PoroelasticPropertyFunctions)
                model.TC_PoroelasticPropertyFunctions = TC_PoroelasticPropertyFunctions(model); 
            end
            return
        end

        function state = validateState(model, state)
           state = model.fluidModel.validateState(state);
           state = model.mechModel.validateState(state);
        end

        function containers = getStateFunctionGroupings(model)
            % called in getEquations when calling
            % model.initStateFunctionContainers() method
            containers = getStateFunctionGroupings@PhysicalModel(model);
            extra = {model.TC_PoroelasticPropertyFunctions};
            if ~isempty(model.FacilityModel)
                fm_props = model.FacilityModel.getStateFunctionGroupings();
                extra = [extra, fm_props];
            end
            extra = extra(~cellfun(@isempty, extra));
            containers = [containers, extra];
        end
    end
end
