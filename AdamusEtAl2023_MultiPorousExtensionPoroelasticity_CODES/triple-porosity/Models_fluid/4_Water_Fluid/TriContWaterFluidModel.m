classdef TriContWaterFluidModel < TriContinuumWaterModel

% DESCRIPTION: 
%   This model is for the fluid part of a fully coupled triple continuum 
%   poroelastic simulation. A few functionalities that are needed for the coupled solver are added.

% PARAMETERS:
%   G     - Simulation grid.
%   rock  - Rock cell for fracture / macropore / micropore phases
%   fluid - Fluid cell for fracture / macropore / micropore phases

% RETURNS:
%   class instance

    properties
        primaryVarNames;
    end

    methods
        function model = TriContWaterFluidModel(G, rock, fluid, varargin)
            model = model@TriContinuumWaterModel(G, rock, fluid);
            model.primaryVarNames = {'pressure', 'pressure_matrix', 'pressure_matrix2'}; % well variables
                                                                     
            model = merge_options(model, varargin{:});
        end

        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, ...
                                                        varargin)
            error(['The TriContWaterFluidModel is not meant for being called directly. ' ...
                   'It should be used coupled with a mechanical model, derived ' ...
                   'from TriContMechFluidModel.'])
        end

        function varnames = getAllPrimaryVariables(model, state)
        % List all the primary variables that are recognized and can be
        % handled by the model, used by the updateState member function.
            varnames = model.primaryVarNames;
            [~, wellVarNames, ~] = ...
                model.FacilityModel.getAllPrimaryVariables(state.wellSol);
            varnames = {varnames{:}, wellVarNames{:}};
        end

        function fds = getAllVarsNames(model)
        % List all the variables that are recognized and can be handled by the model
            fds = {'pressure', 'pressure_matrix', 'pressure_matrix2','wellSol'};
        end

        function [state, report] = updateState(model, state, problem, dx, ...
                                               drivingForces)
        % Updates the state variables that belong to the model, that is the fluid
        % variables. The mechanical variables in the state will be updated
        % by the mechanical model, see TriContMechanicMechModel.
            vars = problem.primaryVariables;
            ind = false(size(vars));
            fluidvars = model.getAllPrimaryVariables(state);

            [lia, loc] = ismember(fluidvars, vars);
            assert(all(lia), 'The primary variables are not set correctly.');
            ind(loc) = true;
            problem.primaryVariables = vars(ind);
            dx   = dx(ind);

            [state, report] = updateState@TriContinuumWaterModel(model, state, problem, ...
                                                                  dx, drivingForces);
        end

    end

end
