classdef TriContMechMechanicModel < TriContinuumMechanicModel

% DESCRIPTION: 
% This model is for the mechanical part of a fully coupled model. 
% It adds some functionalities that are needed to couple the solver
% with a fluid model.

% PARAMETERS:
%   G            - Grid structure
%   rock         - Rock structure for fracture phase
%   rock_matrix  - Rock structure for macropore phase
%   rock_matrix2 - Rock structure for micropore phase
%   mech_problem - Structure that contains the mechanical parameters of the system

% RETURNS:
%   class instance

    properties
        % Primary variables that are used in the mechanical system
        primaryVarNames;
    end

    methods
        function model = TriContMechMechanicModel(G, rock, rock_matrix, rock_matrix2, mech_problem, varargin)
        % constructor
            model = model@TriContinuumMechanicModel(G, rock, rock_matrix, rock_matrix2, mech_problem, varargin{:});
            model.primaryVarNames = {'xd'};
            model = merge_options(model, varargin{:});
        end

        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, varargin)
            error(['The TriContMechMechanicalModel is not meant for being called directly. ' ...
                   'It should be used coupled with a mechanical model, derived ' ...
                   'from TriContMechFluidModel.'])
        end

        function  fds = getAllVarsNames(model) % get function for the varsnames
        % list all the variables that are recognized and can be handled by the model
            fds = {'xd', 'uu', 'u', 'stress', 'strain', 'vdiv'};
        end

        function varnames = getAllPrimaryVariables(model) % get function
        % list all the primary variables that are recognized and can be
        % handled by the model, used by the updateState member function.
            varnames = model.primaryVarNames;
        end

        function [state, report] = updateState(model, state, problem, dx, ...
                                               drivingForces)
            % updates the state variable that belongs to the model, that is
            % the mechanical variables. The fluid variables in the state will
            % be updated by the fluid model. 
            vars = problem.primaryVariables;
            ind = false(size(vars));
            mechvars = model.getAllPrimaryVariables();
            [lia, loc] = ismember(mechvars, vars);
            assert(all(lia), 'The primary variables are not set correctly');

            ind(loc) = true;
            problem.primaryVariables = vars(ind);
            dx   = dx(ind);

            [state, report] = updateState@TriContinuumMechanicModel(model, state, problem, ...
                                                        dx, drivingForces);                        
            state = addDerivedQuantities(model, state);
        end

    end
end
