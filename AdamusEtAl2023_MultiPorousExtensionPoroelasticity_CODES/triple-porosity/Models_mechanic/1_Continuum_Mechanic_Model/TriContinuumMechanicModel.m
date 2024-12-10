classdef TriContinuumMechanicModel < MechanicModel

% DESCRIPTION:
% Model that can handle a triple-porosity mechanical system. The inputs can be standard
% boundary conditions and body forces for mechanics, but also a fluid
% pressure, which enters the mechanical equation as grad(pressure).

% PARAMETERS:
%   G            - Grid structure
%   rock         - Rock structure for the fracture phase
%   rock_matrix  - Rock structure for the macropore phase
%   rock_matrix2 - Rock structure for the micropore phase
%   mech_problem - Struct containing parameters required to instatiate
%                  mechanical problem

% RETURNS:
%   class instance

    properties
        rock_matrix; 
        rock_matrix2; 
        constitutive_coefficients_object;
    end

    methods
        function model = TriContinuumMechanicModel(G, rock, rock_matrix, rock_matrix2,  mech_problem, varargin)
        % Constructor for TriContinuumMechanicModel
            model = model@MechanicModel(G, rock, mech_problem, varargin{:});
            
            % Physical properties of matrix rock
            model.rock_matrix  = rock_matrix;
            model.rock_matrix2  = rock_matrix2;

             % Create required coefficients
            coefficient_model_type = evalCoefficientModelType(model);          
            switch coefficient_model_type
                case 'requireInformation'
                    disp('Missing mechanical properties')
                case 'TriIsotropicStiffCoefficientModels'
                    coefficient_model_handle = str2func(coefficient_model_type);
                    model.constitutive_coefficients_object = coefficient_model_handle(model);
            end   
        end

        
        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, ...
                                                        varargin)
            % Assemble the equations for the mechanical system
            opt = struct('Verbose'       , mrstVerbose , ...
                         'reverseMode'   , false       , ...
                         'scaling'       , []          , ...
                         'resOnly'       , false       , ...
                         'history'       , []          , ...
                         'iteration'     , -1          , ...
                         'stepOptions'   , []          , ...
                         'addflux'       , false); % Compatibility only

            opt = merge_options(opt, varargin{:});
                        
            % The fluid pressure stimulates the mechanical system. It is
            % given as a driving force.
            fluidp = drivingForces.fluidp;
            fluidp_matrix = drivingForces.fluidp_matrix;
            fluidp_matrix2 = drivingForces.fluidp_matrix2;
            
            xd = model.getProps(state, 'xd'); 

            if ~opt.resOnly
                xd = initVariablesADI(xd);  
            end

            eqs = equationsTCPoroMechanics(xd, model, fluidp, fluidp_matrix, fluidp_matrix2); 

            primaryVars = {'xd'};
            names = {'disp'};
            types = {'disp_dofs'};

            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt); 
            problem.iterationNo = opt.iteration;
        end

        function forces = getValidDrivingForces(model)

            forces.fluidp = [];
            forces.fluidp_matrix = [];
            forces.fluidp_matrix2 = [];
        end

        function [fn, index] = getVariableField(model, name, varargin)
        % Get the index/name mapping for the model (such as where
        % pressure or water saturation is located in state)
            switch(lower(name))
              case {'uu'}
                % Displacement field as a matrix (one column per Cartesian direction)
                fn = 'uu';
                index = ':';
              case {'u'}
                % Displacement field given as a column vector where the cartesian components are
                % stabbed.
                fn = 'u';
                index = ':';
              case {'xd'}
                % Same as 'u' but the degree of freedom where the Dirichlet conditions (fixed
                % displacement) are removed
                fn = 'xd';
                index = 1;
              case {'stress'}
                % Stress tensor (Voigt notation, one column per component)
                fn = 'stress';
                index = ':';
              case {'strain'}
                % Strain tensor (Voigt notation, one column per component)
                fn = 'strain';
                index = ':';
              case {'vdiv'}
                % Volume weighted divergence field of displacement.
                fn = 'vdiv';
                index = ':';
              otherwise
                % This will throw an error 
                [fn, index] = getVariableField@MechanicModel(model, name);
            end
        end

        function [primaryVars, fds] = getAllVarsNames(model)
        % List the variables that are used by the model. Used when 
        % coupling the mechanical model with a fluid model.
            primaryVars = {'xd'};
            fds = {'xd', 'uu', 'u', 'local_strain', 'stress', 'strain', ...
                   'strain_frac', 'strain_mat', 'strain_mat2', 'vdiv'};
        end
        
        function [coefficient_model_type] = evalCoefficientModelType(model)
        % Function to evaluate the type of coefficient model that we need
        % to use. Returns a string, which will be converted to a function
        % handle used to identify the appropriate models to be used in the
        % calculation of the coefficient models.
            if ~isfield(model.mech, 'E_m') && ~isfield(model.mech, 'nu_m')
                coefficient_model_type = 'requireInformation';     
            else                  
                coefficient_model_type = 'TriIsotropicStiffCoefficientModels'; 
            end
        end
        
    end
end
