classdef TriContMechWaterModel < TriContMechFluidModel

% DESCRIPTION: 
%   Model for fully coupled mechanical single-phase (water) fluid simulation.

% PARAMETERS:
%   G     - Simulation grid
%   rock  - Rock cell for fracture / macropore / micropore phases
%   fluid - Fluid cell for fracture / macropore / micropore phases

% RETURNS:
%   class instance

    properties

    end

    methods
        function model = TriContMechWaterModel(G, rock, fluid, mech_problem, varargin)
            model = model@TriContMechFluidModel(G, rock, fluid, mech_problem, ...
                                                 varargin{:});              
        end

        function fluidModel = setupFluidModel(model)
            fluidModel =TriContWaterFluidModel(model.G, {model.rock, model.rock_matrix, model.rock_matrix2},...
                                                {model.fluid, model.fluid_matrix, model.fluid_matrix2});
        end

        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, varargin)
            % Simulation options
            opt = struct('Verbose', mrstVerbose, ...
                         'reverseMode', false,...
                         'resOnly', false,...
                         'iteration', -1);  
            opt = merge_options(opt, varargin{:});

            % Shortcuts to sub-models
            fluidModel = model.fluidModel; 
            mechModel  = model.mechModel;  
                        
            % Extract variables and initialise as AD  properties at current timestep
            if ~opt.reverseMode
                [p, pm, pm2, wellSol, xd] = model.getProps(state, 'pressure', 'pressure_matrix',... 
                                                     'pressure_matrix2','wellSol', 'xd');
                [wellVars, wellVarNames, ~] = fluidModel.FacilityModel.getAllPrimaryVariables(wellSol);                                 
                [p, pm, pm2, wellVars{:}, xd] = initVariablesADI(p, pm, pm2, wellVars{:}, xd);
            else
                error('Reverse mode AD currently not implemented for TC-mech module')
            end
            
            % well variables
            [frac_index, mat_index, mat2_index] = fluidModel.findDCWells(wellSol);   
                        
            % Assemble linearised problem
            [w_eqs, w_eqsnames, w_eqstypes, state] = equationsTCWaterMech(state0, ...
                                                                          state, model, dt, ... 
                                                                          drivingForces, ...
                                                                         'iteration', ...
                                                                          opt.iteration);
                                                          
            [mech_eqs, mech_eqsnames, mech_eqstypes] = equationsTCPoroMechanics(xd, ...
                                                              mechModel, p, pm, pm2);                                                         
            eqs = horzcat(w_eqs, mech_eqs);
            names = {w_eqsnames{:}, mech_eqsnames{:}};
            types = {w_eqstypes{:}, mech_eqstypes{:}};
            primaryVars = {'pressure', 'pressure_matrix', 'pressure_matrix2', wellVarNames{:}, 'xd'};
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
        end       
        
        function mechTerm = computeStrainTerms(model, xd)
            % Method to compute strain terms (alpha_f:strain, alpha_m:strain, and alpha_m2:strain) used in the fluid content change equations
            G = model.G;
            s = model.mechModel.operators;
            cM = model.mechModel.constitutive_coefficients_object;
            
                alpha_f = cM.alpha_f(:,1);
                alpha_m = cM.alpha_m(:,1); 
                alpha_m2= cM.alpha_m2(:,1); 
            
            mechTerm.fracture = (s.div.*alpha_f)*xd./(G.cells.volumes);
            mechTerm.matrix = (s.div.*alpha_m)*xd./(G.cells.volumes);
            mechTerm.matrix2= (s.div.*alpha_m2)*xd./(G.cells.volumes);     
            
        end
       
    end
end
