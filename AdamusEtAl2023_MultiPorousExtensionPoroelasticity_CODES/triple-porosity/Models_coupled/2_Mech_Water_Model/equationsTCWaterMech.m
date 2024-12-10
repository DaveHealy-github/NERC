function [eqs, names, types, state] = equationsTCWaterMech(state0, state, model, dt, drivingForces, varargin)

% DESCRIPTION: 
%   Assembles the residuals for the triple-porosity fluid continuity
%   equations, akin to equations (33)-(34) in Adamus et al. (2023)

% PARAMETERS:
%   state0                   - state    (previous time step)
%   state                    - state (current time step)
%   model                    - model class instance that is used.
%                              fluidModel in this case.
%   dt                       - time step
%   drivingForces            - structure that gathers the well parameters and boundary conditions.

% RETURNS:
%   eqs   - The residual values as ADI variables (that is with the Jacobian)
%           if the inputs were also ADI.
%   names - The name of each equations
%   types - The type of each equations
%   state - Some field related to well control of the state variables may be updated.

% Note that state is given only for output
opt = struct('Verbose', mrstVerbose, ...
         'reverseMode', false,...
         'resOnly', false,...
         'iteration', -1);  % Compatibility only
opt = merge_options(opt, varargin{:});

% Inovking commonly used parts of the model and forces
f = model.fluidModel.fluid;
s = model.fluidModel.operators;
f_m = model.fluidModel.fluid_matrix;
s_m = model.fluidModel.operators_matrix;   
f_m2 = model.fluidModel.fluid_matrix2;
s_m2 = model.fluidModel.operators_matrix2; 
T_f = s.T;


%% -------------------------------------------------------------------------
% Initialization of independent variables
% properties at current timestep, and current iteration level
[p, pm, pm2, wellSol, xd] = model.getProps(state, 'pressure', 'pressure_matrix', 'pressure_matrix2', 'wellSol', 'xd');
% properties at previous timestep
[p0, pm0, pm20, wellSol0, xd0] = model.getProps(state0, 'pressure', 'pressure_matrix', 'pressure_matrix2', 'wellSol', 'xd');

 [wellVars, ~, wellMap] = ...
     model.fluidModel.FacilityModel.getAllPrimaryVariables(wellSol);
if ~opt.reverseMode
    [p, pm, pm2, wellVars{:}, xd] = initVariablesADI(p, pm, pm2, wellVars{:}, xd);
else
    error('Reverse mode AD currently not implemented for TC-mech module')
end

gdz = s.Grad(model.G.cells.centroids) * model.fluidModel.getGravityVector()';

% update state with AD-variables (required for state functions)
state = model.setProps(state, {'pressure', 'pressure_matrix', 'pressure_matrix2', 'xd'}, {p, pm, pm2, xd});
state0 = model.setProps(state0, {'pressure', 'pressure_matrix', 'pressure_matrix2', 'xd'}, {p0, pm0, pm20, xd0});
% set up properties
state = model.initStateFunctionContainers(state);


%% -------------------------------------------------------------------------
% Rock props
[poro_m, poro_m2, poro_f] = model.getProps(state, 'PoroelasticMatrixPoro', 'PoroelasticMatrix2Poro', 'PoroelasticFracturePoro');
[poro_m0, poro_m20, poro_f0] = model.getProps(state0, 'PoroelasticMatrixPoro', 'PoroelasticMatrix2Poro', 'PoroelasticFracturePoro');

%% -------------------------------------------------------------------------
% Fluid props
% evaluate fracture-water properties
bW     = f.bW(p);
bW0 = f.bW(p0);
rhoW   = bW.*f.rhoWS;
% rhoW on face, average of neighboring cells 
rhoWf  = s.faceAvg(rhoW);
mobW   = 1./f.muW(p);
dpW     = s.Grad(p) - rhoWf.*gdz;
% water upstream-index
upcw = (value(dpW)<=0);
 vW = - s.faceUpstr(upcw, mobW).*T_f.*dpW;
 bWvW = s.faceUpstr(upcw, bW).*vW; 

% evaluate macropores-water properties
bWm     = f_m.bW(pm);
bWm0 = f_m.bW(pm0);
rhoWm   = bWm.*f_m.rhoWS;
% rhoW on face, avarge of neighboring cells 
rhoWfm  = s_m.faceAvg(rhoWm);
mobWm   = 1./f_m.muW(pm);
dpWm     = s_m.Grad(pm) - rhoWfm.*gdz;
% water upstream-index
upcwm = (value(dpWm)<=0);
 vWm = - s_m.faceUpstr(upcwm, mobWm).*s_m.T.*dpWm;
 bWvWm = s_m.faceUpstr(upcwm, bWm).*vWm;   

% evaluate micropores-water properties
bWm2     = f_m2.bW(pm2);
bWm20 = f_m2.bW(pm20);
rhoWm2   = bWm2.*f_m2.rhoWS;
% rhoW on face, avarge of neighboring cells 
rhoWfm2  = s_m2.faceAvg(rhoWm);
mobWm2   = 1./f_m2.muW(pm2);           
dpWm2     = s_m2.Grad(pm2) - rhoWfm2.*gdz;
% water upstream-index
upcwm2 = (value(dpWm2)<=0);
 vWm2 = - s_m2.faceUpstr(upcwm2, mobWm2).*s_m2.T.*dpWm2; 
 bWvWm2 = s_m2.faceUpstr(upcwm2, bWm2).*vWm2; 

%%% interflux (flux between sets, e.g., fractures and macroporosity)
k13=model.mechModel.mech.k13;
k23=model.mechModel.mech.k23;
k12=model.mechModel.mech.k12;
T_f_m=k13;
T_f_m2=k23;
T_m_m2=k12;
vW_f_m = - s.faceUpstr(upcw, mobW).*T_f_m.*dpWm2;
bWvW_f_m = s.faceUpstr(upcw, bW).*vW_f_m; 
vW_f_m2 = - s.faceUpstr(upcw, mobW).*T_f_m2.*dpWm2;
bWvW_f_m2 = s.faceUpstr(upcw, bW).*vW_f_m2; 
vW_m_m2 = - s_m.faceUpstr(upcwm, mobWm).*T_m_m2.*dpWm2;
bWvW_m_m2 = s_m.faceUpstr(upcwm, bWm).*vW_m_m2; 

%% -------------------------------------------------------------------------
% Leakage calculation
Gamma12=model.mechModel.mech.Gamma12;
Gamma13=model.mechModel.mech.Gamma13;
Gamma23=model.mechModel.mech.Gamma23;
Leak12=Gamma12*(pm2-pm);
Leak13=Gamma13*(p-pm);
Leak23=Gamma23*(p-pm2);
Leak12 = model.G.cells.volumes.*Leak12;
Leak13 = model.G.cells.volumes.*Leak13;
Leak23 = model.G.cells.volumes.*Leak23;

% Output Additional Information
if model.outputFluxes
    state = model.storeFluxes(state, vW, vWm, vWm2);                    % original store function (state, vW, vo, vg) 
end

if model.extraStateOutput
    state = model.storebfactors(state, bW, bWm, bWm2);                  
    state = model.storeMobilities(state, mobW, mobWm, mobWm2);              
    state = model.storeUpstreamIndices(state, upcw, upcwm, upcwm2);     
    state.Leak12 = Leak12;
    state.Leak13 = Leak13;
    state.Leak23 = Leak23;
end

%% EQUATIONS ---------------------------------------------------------------
names = {'water', 'water_matrix', 'water_matrix2'};
types = {'cell', 'cell', 'cell'};
mob = {mobW, mobWm, mobWm2}; 
rho = {rhoW, rhoWm, rhoWm2}; 

% Upstream weight b factors (THAT CAN BE OMITTED) and multiply by interface fluxes to obtain the fluxes at standard conditions.  
% Fractures
eqs{1} = (model.G.cells.volumes./dt).*(poro_f.*bW - poro_f0.*bW0) + s.Div(bWvW) + s.Div(bWvW_f_m) + s.Div(bWvW_f_m2) + Leak13 + Leak23 ;  

% Macroporosity
eqs{2} = (model.G.cells.volumes./dt).*(poro_m.*bWm - poro_m0.*bWm0) + s.Div(bWvWm) + s.Div(bWvW_f_m) + s.Div(bWvW_m_m2) - Leak13 - Leak12;

% Microporosity
eqs{3} = (model.G.cells.volumes./dt).*( poro_m2.*bWm2 - poro_m20.*bWm20 ) + s.Div(bWvWm2) + s.Div(bWvW_f_m2) + s.Div(bWvW_m_m2) - Leak23 + Leak12; 

% Add BCs
% dummy saturation
sW = ones(model.G.cells.num, 1);
sWm = ones(model.G.cells.num, 1);
sWm2 = ones(model.G.cells.num, 1);

[eqs(1), state] = addBoundaryConditionsAndSources(model.fluidModel, ...
                                                  eqs(1), names(1), types(1), state, ...
                                                  {p,pm,pm2}, {sW}, mob(1), rho(1), ...
                                                  {}, {}, ...
                                                  drivingForces);
[eqs(2), state] = addBoundaryConditionsAndSources(model.fluidModel, ...
                                                  eqs(2), names(2), types(2), state, ...
                                                  {p,pm,pm2}, {sWm}, mob(2), rho(2), ...
                                                  {}, {}, ...
                                                  drivingForces);
[eqs(3), state] = addBoundaryConditionsAndSources(model.fluidModel, ...
                                                  eqs(3), names(3), types(3), state, ...
                                                  {p,pm,pm2}, {sWm2}, mob(3), rho(3), ...
                                                  {}, {}, ...
                                                  drivingForces);

% Well equations, for now, just implemented from fracture  
% these need to come from the fluid model
% [eqs(1), names, types, state.wellSol] = model.fluidModel.insertWellEquations(...
%                                                           eqs(1), names, types, wellSol0,...
%                                                           wellSol, wellVars, wellMap,...
%                                                           {p}, mob(1), rho(1), {}, {}, dt, opt); 

% [eqs(2), names, types, state.wellSol] = model.fluidModel.insertWellEquations(...
%                                                           eqs(2), names, types, wellSol0,...
%                                                           wellSol, wellVars, wellMap,...
%                                                           {pm}, mob(2), rho(2), {}, {}, dt, opt); 
% 
% [eqs(3), names, types, state.wellSol] = model.fluidModel.insertWellEquations(...
%                                                           eqs(3), names, types, wellSol0,...
%                                                           wellSol, wellVars, wellMap,...
%                                                           {pm2}, mob(3), rho(3), {}, {}, dt, opt); 


end
