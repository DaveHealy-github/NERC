clear all; close all; clc;

%  %  %  %  %  %
% DESCRIPTION: %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
% This code simulates a uniaxial single-phase consolidation problem,    %
% such that the bottom boundary is fixed, and the left, right and top   %
% boundaries admit only vertical displacement. All of the boundaries    %
% are no flux boundaries, apart from the top boundary which is exposed  %
% to time-dependent drainage.                                           %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %


%%% Load required modules %%%
mrstModule add triple-porosity dual-continuum-mech ad-core ad-mechanics dual-porosity ad-props vemmech

%%% Setup default options %%%
opt = struct('cartDims'            , [20, 20], ...
             'L'                  , [2, 2], ...
             'fluid_model'        , 'water', ...
             'verbose'            , false);

%%% Setup Grid %%%
G = cartGrid(opt.cartDims, opt.L);
G = computeGeometry(G);
 plotGrid(G);


%% Setup rock parameters (fluid description and volume fractions - third and forth column of Table 2 in Adamus et al (2023))
% multipliers (>0)
x_coeff=0.8; y_coeff=0.5; z_coeff=10;

% macropore phase (1)
perm_matrix = 0*milli*darcy*ones(G.cartDims);                      % k(1)
rock_matrix = struct('perm', reshape(perm_matrix, [], 1), ...      
              'poro', ones(G.cells.num, 1)*0.095,...               % phi(1)
              'vol_fraction', ones(G.cells.num, 1)*0.49525);       % v(1)

% micropore phase (2)
 perm_matrix2 = 0*milli*darcy*ones(G.cartDims);                    % k(2)
 rock_matrix2 = struct('perm', reshape(perm_matrix2, [], 1), ...
               'poro', ones(G.cells.num, 1)*1*0.05,...             % phi(2)
               'vol_fraction', ones(G.cells.num, 1)*0.49525);      % v(2)

 % fracture phase (3)
perm_fracture = z_coeff*milli*darcy*ones(G.cartDims);              % k(3)
rock_fracture = struct('perm', reshape(perm_fracture, [], 1), ...
              'poro', ones(G.cells.num, 1)*0.009,...               % phi(3)
              'vol_fraction', ones(G.cells.num, 1)-...
              rock_matrix.vol_fraction-rock_matrix2.vol_fraction); % v(3)    

                        
%% Setup fluid model and choose interflow parameters (third column of Table 2 in Adamus et al (2023))
mu=1*centi*poise;                                  % fluid viscosity 
L=0.1;                                             % fracture spacing (concentration parameter)
delta=pi^2;                                        % shape factor

k12=0;                                             % permeability at the interset (macropores and micropores)
k13=(1/z_coeff)*milli*darcy;                       % permeability at the interset (macropores and fractures)
k23=(1/z_coeff)*milli*darcy;                       % permeability at the interset (micropores and fractures)

Gamma12=0;                                         % leakage coeff. (between macropores and micropores)
Gamma13=(delta*(milli*darcy/z_coeff))/(L^2*mu);    % leakage coeff. (between macropores and fractures)
Gamma23=(delta*(milli*darcy/z_coeff))/(L^2*mu);    % leakage coeff. (between micropores and fractures)

fluid = initSimpleADIFluid('phases', 'W', 'mu', mu, 'rho', 1000*kilogram/meter^3); 
fluid_fracture = fluid; fluid_matrix = fluid; fluid_matrix2 = fluid;


%% Setup rock parameters (stiffnesses and poroel. coeff. - first and second column of Table 2)
a22=0.0993*1e-9;     % storage modulus (macropore phase)
a33=x_coeff*a22;     % storage modulus (micropore phase)
a44=0.145*1e-9;      % storage modulus (fracture phase)
a23=0;               % interset storage modulus (between macropores and micropores)
a24=1.35e-12;        % interset storage modulus (between macropores and fractures)
a34=1.351e-12;       % interset storage modulus (between micropores and fractures)
b1=0.0253*1e-9;      % poroelastic expansion (macropore phase)
b2=y_coeff*b1;       % poroelastic expansion (micropore phase)
b3=0.049*1e-9;       % poroelastic expansion (fracture phase)

nu = 0.15;           % Poisson's ratio (intact solid)
E_m = 36e9;          % Young's modulus (macropore phase)
nu_m = nu;           % Poisson's ratio (macropore phase)
E_m2 = 50e9;         % Young's modulus (micropore phase)
nu_m2 = nu;          % Poisson's ratio (micropore phase)
E_f = 0.15E9;        % Young's modulus (fracture phase)
nu_f = 0.12;         % Poisson's ratio (fracture phase)
K_s = 37E9;          % Bulk modulus (intact solid)

[E,K,shear,nu] = RV_bounds_choose(E_m, E_m2, E_f, nu_m, nu_m2, nu_f, rock_matrix.vol_fraction(1),...
              rock_matrix2.vol_fraction(1), rock_fracture.vol_fraction(1), 'lower'); % stiffnesses of an upscaled rock with Reuss bound

a22 = repmat(a22, G.cells.num, 1);   
a33 = repmat(a33, G.cells.num, 1);
a44 = repmat(a44, G.cells.num, 1);
a23 = repmat(a23, G.cells.num, 1);
a24 = repmat(a24, G.cells.num, 1);
a34 = repmat(a34, G.cells.num, 1);
b1 = repmat(b1, G.cells.num, 1);
b2 = repmat(b2, G.cells.num, 1);
b3 = repmat(b3, G.cells.num, 1);

E = repmat(E, G.cells.num, 1);   
nu = repmat(nu, G.cells.num, 1);
E_m = repmat(E_m, G.cells.num, 1);
nu_m = repmat(nu_m, G.cells.num, 1);
E_m2 = repmat(E_m2, G.cells.num, 1);
nu_m2 = repmat(nu_m2, G.cells.num, 1);
E_f = repmat(E_f, G.cells.num, 1);
nu_f = repmat(nu_f, G.cells.num, 1);
K_s = repmat(K_s, G.cells.num, 1);


%% Setup boundary conditions for mechanics
% we first want to create a structure 'bc', which we can fudge by initialising the bc's using pside (see other modules). 
oside = {'WEST', 'EAST', 'SOUTH', 'NORTH'};
bc = cell(4,1);
for i = 1:numel(oside)
    bc{i} = pside([], G, oside{i}, 0);
    bc{i} = rmfield(bc{i}, 'type'); 
    bc{i} = rmfield(bc{i}, 'sat');    
end

% Displacement BCs
% Find the nodes for the different sides and set the boundaray conditions for elasticity.
for i = 1 : 4
    inodes = mcolon(G.faces.nodePos(bc{i}.face), G.faces.nodePos(bc{i}.face + 1) - 1);
    nodes = unique(G.faces.nodes(inodes));
    disp_bc = struct('nodes'   , nodes,      ...
                     'uu'      , 0,          ...
                     'faces'   , bc{i}.face, ...
                     'uu_face' , 0,          ...          
                     'mask'    , true(numel(nodes), G.griddim));
    bc{i}.el_bc = struct('disp_bc', disp_bc, 'force_bc', []);
end
bcdisp_zero = @(x) x*0.0; % Boundary displacement function set to zero.
bc_el_sides{1} = bc{1}; 
bc_el_sides{1}.el_bc.disp_bc.mask(:, 2) = false;   % x fixed, y free
bc_el_sides{2} = bc{2}; 
bc_el_sides{2}.el_bc.disp_bc.mask(:, 2) = false;   % x fixed, y free
bc_el_sides{3} = bc{3}; 
bc_el_sides{3}.el_bc.disp_bc.mask(:, :) = true;    % x fixed, y fixed
bc_el_sides{4} = bc{4}; 
bc_el_sides{4}.el_bc.disp_bc.mask(:, 2) = false;   % x fixed, y free

% collect the displacement boundary conditions
nodes = [];
faces = [];
mask = [];
for i = 1 : numel(bc)
    if(~isempty(bc_el_sides{i}))
        nodes = [nodes; bc_el_sides{i}.el_bc.disp_bc.nodes]; 
        faces = [faces; bc_el_sides{i}.el_bc.disp_bc.faces]; 
        mask  = [mask; bc_el_sides{i}.el_bc.disp_bc.mask]; 
    end
end

disp_node = bcdisp_zero(G.nodes.coords(nodes, :));
disp_faces = bcdisp_zero(G.faces.centroids(faces, :)); 
disp_bc = struct('nodes', nodes, 'uu', disp_node, 'faces', faces, 'uu_face', disp_faces, 'mask', mask); 

% Force BCs
facesf = bc{4}.face;
force = 1e6;             % compressional force 
force_bc = struct('faces', facesf, 'force', force*ones(G.cartDims(1),1)*[0 -1]);
el_bc = struct('disp_bc', disp_bc, 'force_bc', force_bc);


%% Gravity
% the gravity in this option affects only the fluid behaviour
gravity off;
    

%% Setup load for mechanics
% in this example we do not impose any volumetric force
load = @(x) (0*x);


%% Gather all parameters in a struct that will be used in mechanical and combined fluid-mechanical models
mech = struct('a22',a22 ,'a33', a33, 'a44', a44, 'a23', a23, 'a24', a24, 'a34', a34, 'b1', b1, 'b2', b2, 'b3', b3, ...
              'E', E, 'nu', nu, 'E_m', E_m, 'nu_m', nu_m, 'E_m2', E_m2, 'nu_m2', nu_m2, 'E_f', E_f, 'nu_f', nu_f,...
              'K_s', K_s, 'el_bc', el_bc, 'load', load, 'x_coeff', x_coeff, 'y_coeff', y_coeff, ...
              'k12', k12, 'k13', k13, 'k23', k23, 'Gamma12', Gamma12, 'Gamma13', Gamma13, 'Gamma23', Gamma23);


%% Setup fully coupled and fixed stress splitting models
fullycoupledOptions = {'verbose', opt.verbose};
TC_model = TriContMechWaterModel(G, {rock_fracture, rock_matrix, rock_matrix2}, {fluid_fracture, fluid_matrix, fluid_matrix2}, mech, fullycoupledOptions{:});
TC_model.transfer_model_object = TriTransferFunction();
TC_model = TC_model.validateModel();


%% Setup initial state and fluid BCs
pressure = zeros(G.cells.num,1);
state0 = struct('pressure', pressure, 'pressure_matrix', pressure, 'pressure_matrix2', pressure, 's', ones(G.cells.num, 1), 'swm', ones(G.cells.num, 1), 'swm2', ones(G.cells.num, 1)); 

% need to initiate the fluid bc's
bc_f0 = fluxside([], G, 'WEST', 0, 'sat', 1);
bc_f0 = fluxside(bc_f0, G, 'EAST', 0,'sat', 1);
bc_f0 = fluxside(bc_f0, G, 'SOUTH', 0, 'sat', 1);
bc_f0 = fluxside(bc_f0, G, 'NORTH', 0, 'sat', 1);

%% Simulate 
t=100;
time = [0, logspace(-6,5,t)];
dt = diff(time);
[p_m, p_m2, p_f, ~, states] = simTC_mech(state0, dt, TC_model, bc_f0);   

%% Plot results
figure
semilogx(time, 1e-6*sum(p_m,1)./G.cells.num, '-blue', 'linewidth', 2); hold on
semilogx(time, 1e-6*sum(p_m2,1)./G.cells.num, '-red', 'linewidth', 2); hold on
semilogx(time, 1e-6*sum(p_f,1)./G.cells.num, '-black', 'linewidth', 2); hold on
axis square; grid on
xlim([0.9999e-5 1.00e5]); ylim([0 0.5])
xlabel('time [s]'); ylabel('average pressure [MPa]')
legend('macropores', 'micropores', 'fractures')
title('Triple-porosity: average pressure vs. time')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DELTE THE COMMENT SIGNS TO SEE MORE CONSOLIDATION PLOTS FOR:     %%%
%%% LOCAL PRESSURES, FLUID CONTENT CHANGES, AND STRAINS VERSUS TIME  %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{

%% Additional data for additional plots
p_m_bottom=p_m((G.cartDims(1,2)*G.cartDims(1,1)-G.cartDims(1,1)+1):(G.cartDims(1,2)*G.cartDims(1,1)),:);
p_m2_bottom=p_m2((G.cartDims(1,2)*G.cartDims(1,1)-G.cartDims(1,1)+1):(G.cartDims(1,2)*G.cartDims(1,1)),:);
p_f_bottom=p_f((G.cartDims(1,2)*G.cartDims(1,1)-G.cartDims(1,1)+1):(G.cartDims(1,2)*G.cartDims(1,1)),:);
p_m_top=p_m(1:G.cartDims(1,1),:);
p_m2_top=p_m2(1:G.cartDims(1,1),:);
p_f_top=p_f(1:G.cartDims(1,1),:);

if mod(G.cartDims(1,1),2) == 1 && mod(G.cartDims(1,2),2) == 0
p_m_middle=p_m(0.5*G.cartDims(1,2)*G.cartDims(1,1)-0.5*G.cartDims(1,1)+0.5:0.5*G.cartDims(1,2)*G.cartDims(1,1)+0.5*G.cartDims(1,1)-0.5,:);
p_m2_middle=p_m2(0.5*G.cartDims(1,2)*G.cartDims(1,1)-0.5*G.cartDims(1,1)+0.5:0.5*G.cartDims(1,2)*G.cartDims(1,1)+0.5*G.cartDims(1,1)-0.5,:);
p_f_middle=p_f(0.5*G.cartDims(1,2)*G.cartDims(1,1)-0.5*G.cartDims(1,1)+0.5:0.5*G.cartDims(1,2)*G.cartDims(1,1)+0.5*G.cartDims(1,1)-0.5,:);
else    
p_m_middle=p_m(0.5*G.cartDims(1,2)*G.cartDims(1,1)-0.5*G.cartDims(1,1)+1:0.5*G.cartDims(1,2)*G.cartDims(1,1)+0.5*G.cartDims(1,1),:);
p_m2_middle=p_m2(0.5*G.cartDims(1,2)*G.cartDims(1,1)-0.5*G.cartDims(1,1)+1:0.5*G.cartDims(1,2)*G.cartDims(1,1)+0.5*G.cartDims(1,1),:);
p_f_middle=p_f(0.5*G.cartDims(1,2)*G.cartDims(1,1)-0.5*G.cartDims(1,1)+1:0.5*G.cartDims(1,2)*G.cartDims(1,1)+0.5*G.cartDims(1,1),:);
end

alpha_f=TC_model.mechModel.constitutive_coefficients_object.alpha_f;       
alpha_m=TC_model.mechModel.constitutive_coefficients_object.alpha_m; 
alpha_m2=TC_model.mechModel.constitutive_coefficients_object.alpha_m2; 
invM_f=TC_model.mechModel.constitutive_coefficients_object.invM_f;
invM_m=TC_model.mechModel.constitutive_coefficients_object.invM_m; 
invM_m2=TC_model.mechModel.constitutive_coefficients_object.invM_m2;
invM13=TC_model.mechModel.constitutive_coefficients_object.invM13;
invM23=TC_model.mechModel.constitutive_coefficients_object.invM23; 
invM12=TC_model.mechModel.constitutive_coefficients_object.invM12;

for i=1:(t+1)
zeta_f(:,i)=(TC_model.mechModel.operators.div.*alpha_f(1,1))*states{i,1}.xd./(G.cells.volumes) + invM_f(1,:).*p_f(:,i) + invM13(1,:).*p_m(:,i) + invM23(1,:).*p_m2(:,i);
zeta_m(:,i)=(TC_model.mechModel.operators.div.*alpha_m(1,1))*states{i,1}.xd./(G.cells.volumes) + invM_m(1,:).*p_m(:,i) + invM13(1,:).*p_f(:,i) + invM12(1,:).*p_m2(:,i);
zeta_m2(:,i)=(TC_model.mechModel.operators.div.*alpha_m2(1,1))*states{i,1}.xd./(G.cells.volumes) + invM_m2(1,:).*p_m2(:,i) + invM12(1,:).*p_m(:,i) + invM23(1,:).*p_f(:,i);
strain(:,i)=(TC_model.mechModel.operators.div.*1)*states{i,1}.xd./(G.cells.volumes);
end


%% Additional plots
figure
semilogx(time, 1e-6*sum(p_m_top,1)./G.cartDims(1,1), '.blue', 'linewidth', 2); hold on
semilogx(time, 1e-6*sum(p_m2_top,1)./G.cartDims(1,1), '.red', 'linewidth', 2); hold on
semilogx(time, 1e-6*sum(p_f_top,1)./G.cartDims(1,1), '.black', 'linewidth', 2); hold on
semilogx(time, 1e-6*sum(p_m_middle,1)./(G.cartDims(1,1)), '-blue', 'linewidth', 2); hold on
semilogx(time, 1e-6*sum(p_m2_middle,1)./(G.cartDims(1,1)), '-red', 'linewidth', 2); hold on
semilogx(time, 1e-6*sum(p_f_middle,1)./(G.cartDims(1,1)), '-black', 'linewidth', 2); hold on
semilogx(time, 1e-6*sum(p_m_bottom,1)./G.cartDims(1,1), '--blue', 'linewidth', 2); hold on
semilogx(time, 1e-6*sum(p_m2_bottom,1)./G.cartDims(1,1), '--red', 'linewidth', 2); hold on
semilogx(time, 1e-6*sum(p_f_bottom,1)./G.cartDims(1,1), '--black', 'linewidth', 2); hold on
grid on; axis square
xlim([0.9999e-5 1.00e5]); ylim([0 0.5])
xlabel('time [s]'); ylabel('local pressure [MPa]')
legend('macro. top','micro. top', 'frac. top', 'macro. mid.','micro. mid.', 'frac. mid.', 'macro. bot.', 'micro. bot.', 'frac. bot.')
title('Triple-porosity: local pressure vs. time')

figure
semilogx(time, sum(zeta_m,1)./G.cells.num, '-', 'linewidth', 2); hold on
semilogx(time, sum(zeta_m2,1)./G.cells.num, '--', 'linewidth', 2); hold on
semilogx(time, sum(zeta_f,1)./G.cells.num, '-black', 'linewidth', 2); hold on
axis square; grid on
xlim([0.9999e-5 1.00e5])
xlabel('time [s]'); ylabel('average fluid content change')
legend('macropores', 'micropores', 'fractures')
title('Triple-porosity: fluid content vs. time')

figure
semilogx(time, sum(strain,1)./G.cells.num, '-', 'linewidth', 2)
axis square; grid on
xlim([0.9999e-5 1.00e5])
xlabel('time [s]'); ylabel('average strain')
title('Triple-porosity: strain vs. time')

%}