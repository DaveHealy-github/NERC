function [eqs, names, types] = equationsTCPoroMechanics(x, model, fluidp, fluidp_matrix, fluidp_matrix2)

% DESCRIPTION: 
%   Assembles the residual for triple-permeability momentum balance. The
%   function takes fluid input, given by the pore pressures.

% PARAMETERS:
%   x                 - Displacement
%   model             - Model class instance that is used
%   fluidp            - Fluid pressure of the fracture phase
%   fluidp_matrix     - Fluid pressure of the macropore phase
%   fluidp_matrix2    - Fluid pressure of the micropore phase

% RETURNS:
%   eqs   - The residual values as ADI variables (that is with the Jacobian)
%           if the inputs were also ADI.
%   names - The name of each equations
%   types - The type of each equations

    
    G = model.G;
    d = G.griddim;
    s = model.operators;
    cM = model.constitutive_coefficients_object;
    alpha_f = repmat(cM.alpha_f(:,1:d), 1, size(s.ovol_div,2)/d);
    alpha_f = alpha_f(:, ~s.isdirdofs);
    alpha_m = repmat(cM.alpha_m(:,1:d), 1, size(s.ovol_div,2)/d);
    alpha_m = alpha_m(:, ~s.isdirdofs);
    alpha_m2= repmat(cM.alpha_m2(:,1:d), 1, size(s.ovol_div,2)/d);
    alpha_m2= alpha_m2(:, ~s.isdirdofs);

    eqs{1} = s.A * x - (s.gradP .* alpha_m') * fluidp_matrix - (s.gradP .*  alpha_f') * fluidp - (s.gradP .*  alpha_m2') * fluidp_matrix2 - s.rhs;

    % normalization constant
    fac =  1 / (1e6 * mean(G.cells.volumes)); 
    eqs{1} = eqs{1} * fac;
    names = {'disp'};
    types = {'disp_dofs'};

end