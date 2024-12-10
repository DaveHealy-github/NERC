function [alpha_m, alpha_m2, alpha_f, invM_m, invM_m2, invM_f, invM12, invM13, invM23] = TriCalcIsotropicCoefficientModels(coefficients, model)

% DESCRIPTION: 
%   Computes poroelastic coefficients from equations (24)-(25) in 
%   Adamus et al. (2023)

% PARAMETERS:
%   coefficients - class instance
%   model 

% RETURNS:
%   alpha_f    - Biot coefficient of the fracture phase
%   alpha_m    - Biot coefficient of the macropore phase
%   alpha_m2   - Biot coefficient of the micropore phase
%   invN_f     - Storage modulus of the fracture phase
%   invN_m     - Storage modulus of the macropore phase
%   invN_m2    - Storage modulus of the micropore phase
%   invM13     - Interset storage modulus between fractures and macropores 
%   invM23     - Interset storage modulus between fractures and micropores
%   invM12     - Interset storage modulus between macropores and micropores

d = model.G.griddim;
if d == 3
    rank2 = @(x) bsxfun(@times, x, [1, 1, 1, 0, 0, 0]);
else
    rank2 = @(x) bsxfun(@times, x, [1, 1, 0]);
end
 
K_f = model.mech.E_f./(3.*(1-2.*model.mech.nu_f));    % bulk modulus of a fracture phase
K_m = model.mech.E_m./(3.*(1-2.*model.mech.nu_m));    % bulk modulus of a macropore phase
K_m2= model.mech.E_m2./(3.*(1-2.*model.mech.nu_m2));  % bulk modulus of a micropore phase
K=1./((1./K_m).*model.rock_matrix.vol_fraction + (1./K_m2).*model.rock_matrix2.vol_fraction + (1./K_f).*model.rock.vol_fraction); % bulk modulus of upscaled medium

a22=model.mech.a22;
a33=model.mech.a33;
a44=model.mech.a44;
b1=model.mech.b1;
b2=model.mech.b2;
b3=model.mech.b3;

alpha_m  = 3.*b1.*K;       
alpha_m2 = 3.*b2.*K;       
alpha_f  = 3.*b3.*K;       
invM_m = a22-K.*b1.^2;     
invM_m2 = a33-K.*b2.^2;    
invM_f  = a44-K.*b3.^2;    
invM12  = model.mech.a23;  
invM13  = model.mech.a24;  
invM23  = model.mech.a34;  

alpha_m  = rank2(alpha_m);
alpha_m2 = rank2(alpha_m2);
alpha_f  = rank2(alpha_f);
end

