function [E, K, shear, nu] = RV_bounds_choose(E_m, E_m2, E_f, nu_m, nu_m2, nu_f, vol_m, vol_m2, vol_f, type)

% DESCRIPTION:
%   Function to calculate the upscaled properties of a three phase, isotropic
%   material according to Reuss (lower) bound or Voigt (upper) bound.
%
% PARAMETERS:
%   E_m     - matrix Young's modulus
%   E_f     - fracture Young's modulus
%   nu_m    - matrix Poisson's ratio
%   nu_f    - fracture Poisson's ratio
%   vol_m   - macropore phase volume fraction
%   vol_m2  - micropore phase volume fraction
%   vol_f   - fracture phase volume fraction
%   type    - lower or upper bound
%
% RETURNS:
%   E      - upscaled Young's modulus (2D)
%   K      - upscaled bulk modulus
%   shear  - upscaled shear modulus
%   nu     - upscaled Poisson's modulus (2D)

% Shear moduli
shear_f = E_f/(2*(1+nu_f));
shear_m = E_m/(2*(1+nu_m));
shear_m2= E_m2/(2*(1+nu_m2));

% Bulk moduli, with plane-stress
K_f = E_f/(3*(1-2*nu_f)); 
K_m = E_m/(3*(1-2*nu_m)); 
K_m2= E_m2/(3*(1-2*nu_m2));

% Calculate 
if strcmp(type, 'upper') % Voigt
    K = K_m*vol_m + K_m2*vol_m2 + K_f*vol_f;
    shear = shear_m*vol_m + shear_m2*vol_m2 + shear_f*vol_f;
elseif strcmp(type, 'lower') % Reuss
    K = 1/((1/K_m)*vol_m + (1/K_m2)*vol_m2 + (1/K_f)*vol_f);
    shear = 1/((1/shear_m)*vol_m + (1/shear_m2)*vol_m2 + (1/shear_f)*vol_f);
end

K_2d = 9*K*shear/(3*K+4*shear);
nu = (K_2d - shear)/(K_2d + shear);
E = K_2d*2*(1-nu);
end
