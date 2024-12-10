clear all; close all; clc;

%  %  %  %  %  %
% DESCRIPTION: %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
% This code computes the poroelastic coefficients and relates them to the
% imposed stress tensor. The model has three vertical subsets of cracks
% that are either connected (forming a single connected set) or isolated
% (forming three isolated sets). Both possibilities are ploted. Functions
% designed to compute dry excess compliances (effective medium theory) and 
% then poroelastic coefficients (Shapiro and Kachanov, 1997) are utilised.
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %


%%%_required_parameters_%%%

m=100;                        % number of pores in a subset 
y=0.01;                       % aspect ratio of pores in each subset
theta_p1=10;                  % azimuthal orientation of a first vertical subset
theta_p2=130;                 % azimuthal orientation of a second vertical subset
theta_p3=250;                 % azimuthal orientation of a third vertical subset
Kf=1/0.45;                    % fluid compressibility
E=87;                         % Young's modulus of a background solid
v=0.11;                       % Poisson's ratio of a background solid

[H_c, phi, Kf, K] = fun_emt_eachpore(m, y, theta_p1, theta_p2, theta_p3, Kf, E, v);     % compute dry excess compliance and volume fraction for each pore, also provide bulk modulus of fluid and solid background 
[S, B] = fun_poroelast_singleset(H_c, phi, Kf, K)                                      % compute poroelast. coefficients for a single connected set composed of three subsets
[S_p1, S_p2, S_p3, B_p1, B_p2, B_p3] = fun_poroelast_isolatedsets(H_c, phi, Kf, K)     % compute analogous poroelast. coefficients for the case of isolated subsets

sigma11 = linspace(1,100,100);                                                          % uniaxial stress - the above functions are general enough to support any stress tensor choice

for i=1:length(sigma11)
pf(i)=-(1/3).*B(1,1).*sigma11(i);
pf_p1(i)=-(1/3).*B_p1(1,1).*sigma11(i);
pf_p2(i)=-(1/3).*B_p2(1,1).*sigma11(i);
pf_p3(i)=-(1/3).*B_p3(1,1).*sigma11(i);
zeta(i)=(1/3).*S.*B(1,1).*sigma11(i);
zeta_p1(i)=(1/3).*S_p1.*B_p1(1,1).*sigma11(i); 
zeta_p2(i)=(1/3).*S_p2.*B_p2(1,1).*sigma11(i);
zeta_p3(i)=(1/3).*S_p3.*B_p3(1,1).*sigma11(i); 
end

%% PLOTS

figure; 
plot(sigma11,-pf,'black'); hold on;
plot(sigma11,-pf_p1,'black--'); hold on;
plot(sigma11,-pf_p2,'black-.'); hold on;
plot(sigma11,-pf_p3,'black:'); hold on;
plot(sigma11,-(pf_p1+pf_p2+pf_p3)/3,'black'); hold on;                     % check: this line should be equal to -pf
legend('p_f^{connected}', 'p_f^{(1)}', 'p_f^{(2)}', 'p_f^{(3)}'); 
grid on; axis square;
ylim([0 80])
xlabel('\sigma_{11} [MPa]'); ylabel('pore pressure [MPa]');

figure; 
plot(sigma11,0.001*zeta,'black'); hold on;                                 % sigma in MPa but S is in 1/GPa multiplication by 0.001 neded for dimensionless result
plot(sigma11,0.001*zeta_p1,'black--'); hold on;
plot(sigma11,0.001*zeta_p2,'black-.'); hold on;
plot(sigma11,0.001*zeta_p3,'black:'); hold on;
plot(sigma11,0.001*(zeta_p1+zeta_p2+zeta_p3),'black'); hold on;            % check: this line should be equal to 0.001*zeta
legend('\zeta^{connected}', '\zeta^{(1)}', '\zeta^{(2)}', '\zeta^{(3)}'); 
grid on; axis square;
xlabel('\sigma_{11} [MPa]'); ylabel('fluid content change');

