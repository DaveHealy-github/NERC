clear all; close all; clc;

%  %  %  %  %  %  
%  DESCRIPTION %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  % 
%  This code compares two micromechanical approaches that predict       %
%  the effect of fluid (difference between undrained and drained        %
%  compliances); without loss of generalisation, we assume one set      %
%  composed of three connected subsets containing m pores each.         %
%  Using this code, Figure 6 and Table 1 of Adamus et al. (2023),       %
%  ''Multi-porous extension of anisotropic poroelasticity: linkage      %
%  with micormechanics'', JGR, can be generated. User MUST choose the   % 
%  microstrucural setting/case, the range of the aspect ratio,          %
%  and the stiffness of the background (e.g., Berea sandstone).         %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %


%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %      
%  CASES OF MICROSTRUCTURE FROM TABLE 1 (shape, orientation, size)      %                                                                      %
%  i-identical, s-slightly varying, p-nonrand pattern, r-random         %
%                                                                       %
%  1) random shape, orientation, size                   (r,r,r)         %
%  2) patterned shape, orientation, size                (p,p,p)         %
%  3) slightly-varying shape, orientation, size         (s,s,s)         %
%  4) random shapes, identical orientation, size        (r,i,i)         %
%  5) random orientation, identical shape, size         (i,r,i)         %
%  6) random size, identical shape, orientation         (i,i,r)         %
%  7) patterned shape, identical orientation, size      (p,i,i)         %
%  8) patterned orientation, identical shape, size      (i,p,i)         %
%  9) patterned size, identical shape, orientation      (i,i,p)         %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %


%%%_set_parameters_%%%

m=linspace(1,1000,1000);      % changing number of pores in a set
k=1;                        % number of iterations for each m pores scenario (can be set to 1)
ym=2;                       % maximum aspect ratio considered

for jj=1:length(m)
x=3.*m(jj);                 % number of pores in a set (changes with m)
Kf=1/0.45*ones(x,k);        % fluid compressibility 
E=87*ones(x,k);             % Young's modulus of a background solid
v=0.11*ones(x,k);           % Poisson's ratio of a background solid
l=E.*v./((1+v).*(1-2.*v));  % Lame' parameter of a background solid
G=E./(2+2.*v);              % rigidity of a backround solid
K=E./(3.*(1-2.*v));         % bulk modulus of a backround solid

for j=1:k
for i=1:x
Cb(:,:,i,j)= [ l(i,j)+2.*G(i,j), l(i,j), l(i,j), 0, 0, 0;
     l(i,j),  l(i,j)+2.*G(i,j), l(i,j), 0, 0, 0;
     l(i,j), l(i,j),  l(i,j)+2.*G(i,j), 0, 0, 0;
     0, 0, 0, G(i,j), 0, 0;
     0, 0, 0, 0, G(i,j), 0; 
     0, 0, 0, 0, 0, G(i,j); ]; % background stiffness
end
end
F=[1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 2 0 0;
    0 0 0 0 2 0;   % correction matrix due to double-dot product (P:Cb or Cb:e)
    0 0 0 0 0 2;];


%%%_choose_pore_shapes_by_UNCOMMENTING_ONE_OUT_OF_FOUR_opions_in_a_box_below_%%%

x0=rand(1,x).*ones(x,k);      % x_n are random numbers used for shape, orient, size calculations
x1=rand(1,x).*ones(x/3,k); 
x2=rand(1,x).*ones(x/3,k); 
x3=rand(1,x).*ones(x/3,k);

for j=1:k
for i=1:x

% UNCOMMENT ONE !!! %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  y(i,j)=ym.*rand;                               % random shape
%  y=ym.*x0;                                      % identical shape
%  y(i,j)=ym.*x0(i,j).*(1+0.2.*rand-0.1);         % almost identical shape
%  y=ym.*[x1; x2; x3];                            % pattern for each subset 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%_Eshelby matrix computation_%%%

if y(i,j)<1
g=(y./((1-y.^2).^(3./2))).*(acos(y)-y.*sqrt(1-y.^2));
else
g=(y./((y.^2-1).^(3./2))).*(y.*(y.^2-1).^(1./2)-acosh(y));
end

 e11=(-3.*y.^2)./(8.*(1-v(i,j)).*(1-y.^2))+(g./(4.*(1-v(i,j)))).*(1-2.*v(i,j)+(9)./(4.*(1-y.^2)));
 e33=(1./(1-v(i,j))).*(2-v(i,j)-1./(1-y.^2))+(g./(2.*(1-v(i,j)))).*(-4+2.*v(i,j)+3./(1-y.^2));
 e12=(1./(8.*(1-v(i,j)))).*(1-1./(1-y.^2))+(g./(16.*(1-v(i,j)))).*(-4+8.*v(i,j)+3./(1-y.^2));
 e13=(y.^2)./(2.*(1-v(i,j)).*(1-y.^2))-(g./(4.*(1-v(i,j)))).*(1-2.*v(i,j)+(3.*y.^2)./((1-y.^2)));
 e31=(1./(2.*(1-v(i,j)))).*(-1+2.*v(i,j)+1./(1-y.^2))+(g./(4.*(1-v(i,j)))).*(2-4.*v(i,j)-3./(1-y.^2));
 e66=(-1.*y.^2)./(8.*(1-v(i,j)).*(1-y.^2))+(g./(16.*(1-v(i,j)))).*(4-8.*v(i,j)+(3)./(1.*(1-y.^2)));
 e44=(1./(4.*(1-v(i,j)))).*(1-2.*v(i,j)+(1+y.^2)./(1-y.^2))-(g./(8.*(1-v(i,j)))).*(1-2.*v(i,j)+3.*(1+y.^2)./(1-y.^2));
 e22=e11; e21=e12; e23=e13; e32=e31; e55=e44;
 e(:,:,i,j)=[e11(i,j) e12(i,j) e13(i,j) 0 0 0; e21(i,j) e22(i,j) e23(i,j) 0 0 0; e31(i,j) e32(i,j) e33(i,j) 0 0 0;
             0 0 0 e44(i,j) 0 0;    0 0 0 0 e55(i,j) 0;     0 0 0 0 0 e66(i,j); ];
end
end


%%%_choose_pore_orientations_by_UNCOMMENTING_ONE_OUT_OF_FOUR_opions_in_a_box_below_%%%

x4=rand(3,1,k).*ones(3,x,k); x4=x4./vecnorm(x4); x5=360*rand(1,k).*ones(x,k);   
x6=rand(3,1,k).*ones(3,x/3,k); x6=x6./vecnorm(x6); 
x7=rand(3,1,k).*ones(3,x/3,k); x7=x7./vecnorm(x7);
x8=rand(3,1,k).*ones(3,x/3,k); x8=x8./vecnorm(x8); x9=cat(2,x6,x7,x8);
x10=360.*rand(1,k).*ones(x/3,k); x11=360.*rand(1,k).*ones(x/3,k); x12=360.*rand(1,k).*ones(x/3,k); 

for j=1:k
for i=1:x

% UNCOMMENT ONE !!! %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    w(:,i,j)=rand(3,1,1); w(:,i,j)=w(:,i,j)/norm(w(:,i,j)); theta(i,j)=rand*360;              % random orientations  (w-symmetry axis, theta-rotation about w)   
%    w=x4;   theta=x5;                                                                         % identical orientations of every pore                    
%    w=x4; theta(i,j)=x5(i,j)+rand*20-10;                                                      % almost identical (up to 10 degrees difference)                                  
%    w=x9; theta=[x10;x11;x12;];                                                               % patterns for each subset (identical orientations) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                                                                                     
bw(:,:,i,j) = [0, -w(3,i,j), w(2,i,j); w(3,i,j), 0, -w(1,i,j); -w(2,i,j), w(1,i,j), 0];
A(:,:,i,j) = eye(3) + sind(theta(i,j))*bw(:,:,i,j) + (1-cosd(theta(i,j)))*bw(:,:,i,j)*bw(:,:,i,j);

A11(i,j)=A(1,1,i,j); A12(i,j)=A(1,2,i,j); A13(i,j)=A(1,3,i,j); A21(i,j)=A(2,1,i,j); A22(i,j)=A(2,2,i,j); A23(i,j)=A(2,3,i,j); A31(i,j)=A(3,1,i,j); A32(i,j)=A(3,2,i,j); A33(i,j)=A(3,3,i,j);   
MA(:,:,i,j)= [A11(i,j).*A11(i,j) A12(i,j).*A12(i,j) A13(i,j).*A13(i,j) A12(i,j).*A13(i,j) A11(i,j).*A13(i,j) A11(i,j).*A12(i,j);
     A21(i,j).*A21(i,j) A22(i,j).*A22(i,j) A23(i,j).*A23(i,j) A22(i,j).*A23(i,j) A21(i,j).*A23(i,j) A21(i,j).*A22(i,j);
     A31(i,j).*A31(i,j) A32(i,j).*A32(i,j) A33(i,j).*A33(i,j) A32(i,j).*A33(i,j) A31(i,j).*A33(i,j) A31(i,j).*A32(i,j);
     2.*A21(i,j).*A31(i,j) 2.*A22(i,j).*A32(i,j) 2.*A23(i,j).*A33(i,j) A22(i,j).*A33(i,j)+A23(i,j).*A32(i,j) A21(i,j).*A33(i,j)+A23(i,j).*A31(i,j) A21(i,j).*A32(i,j)+A22(i,j).*A31(i,j);
     2.*A11(i,j).*A31(i,j) 2.*A12(i,j).*A32(i,j) 2.*A13(i,j).*A33(i,j) A12(i,j).*A33(i,j)+A13(i,j).*A32(i,j) A11(i,j).*A33(i,j)+A13(i,j).*A31(i,j) A11(i,j).*A32(i,j)+A12(i,j).*A31(i,j);
     2.*A11(i,j).*A21(i,j) 2.*A12(i,j).*A22(i,j) 2.*A13(i,j).*A23(i,j) A12(i,j).*A23(i,j)+A13(i,j).*A22(i,j) A11(i,j).*A23(i,j)+A13(i,j).*A21(i,j) A11(i,j).*A22(i,j)+A12(i,j).*A21(i,j);];

Erot(:,:,i,j)=MA(:,:,i,j)'*e(:,:,i,j)*MA(:,:,i,j);
H_c(:,:,i,j)=inv(Cb(:,:,i,j)-F*Cb(:,:,i,j)*Erot(:,:,i,j)); % dry excess compliance of each single pore - pore impact approach
R1(i,j)=H_c(1,1,i,j)+H_c(1,2,i,j)+H_c(1,3,i,j); R2(i,j)=H_c(2,1,i,j)+H_c(2,2,i,j)+H_c(2,3,i,j); R3(i,j)=H_c(3,1,i,j)+H_c(3,2,i,j)+H_c(3,3,i,j); invKd(i,j)=R1(i,j)+R2(i,j)+R3(i,j); % inverse of bulk modulus of a dry pore
del(i,j)=(1./Kf(i,j)-1./K(i,j))./(invKd(i,j)); % parameter delta to compute saturated excess compliance below 
dH_c(:,:,i,j)=-(1./invKd(i,j)).*(1./(1+del(i,j))).*[R1(i,j).^2 R1(i,j).*R2(i,j) R1(i,j).*R3(i,j) 0 0 0; ...
                                                R1(i,j).*R2(i,j) R2(i,j).^2 R2(i,j).*R3(i,j) 0 0 0; ...
                                                R1(i,j).*R3(i,j) R2(i,j).*R3(i,j) R3(i,j).^2 0 0 0; ...
                                                0 0 0 0 0 0;0 0 0 0 0 0; 0 0 0 0 0 0;];  % saturated excess compliance for each single pore
end
end


%%%_choose_pore_sizes_by_UNCOMMENTING_ONE_OUT_OF_FOUR_opions_below_%%%

x13=rand(1,k).*ones(x,k); 
x14=rand(1,k).*ones(x/3,k); x15=rand(1,k).*ones(x/3,k); x16=rand(1,k).*ones(x/3,k);

for j=1:k  
for i=1:x
dm=3./(x);                                         % max pore density parameter --> a^3/V where 'a' is responsible for size


% UNCOMMENT ONE !!! %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   d(i,j)=dm.*rand;                                 % random sizes
%   d=dm*x13;                                        % identical sizes
%   d(i,j)=dm.*x13(i,j).*(1+0.2*rand-0.1);           % almost identical sizes
%   d=dm.*[x14; x15; x16];                           % pattern of ident sizes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi(i,j)=(1./3).*4.*pi.*y(i,j).*d(i,j);             % volume fraction of a single pore
Delta_I_c(:,:,i,j)=phi(i,j).*dH_c(:,:,i,j);         % fluid effect for single pore
pH_c(:,:,i,j)=phi(i,j).*H_c(:,:,i,j);
end
end

Delta_I=sum(Delta_I_c,[3 4])./k;                    % fluid effect - pore impact approach
phi_p=sum(phi,1);                                   % total volume fraction

for j=1:k   % set impact approach
H_p(:,:,j)=sum(pH_c(:,:,:,j),3)./phi_p(j);
R1=H_p(1,1,j)+H_p(1,2,j)+H_p(1,3,j); R2=H_p(2,1,j)+H_p(2,2,j)+H_p(2,3,j); R3=H_p(3,1,j)+H_p(3,2,j)+H_p(3,3,j); invKd_p=R1+R2+R3; % inverse of bulk modulus of a dry pore set
del=(1./Kf(1)-1./K(1))./(invKd_p);  % parameter delta to compute saturated excess compliance of a pore set below 
dH_p(:,:,j)=-(1./invKd_p).*(1./(1+del)).*[R1^2 R1*R2 R1*R3 0 0 0; R1*R2 R2^2 R2*R3 0 0 0; R1*R3 R2*R3 R3^2 0 0 0; 0 0 0 0 0 0;0 0 0 0 0 0; 0 0 0 0 0 0;];
Delta_II(:,:,j)=phi_p(j).*dH_p(:,:,j);
end

Delta_II=sum(Delta_II,3)./k;                            % fluid effect - set impact approach
Diff=Delta_I-Delta_II;                                  % absolute discrepancy
Diffro(jj)=100.*norm(Diff, 'fro')./norm(Delta_I,'fro'); % measurement of relative discrepancy
end

%% PLOT

plot(3*m,Diffro,'black');  
xlabel('Number of pores (x)'); ylabel('Discrepancy [%]');
ylim([0,100]); grid on;


%%%_outputs_for_Table1_%%%

mean(Diffro)
max(Diffro)
std(Diffro)
