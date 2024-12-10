clear all; close all; clc;

%  %  %  %  %  %  
%  DESCRIPTION %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  % 
%  This code compares two micromechanical approaches that predict       %
%  the effect of fluid (difference between undrained and drained        %
%  compliances) and its sensitivity to aspect ratio; we assume one set  %
%  composed of three connected subsets containing m=100 pores each.     %
%  Using this code, Figure 7 in Adamus et al. (2023), ''Multi-porous    %
%  extension of anisotropic poroelasticity: linkage with                %
%  micormechanics'', JGR, can be generated.                             %
%  Shape of each pore is assumed identical but the size and orientation %
%  MUST be chosen by the user (by twice uncommenting one out of four    %
%  options outlined below). Various measures of discrepancies are       %
%  provided. The discrepancy for each aspect ratio simulation can be    %
%  averaged over k iterations for curve smoothing purposes.             %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %


%%%_set_aspect ratios_%%%

n=logspace(-3.001,1,100);           % array of aspect ratios 
k=10;                              % number of iterations for each aspect ratio (curve smooting)
m=100;                              % pores in a subset
x=3*m;                              % total number of pores  

for jj=1:length(n)                  % loop calculations for each aspect ratio 'y' chosen ('n' array)
y=n(jj);


%%%_set_rock_background_%%%
          
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
     0, 0, 0, 0, 0, G(i,j); ]; 
end
end
F=[1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 2 0 0;
    0 0 0 0 2 0;   % correction matrix due to double-dot product (P:Cb or Cb:e)
    0 0 0 0 0 2;];


%%%_Eshelby matrix computation_%%%

if y==1
    y=1.001;
end
if y<1
g=(y./((1-y.^2).^(3./2))).*(acos(y)-y.*sqrt(1-y.^2));
else
g=(y./((y.^2-1).^(3./2))).*(y.*(y.^2-1).^(1./2)-acosh(y));
end

for j=1:k
for i=1:x

 e11(i,j)=(-3.*y.^2)./(8.*(1-v(i,j)).*(1-y.^2))+(g./(4.*(1-v(i,j)))).*(1-2.*v(i,j)+(9)./(4.*(1-y.^2)));
 e33(i,j)=(1./(1-v(i,j))).*(2-v(i,j)-1./(1-y.^2))+(g./(2.*(1-v(i,j)))).*(-4+2.*v(i,j)+3./(1-y.^2));
 e12(i,j)=(1./(8.*(1-v(i,j)))).*(1-1./(1-y.^2))+(g./(16.*(1-v(i,j)))).*(-4+8.*v(i,j)+3./(1-y.^2));
 e13(i,j)=(y.^2)./(2.*(1-v(i,j)).*(1-y.^2))-(g./(4.*(1-v(i,j)))).*(1-2.*v(i,j)+(3.*y.^2)./((1-y.^2)));
 e31(i,j)=(1./(2.*(1-v(i,j)))).*(-1+2.*v(i,j)+1./(1-y.^2))+(g./(4.*(1-v(i,j)))).*(2-4.*v(i,j)-3./(1-y.^2));
 e66(i,j)=(-1.*y.^2)./(8.*(1-v(i,j)).*(1-y.^2))+(g./(16.*(1-v(i,j)))).*(4-8.*v(i,j)+(3)./(1.*(1-y.^2)));
 e44(i,j)=(1./(4.*(1-v(i,j)))).*(1-2.*v(i,j)+(1+y.^2)./(1-y.^2))-(g./(8.*(1-v(i,j)))).*(1-2.*v(i,j)+3.*(1+y.^2)./(1-y.^2));
 e22(i,j)=e11(i,j); e21(i,j)=e12(i,j); e23(i,j)=e13(i,j); e32(i,j)=e31(i,j); e55(i,j)=e44(i,j);
 e(:,:,i,j)=[e11(i,j) e12(i,j) e13(i,j) 0 0 0; e21(i,j) e22(i,j) e23(i,j) 0 0 0; e31(i,j) e32(i,j) e33(i,j) 0 0 0;
             0 0 0 e44(i,j) 0 0;    0 0 0 0 e55(i,j) 0;     0 0 0 0 0 e66(i,j); ];
end
end


%%%_choose_pore_orientations_by_UNCOMMENTING_ONE_OUT_OF_FOUR_opions_in_a_box_below_%%%

x4=rand(3,1,k).*ones(3,x,k); x4=x4./vecnorm(x4); x5=360*rand(1,k).*ones(x,k);  % helpers to be used below 
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

x13=rand(1,k).*ones(x,k); % helpers to be used below 
x14=rand(1,k).*ones(x/3,k); x15=rand(1,k).*ones(x/3,k); x16=rand(1,k).*ones(x/3,k);

for j=1:k  
for i=1:x
dm=1./(m);                                         % max pore density parameter --> a^3/V where 'a' is responsible for size

% UNCOMMENT ONE !!! %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   d(i,j)=dm.*rand;                                 % random sizes
%   d=dm*x13;                                        % identical sizes
%   d(i,j)=dm.*x13(i,j).*(1+0.2*rand-0.1);           % almost identical sizes
%   d=dm.*[x14; x15; x16];                           % pattern of ident sizes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi(i,j)=(1./3).*4.*pi.*y.*d(i,j);                   % volume fraction of a single pore
Delta_I_c(:,:,i,j)=phi(i,j).*dH_c(:,:,i,j);          % fluid effect for single pore
pH_c(:,:,i,j)=phi(i,j).*H_c(:,:,i,j);
end
end

Delta_I=sum(Delta_I_c,[3 4])./k;                   % fluid effect - pore impact approach
phi_p=sum(phi,1);                                    % total volume fraction

for j=1:k % set impact approach
H_p(:,:,j)=sum(pH_c(:,:,:,j),3)./phi_p(j);
R1=H_p(1,1,j)+H_p(1,2,j)+H_p(1,3,j); R2=H_p(2,1,j)+H_p(2,2,j)+H_p(2,3,j); R3=H_p(3,1,j)+H_p(3,2,j)+H_p(3,3,j); invKd_p=R1+R2+R3; % inverse of bulk modulus of a dry pore set
del=(1./Kf(1)-1./K(1))./(invKd_p); 
dH_p(:,:,j)=-(1./invKd_p).*(1./(1+del)).*[R1^2 R1*R2 R1*R3 0 0 0; R1*R2 R2^2 R2*R3 0 0 0; R1*R3 R2*R3 R3^2 0 0 0; 0 0 0 0 0 0;0 0 0 0 0 0; 0 0 0 0 0 0;];
Delta_II(:,:,j)=phi_p(j).*dH_p(:,:,j);
end

Delta_II=sum(Delta_II,3)./k;                      % fluid effect - set impact approach
Diff=Delta_I-Delta_II;                            % absolute discrepancy between two fluid effects
Delta=sum(pH_c,[3 4])./k;

%%%_outputs_various_measures_of_discrepancy_%%%

Diffro(jj)=100.*norm(Diff, 'fro')./norm(Delta_I,'fro');        
Diffro2(jj)=100.*norm(Diff, 'fro')./norm(Delta,'fro');
Diffro2a(jj)=100.*norm(Delta_I, 'fro')./norm(Delta,'fro');
Diffro2b(jj)=100.*norm(Delta_II, 'fro')./norm(Delta,'fro');

end


%% PLOT

loglog(n,Diffro,'black'); hold on 
loglog(n,Diffro2,'blue'); hold on 
loglog(n,Diffro2a,'red'); hold on 
loglog(n,Diffro2b,'red--'); hold on 
grid on
xlabel('Aspect ratio'); ylabel('Average discrepancy [%]');
ylim([0.01,100])



