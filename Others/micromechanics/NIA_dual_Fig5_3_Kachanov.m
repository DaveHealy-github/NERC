clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Choose background and inclusion properties              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vol=30*30*30;
% vol=62.5*pi*12.5^2;         % volume of the sample (in mm3)
E=1;                        % background Young modulus (in GPa)
v=0.23;                     % background Poisson ratio
l=E.*v./((1+v).*(1-2.*v));  % background lame parameter
G=E./(2+2.*v);              % background rigidity
K=E./(3.*(1-2.*v));         % background bulk modulus
Kf=0/0.45;                  % min 0GPa (air) max 8GPa (ice at zero deg) 

% first inclusion type
for j=1:2000
n(1)=1*j;                     % number of defined inclusions
a(1)=2;                     % radius of circular inclusion (in mm)
y(1)=0.125;                    % aspect ratio of inclusion
or(1)=0;                    % orientation of inclusion 0=horizontal 1=vertical

% second inclusion type
n(2)=0;                     % number of defined inclusions
a(2)=8;                     % radius of circular inclusion (in mm)
y(2)=0.2;                    % aspect ratio of inclusion
or(2)=1;                    % orientation of inclusion 0=horizontal 1=vertical

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Compute effective elastic properties                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 Cb= [ l+2.*G, l, l, 0, 0, 0;
     l,  l+2.*G, l, 0, 0, 0;
     l, l,  l+2.*G, 0, 0, 0;
     0, 0, 0, G, 0, 0;
     0, 0, 0, 0, G, 0; 
     0, 0, 0, 0, 0, G; ]; 

 C1= [Kf/3, Kf/3, Kf/3, 0, 0, 0;
     Kf/3,  Kf/3, Kf/3, 0, 0, 0;
     Kf/3, Kf/3,  Kf/3, 0, 0, 0;
     0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0; 
     0, 0, 0, 0, 0, 0; ]; 
Sb=inv(Cb);

I=[1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 1 0 0;
    0 0 0 0 1 0;   
    0 0 0 0 0 1;];

F=[1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 2 0 0;
    0 0 0 0 2 0;   % correction matrix due to double-dot product (P:Cb or Cb:e)
    0 0 0 0 0 2;];

B11=1; B12=0; B13=0; B21=0; B22=0; B23=1; B31=0; B32=-1; B33=0; % cracks become having x2 normal
MB= [B11*B11 B12*B12 B13*B13 B12*B13 B11*B13 B11*B12;
     B21*B21 B22*B22 B23*B23 B22*B23 B21*B23 B21*B22;
     B31*B31 B32*B32 B33*B33 B32*B33 B31*B33 B31*B32;
     2*B21*B31 2*B22*B32 2*B23*B33 B22*B33+B23*B32 B21*B33+B23*B31 B21*B32+B22*B31;
     2*B11*B31 2*B12*B32 2*B13*B33 B12*B33+B13*B32 B11*B33+B13*B31 B11*B32+B12*B31;
     2*B11*B21 2*B12*B22 2*B13*B23 B12*B23+B13*B22 B11*B23+B13*B21 B11*B22+B12*B21;];

H1=[0.5 0.5 0 0 0 0;
    0.5 0.5 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;];

H5=[0.5 -0.5 0 0 0 0;
    -0.5 0.5 0 0 0 0;
         0 0 0 0 0 0;
         0 0 0 0 0 0;
         0 0 0 0 0 0;
        0 0 0 0 0 2;]; % not sure if not 0.5 for ij=66

H2=[0 0 1 0 0 0;
    0 0 1 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;];

H3=[0 0 0 0 0 0;
    0 0 0 0 0 0;
    1 1 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;];

H6=[0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 2 0 0;
    0 0 0 0 2 0;
    0 0 0 0 0 0;]; % not sure if not 0.5 for both

H4=[0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;];

%%% LOOP - y %%%
for i=1:2

if y(i)>1   
 s(i)=1./(1-y(i).^2)-y(i)./(1-y(i).^2).*(1/(sqrt(y(i).^2-1))).*acosh(y(i));
 else
 s(i)=1./(1-y(i).^2)-y(i)./(1-y(i).^2).*(1/(sqrt(1-y(i).^2))).*acos(y(i));
end

t=0.5*(1-s(i));
psi1=(y(i).^2.*(4*t-1)-t)./(4*(1-y(i).^2));
psi2=(y(i).^2.*(1-2*t)-t)./(4*(1-y(i).^2));
psi4=(3.*t-1)./(2*(1-y(i).^2));

  p1=(1./G).*(t+psi1./(1-v));
  p2=psi2./(G.*(1-v));
  p3=p2;
  p4=(1./G).*(1-2*t+psi4./(1-v));
  p5=(1./G).*(t+psi1./(2*(1-v)));
  p6=(1./G).*(0.5*(1-t)+2*psi2./(1-v));

 P=p1*H1+p2*H2+p3*H3+p4*H4+p5*H5+p6*H6; 
 Prot=MB'*P*MB;

if or(i)==0
 N(:,:,i)=((C1-Cb)^(-1)+P)^(-1);
  %N(:,:,i)=((-Cb)^(-1)+P)^(-1); %if dry, to be used in Shapiro analogous
else
 N(:,:,i)=((C1-Cb)^(-1)+Prot)^(-1);
   %N(:,:,i)=((-Cb)^(-1)+Prot)^(-1); %if dry, to be used in Shapiro analogous
end

%  R1=N(1,1,i)+N(1,2,i)+N(1,3,i); R2=N(2,1,i)+N(2,2,i)+N(2,3,i); R3=N(3,1,i)+N(3,2,i)+N(3,3,i); R=R1+R2+R3; del=(Kf-K)./(R); % effect of fluid
%  dN(:,:,i)=-(1./R).*(1./(1+del)).*[R1.^2 R1.*R2 R1.*R3 0 0 0; R1.*R2 R2.^2 R2.*R3 0 0 0; R1.*R3 R2.*R3 R3.^2 0 0 0; 0 0 0 0 0 0;0 0 0 0 0 0; 0 0 0 0 0 0;];
% sigma3=10; Q33(i)=-1/(1+del)*(R3/R); % needed for pore pressures in uniaxial test
% analogous to Shapiro and Kachanov so that set-impact and pore-impact
% approaches can be obtained by manipulating dN (or dH by analogy)...
% C_eff2=Cb+n(1)*p(1)*(N(:,:,1)+dN(:,:,1))+n(2)*p(2)*(N(:,:,2)+dN(:,:,2));

 Q=Cb-F*Cb*(F*P*Cb);

if or(i)==0
H(:,:,i)=inv(Q);
else
H(:,:,i)=inv(Cb-F*Cb*(Prot)*F*Cb);
end

R1=H(1,1,i)+H(1,2,i)+H(1,3,i); R2=H(2,1,i)+H(2,2,i)+H(2,3,i); R3=H(3,1,i)+H(3,2,i)+H(3,3,i); R=R1+R2+R3; del=(1./Kf-1./K)./(R); % effect of fluid
dH(:,:,i)=-(1./R).*(1./(1+del)).*[R1.^2 R1.*R2 R1.*R3 0 0 0; R1.*R2 R2.^2 R2.*R3 0 0 0; R1.*R3 R2.*R3 R3.^2 0 0 0; 0 0 0 0 0 0;0 0 0 0 0 0; 0 0 0 0 0 0;];
sigma3=10; Q33(i)=-1/(1+del)*(R3/R); % needed for pore pressures in uniaxial test

p(i)=(4/3).*(pi.*y(i).*a(i).^3)./vol;  % volume fraction of a defined pore
end

p_total(j)=n(1)*p(1)+n(2)*p(2);
e_total(j)=(n(1)*a(1).^3./vol)+(n(2)*a(2).^3./vol);
n_total(j)=n(1)+n(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  RESULTS  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
%%%% NIA compliance %%%%
%%%%%%%%%%%%%%%%%%%%%%%%
S_eff1=Sb+n(1)*p(1)*(H(:,:,1)+dH(:,:,1))+n(2)*p(2)*(H(:,:,2)+dH(:,:,2));
C_eff1=inv(S_eff1);
E_eff1=1/S_eff1(3,3);    % uniaxial test - assume x3
K_eff1=(C_eff1(1,1)+C_eff1(1,2)+C_eff1(1,3)+C_eff1(2,1)+C_eff1(2,2)+C_eff1(2,3)+C_eff1(3,1)+C_eff1(3,2)+C_eff1(3,3))*(1/9);
K_ratio_NIA1(j)=K_eff1/K;
E_ratio_NIA1(j)=E_eff1/E;

%%%%%%%%%%%%%%%%%%%%%%%%
%%%% NIA stiffness  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%

C_eff2=Cb+n(1)*p(1)*(N(:,:,1))+n(2)*p(2)*(N(:,:,2));
S_eff2=inv(C_eff2);
E_eff2=1/S_eff2(3,3);    % uniaxial test - assume x3
K_eff2=(C_eff2(1,1)+C_eff2(1,2)+C_eff2(1,3)+C_eff2(2,1)+C_eff2(2,2)+C_eff2(2,3)+C_eff2(3,1)+C_eff2(3,2)+C_eff2(3,3))*(1/9);
K_ratio_NIA2(j)=K_eff2/K;
E_ratio_NIA2(j)=E_eff2/E;
end


%%
figure; 
plot(p_total,K_ratio_NIA1,'green'); hold on
plot(p_total,K_ratio_NIA2,'blue')
legend('NIA compl','NIA stiff')
title('Effect of spheres viewed by dual NIA')
yticks([0 0.2 0.4 0.6 0.8 1])
% xticks([0 0.2 0.4 0.6 0.8 1])
ylabel('K/K0')
xlabel('Porosity')
ylim([0 1])
% xlim([0 1])
grid on

figure;
plot(n_total,K_ratio_NIA1,'green'); hold on
plot(n_total,K_ratio_NIA2,'blue')
legend('NIA compl','NIA stiff')
title('Effect of spheres viewed by dual NIA')
yticks([0 0.2 0.4 0.6 0.8 1])
% xticks([0 0.2 0.4 0.6 0.8 1])
ylabel('K/K0')
xlabel('Number of cracks')
ylim([0 1])
% xlim([0 1])
grid on

figure;
plot(e_total,K_ratio_NIA1,'green'); hold on
plot(e_total,K_ratio_NIA2,'blue')
legend('NIA compl','NIA stiff')
title('Effect of spheres viewed by dual NIA')
yticks([0 0.2 0.4 0.6 0.8 1])
% xticks([0 0.2 0.4 0.6 0.8 1])
ylabel('K/K0')
xlabel('Crack density')
ylim([0 1])
% xlim([0 1])
grid on

figure;
plot(p_total,E_ratio_NIA1,'green'); hold on
plot(p_total,E_ratio_NIA2,'blue')
legend('NIA compl','NIA stiff')
title('Effect of spheres viewed by dual NIA')
yticks([0 0.2 0.4 0.6 0.8 1])
% xticks([0 0.2 0.4 0.6 0.8 1])
ylabel('E/E0')
xlabel('Porosity')
ylim([0 1])
% xlim([0 1])
grid on

figure;
plot(n_total,E_ratio_NIA1,'green'); hold on
plot(n_total,E_ratio_NIA2,'blue')
legend('NIA compl','NIA stiff')
title('Effect of spheres viewed by dual NIA')
yticks([0 0.2 0.4 0.6 0.8 1])
% xticks([0 0.2 0.4 0.6 0.8 1])
ylabel('E/E0')
xlabel('Number of cracks')
ylim([0 1])
% xlim([0 1])
grid on

figure;
plot(e_total,E_ratio_NIA1,'green'); hold on
plot(e_total,E_ratio_NIA2,'blue')
legend('NIA compl','NIA stiff')
title('Effect of spheres viewed by dual NIA')
yticks([0 0.2 0.4 0.6 0.8 1])
% xticks([0 0.2 0.4 0.6 0.8 1])
ylabel('E/E0')
xlabel('Crack density')
ylim([0 1])
% xlim([0 1])
grid on

E_ratio_NIA1(300)
E_ratio_NIA2(300)

























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Compute conjectured effective elasticity (poroelastic approach)   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% H=(n(1)*p(1)*H(:,:,1)+n(2)*p(2)*H(:,:,2))/(n(1)*p(1)+n(2)*p(2));  % connected pores treated as one inhomogeneity
% R1=H(1,1)+H(1,2)+H(1,3); R2=H(2,1)+H(2,2)+H(2,3); R3=H(3,1)+H(3,2)+H(3,3); R=R1+R2+R3; del=(1./Kf-1./K)./(R); % effect of fluid
% dH=-(1./R).*(1./(1+del)).*[R1.^2 R1.*R2 R1.*R3 0 0 0; R1.*R2 R2.^2 R2.*R3 0 0 0; R1.*R3 R2.*R3 R3.^2 0 0 0; 0 0 0 0 0 0;0 0 0 0 0 0; 0 0 0 0 0 0;];
% Q33=-1/(1+del)*(R3/R);  
% 
% S_eff_conjecture=Sb+(n(1)*p(1)+n(2)*p(2))*(H+dH);
% C_eff_conjecture=inv(S_eff_conjecture);
% E_eff_conjecture=1/S_eff_conjecture(3,3)
% 
% E_eff-E_eff_conjecture

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Comparison of pore pressures  (both approaches)           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pf1=Q33(1)*sigma3;   % pore pressure in pore 1
% pf2=Q33(2)*sigma3;   % pore pressure in pore 2
% pf=Q33*sigma3;       % uniformed pore pressure