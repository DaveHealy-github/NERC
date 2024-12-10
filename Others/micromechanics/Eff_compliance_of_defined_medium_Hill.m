clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Choose background and inclusion properties              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vol=62.5*pi*12.5^2;         % volume of the sample (in mm3)
E=2.25;                        % background Young modulus (in GPa)
v=0.2;                     % background Poisson ratio
l=E.*v./((1+v).*(1-2.*v));  % background lame parameter
G=E./(2+2.*v);              % background rigidity
K=E./(3.*(1-2.*v));         % background bulk modulus
Kf=0/0.45;                  % min 0GPa (air) max 8GPa (ice at zero deg) 

% first inclusion type
n(1)=1;                     % number of defined inclusions
a(1)=4;                     % radius of circular inclusion (in mm)
y(1)=1.00001;                    % aspect ratio of inclusion
or(1)=0;                    % orientation of inclusion 0=horizontal 1=vertical

% second inclusion type
n(2)=0;                     % number of defined inclusions
a(2)=8;                     % radius of circular inclusion (in mm)
y(2)=0.2;                    % aspect ratio of inclusion
or(2)=1;                     % orientation of inclusion 0=horizontal 1=vertical

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Compute effective elastic properties                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 Cb= [ l+2.*G, l, l, 0, 0, 0;
     l,  l+2.*G, l, 0, 0, 0;
     l, l,  l+2.*G, 0, 0, 0;
     0, 0, 0, G, 0, 0;
     0, 0, 0, 0, G, 0; 
     0, 0, 0, 0, 0, G; ]; 
Sb=inv(Cb);

C1= [Kf/3, Kf/3, Kf/3, 0, 0, 0;
     Kf/3,  Kf/3, Kf/3, 0, 0, 0;
     Kf/3, Kf/3,  Kf/3, 0, 0, 0;
     0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0; 
     0, 0, 0, 0, 0, 0; ]; 

S1= [3/Kf, 3/Kf, 3/Kf, 0, 0, 0;
     3/Kf,  3/Kf, 3/Kf, 0, 0, 0;
     3/Kf, 3/Kf,  3/Kf, 0, 0, 0;
     0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0; 
     0, 0, 0, 0, 0, 0; ]; 

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

%%%%%% above something wrong %%%%%%%

  p1=(1./G).*(t+psi1./(1-v));
  p2=psi2./(G.*(1-v));
  p3=p2;
  p4=(1./G).*(1-2*t+psi4./(1-v));
  p5=(1./G).*(t+psi1./(2*(1-v)));
  p6=(1./G).*(0.5*(1-t)+2*psi2./(1-v));

 P=p1*H1+p2*H2+p3*H3+p4*H4+p5*H5+p6*H6; 
 Prot=MB'*P*MB;

%  Omega=(I+F*P*(C1-Cb))^(-1);
%  Omegarot=(I+F*Prot*(C1-Cb))^(-1);
%  Y=F*(C1*Omega)*F*Sb;
%  Yrot=F*(C1*Omegarot)*F*Sb;
% 
% if or(i)==0
% H(:,:,i)=(S1-Sb)*F*Y;
% else
% H(:,:,i)=(S1-Sb)*F*Yrot;
% end

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

p_total=n(1)*p(1)+n(2)*p(2)

S_eff=Sb+n(1)*p(1)*(H(:,:,1)+dH(:,:,1))+n(2)*p(2)*(H(:,:,2)+dH(:,:,2));
C_eff=inv(S_eff);
E_eff=1/S_eff(3,3)    % uniaxial test - assume x3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Compute conjectured effective elasticity (poroelastic approach)   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H=(n(1)*p(1)*H(:,:,1)+n(2)*p(2)*H(:,:,2))/(n(1)*p(1)+n(2)*p(2));  % connected pores treated as one inhomogeneity
R1=H(1,1)+H(1,2)+H(1,3); R2=H(2,1)+H(2,2)+H(2,3); R3=H(3,1)+H(3,2)+H(3,3); R=R1+R2+R3; del=(1./Kf-1./K)./(R); % effect of fluid
dH=-(1./R).*(1./(1+del)).*[R1.^2 R1.*R2 R1.*R3 0 0 0; R1.*R2 R2.^2 R2.*R3 0 0 0; R1.*R3 R2.*R3 R3.^2 0 0 0; 0 0 0 0 0 0;0 0 0 0 0 0; 0 0 0 0 0 0;];
Q33n=-1/(1+del)*(R3/R);  

S_eff_conjecture=Sb+(n(1)*p(1)+n(2)*p(2))*(H+dH);
C_eff_conjecture=inv(S_eff_conjecture);
E_eff_conjecture=1/S_eff_conjecture(3,3)

E_eff-E_eff_conjecture;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Comparison of pore pressures  (both approaches)           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pf1=Q33(1)*sigma3;   % pore pressure in pore 1
pf2=Q33(2)*sigma3;   % pore pressure in pore 2
pf=Q33n*sigma3;       % uniformed pore pressure
