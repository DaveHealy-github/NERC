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

% first inclusion ty55pe
n(1)=1;                     % number of defined inclusions
a(1)=10;                     % radius of circular inclusion (in mm)
y(1)=0.1;                    % aspect ratio of inclusion
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
Sb=inv(Cb);

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

%%% LOOP - y %%%
for i=1:2

 if y(i)>1   
 g(i)=(y(i)./((y(i).^2-1).^(3./2))).*(y(i).*(y(i).^2-1).^(1./2)-acosh(y(i)));
 else
 g(i)=(y(i)./((1-y(i).^2).^(3./2))).*(acos(y(i))-y(i).*sqrt(1-y(i).^2));
 end

 e11=(-3.*y(i).^2)./(8.*(1-v).*(1-y(i).^2))+(g(i)./(4.*(1-v))).*(1-2.*v+(9)/(4.*(1-y(i).^2)));
 e22=e11;
 e33=(1./(1-v)).*(2-v-1./(1-y(i).^2))+(g(i)./(2.*(1-v))).*(-4+2.*v+3./(1-y(i).^2));
 e12=(1./(8.*(1-v))).*(1-1./(1-y(i).^2))+(g(i)./(16.*(1-v))).*(-4+8.*v+3./(1-y(i).^2));
 e21=e12;
 e13=(y(i).^2)./(2.*(1-v).*(1-y(i).^2))-(g(i)./(4.*(1-v))).*(1-2.*v+(3.*y(i).^2)/((1-y(i).^2)));
 e23=e13;
 e31=(1./(2.*(1-v))).*(-1+2.*v+1./(1-y(i).^2))+(g(i)./(4.*(1-v))).*(2-4.*v-3./(1-y(i).^2));
 e32=e31;
 e66=(-1.*y(i).^2)./(8.*(1-v).*(1-y(i).^2))+(g(i)./(16.*(1-v))).*(4-8.*v+(3)/(1.*(1-y(i).^2)));
 e44=(1./(4.*(1-v))).*(1-2.*v+(1+y(i).^2)./(1-y(i).^2))-(g(i)./(8.*(1-v)))*(1-2.*v+3.*(1+y(i).^2)./(1-y(i).^2));
 e55=e44;

 e=[e11 e12 e13 0 0 0;
    e21 e22 e23 0 0 0;
    e31 e32 e33 0 0 0;
        0 0 0 e44 0 0;    
        0 0 0 0 e55 0;    
        0 0 0 0 0 e66; ]; 
 Erot=MB'*e*MB;

 Q=Cb-F*Cb*e;

if or(i)==0
H(:,:,i)=inv(Q);
else
H(:,:,i)=inv(Cb-F*Cb*(Erot));
end

R1=H(1,1,i)+H(1,2,i)+H(1,3,i); R2=H(2,1,i)+H(2,2,i)+H(2,3,i); R3=H(3,1,i)+H(3,2,i)+H(3,3,i); R=R1+R2+R3; del=(1./Kf-1./K)./(R); % effect of fluid
dH(:,:,i)=-(1./R).*(1./(1+del)).*[R1.^2 R1.*R2 R1.*R3 0 0 0; R1.*R2 R2.^2 R2.*R3 0 0 0; R1.*R3 R2.*R3 R3.^2 0 0 0; 0 0 0 0 0 0;0 0 0 0 0 0; 0 0 0 0 0 0;];
sigma3=10; Q33(i)=-1/(1+del)*(R3/R); % needed for pore pressures in uniaxial test

p(i)=(4/3).*(pi.*y(i).*a(i).^3)./vol;  % volume fraction of a defined pore
end

S_eff=Sb+n(1)*p(1)*(H(:,:,1)+dH(:,:,1))+n(2)*p(2)*(H(:,:,2)+dH(:,:,2));
C_eff=inv(S_eff);
E_eff=1/S_eff(3,3)    % uniaxial test - assume x3

pf1=Q33(1)*sigma3;   % pore pressure in pore 1
pf2=Q33(2)*sigma3;   % pore pressure in pore 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Compute conjectured effective elasticity (poroelastic approach)   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H=(n(1)*p(1)*H(:,:,1)+n(2)*p(2)*H(:,:,2))/(n(1)*p(1)+n(2)*p(2));  % connected pores treated as one inhomogeneity
R1=H(1,1)+H(1,2)+H(1,3); R2=H(2,1)+H(2,2)+H(2,3); R3=H(3,1)+H(3,2)+H(3,3); R=R1+R2+R3; del=(1./Kf-1./K)./(R); % effect of fluid
dH=-(1./R).*(1./(1+del)).*[R1.^2 R1.*R2 R1.*R3 0 0 0; R1.*R2 R2.^2 R2.*R3 0 0 0; R1.*R3 R2.*R3 R3.^2 0 0 0; 0 0 0 0 0 0;0 0 0 0 0 0; 0 0 0 0 0 0;];
Q33=-1/(1+del)*(R3/R); pf=Q33*sigma3; 

S_eff_conjecture=Sb+(n(1)*p(1)+n(2)*p(2))*(H+dH);
C_eff_conjecture=inv(S_eff_conjecture);
E_eff_conjecture=1/S_eff_conjecture(3,3)

E_eff-E_eff_conjecture