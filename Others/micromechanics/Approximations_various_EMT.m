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
Kf=0.0001/0.45; Gf=0.0001;     % min 0GPa (air) max 8GPa (ice at zero deg) 

% first inclusion type
for j=1:100
n(1)=1*j;                     % number of defined inclusions
a(1)=3;                     % radius of circular inclusion (in mm)
y(1)=1;                    % aspect ratio of inclusion
or(1)=0;                    % orientation of inclusion 0=horizontal 1=vertical

% second inclusion type
n(2)=0;                     % number of defined inclusions
a(2)=2;                     % radius of circular inclusion (in mm)
y(2)=2;                    % aspect ratio of inclusion
or(2)=1;                    % orientation of inclusion 0=horizontal 1=vertical

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Compute effective elastic properties                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 Cb= [ l+2.*G, l, l, 0, 0, 0;
     l,  l+2.*G, l, 0, 0, 0;
     l, l,  l+2.*G, 0, 0, 0;
     0, 0, 0, G, 0, 0;
     0, 0, 0, 0, G, 0; 
     0, 0, 0, 0, 0, G; ]; Sb=inv(Cb); 

 I1= [(1/3), (1/3), (1/3), 0, 0, 0;
       (1/3), (1/3), (1/3), 0, 0, 0;
       (1/3), (1/3), (1/3), 0, 0, 0;
            0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0; 
            0, 0, 0, 0, 0, 0; ]; 

 I2= [(2/3), -(1/3), -(1/3), 0, 0, 0;
       -(1/3), (2/3), -(1/3), 0, 0, 0;
       -(1/3), -(1/3), (2/3), 0, 0, 0;
            0, 0, 0, 2, 0, 0;
            0, 0, 0, 0, 2, 0; 
            0, 0, 0, 0, 0, 2; ]; 

 C1=3*Kf*I1+Gf*I2;          % stiffness of inhomogeneity (see Kach book)
 S1=1/(3*Kf)*I1+(1/Gf)*I2;  % compliance of inhomogeneity (see Kach book)

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
        0 0 0 0 0 2;]; 

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
    0 0 0 0 0 0;]; 

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
elseif y(i)<1
 s(i)=1./(1-y(i).^2)-y(i)./(1-y(i).^2).*(1/(sqrt(1-y(i).^2))).*acos(y(i));
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
else
  s(i)=1/3;
  p1=(1-2*v)/(6*G*(1-v));
  p2=(4-5*v)/(15*G*(1-v));
  P=p1*I1+p2*I2;
end

 Prot=MB'*P*MB;
 Omega(:,:,i)=(I+F*P*(C1-Cb))^(-1);

if or(i)==0
 N(:,:,i)=((C1-Cb)^(-1)+P)^(-1);
 else
 N(:,:,i)=((C1-Cb)^(-1)+Prot)^(-1);
 end

 Q=Cb-F*Cb*(F*P*Cb);
 Gamma(:,:,i)=(I+F*Q*(S1-Sb))^(-1);

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

p_total(j)=n(1).*p(1)+n(2).*p(2);
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
E_eff1(j)=1/S_eff1(3,3);    % uniaxial test - assume x3
K_eff1=(C_eff1(1,1)+C_eff1(1,2)+C_eff1(1,3)+C_eff1(2,1)+C_eff1(2,2)+C_eff1(2,3)+C_eff1(3,1)+C_eff1(3,2)+C_eff1(3,3))*(1/9);
K_ratio_NIA1(j)=K_eff1/K;
E_ratio_NIA1(j)=E_eff1(j)/E;
%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Reuss average  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%
S_eff_R=Sb.*(1-p_total(j))+S1.*p_total(j);
C_eff_R=inv(S_eff_R);
E_eff_R=1./(S_eff_R(3,3));
K_eff_R=(C_eff_R(1,1)+C_eff_R(1,2)+C_eff_R(1,3)+C_eff_R(2,1)+C_eff_R(2,2)+C_eff_R(2,3)+C_eff_R(3,1)+C_eff_R(3,2)+C_eff_R(3,3))*(1/9);
% K_eff_R=1./(((1-p_total(j))./K)+p_total(j)./Kf);
% G_eff_R=1./(((1-p_total(j))./G)+p_total(j)./Gf);
% E_eff_R=1./(1./(9*K_eff_R)+1./(3*G_eff_R));       % alternative description
K_ratio_R(j)=K_eff_R./K;
E_ratio_R(j)=E_eff_R./E;

%%%%%%%%%%%%%%%%%%%%%%%%
%%%% NIA stiffness  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%
C_eff2=Cb+n(1)*p(1)*(N(:,:,1))+n(2)*p(2)*(N(:,:,2));
S_eff2=inv(C_eff2);
E_eff2(j)=1/S_eff2(3,3);    % uniaxial test - assume x3
K_eff2=(C_eff2(1,1)+C_eff2(1,2)+C_eff2(1,3)+C_eff2(2,1)+C_eff2(2,2)+C_eff2(2,3)+C_eff2(3,1)+C_eff2(3,2)+C_eff2(3,3))*(1/9);
K_ratio_NIA2(j)=K_eff2/K;
E_ratio_NIA2(j)=E_eff2(j)/E;
%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Voigt average  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%
C_eff_V=Cb.*(1-p_total(j))+C1.*p_total(j);
S_eff_V=inv(C_eff_V);
E_eff_V(j)=1./S_eff_V(3,3);    % uniaxial test - assume x3
K_eff_V=(C_eff_V(1,1)+C_eff_V(1,2)+C_eff_V(1,3)+C_eff_V(2,1)+C_eff_V(2,2)+C_eff_V(2,3)+C_eff_V(3,1)+C_eff_V(3,2)+C_eff_V(3,3))*(1/9);
K_ratio_V(j)=K_eff_V./K;
E_ratio_V(j)=E_eff_V(j)./E;
%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Hashkin bound  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%
K_HSmin=Kf+((1-p_total(j)).*(K-Kf))/(1+p_total(j).*(K-Kf)/(Kf+4*Gf/3));
K_HSmax=Kf+((1-p_total(j)).*(K-Kf))/(1+p_total(j).*(K-Kf)/(Kf+4*G/3));
K_ratio_HSmin(j)=K_HSmin/K;
K_ratio_HSmax(j)=K_HSmax/K;

G_HSmin=Gf+((1-p_total(j)).*(G-Gf))/(1+p_total(j).*(G-Gf)/(Gf+(3/2)*(1/Gf+10/(9*Kf+8*Gf))^(-1)));
G_HSmax=Gf+((1-p_total(j)).*(G-Gf))/(1+p_total(j).*(G-Gf)/(Gf+(3/2)*(1/G+10/(9*K+8*G))^(-1)));
G_ratio_HSmin(j)=G_HSmin/G;
G_ratio_HSmax(j)=G_HSmax/G;

E_HSmin=1./(1./(9*K_HSmin)+1./(3*G_HSmin));
E_HSmax=1./(1./(9*K_HSmax)+1./(3*G_HSmax));
E_ratio_HSmin(j)=E_HSmin/E;
E_ratio_HSmax(j)=E_HSmax/E;
%%%%%%%%%%%%%%%%%%%%%%%%
%%% Self-Cons scheme %%%
%%%%%%%%%%%%%%%%%%%%%%%%
for jj=1:1000000
v_eff(jj)=rand-0.5;
G_eff(jj)=rand*G;

K_eff(jj)=(2.*G_eff(jj).*(1+v_eff(jj)))./(3.*(1-2.*v_eff(jj)));
alpha=3*K_eff(jj)/(3*K_eff(jj)+4*G_eff(jj));
beta=(6/5)*(K_eff(jj)+2.*G_eff(jj))./(3.*K_eff(jj)+2.*G_eff(jj));

K_SC=K_eff(jj)-p_total(j).*K_eff(jj).*(Kf-K)./(K_eff(jj)+alpha.*(Kf-K_eff(jj)));
G_SC=G_eff(jj)-p_total(j).*G_eff(jj).*(Gf-G)./(G_eff(jj)+beta.*(Gf-G_eff(jj)));

err1(jj)=(K_SC-K)^2;
err2(jj)=(G_SC-G)^2;
err(jj)=sqrt(err1(jj)+err2(jj));
end
[row1] = find (err==min(err));

K_true=K_eff(row1);
G_true=G_eff(row1);
E_true=9*K_true*G_true/(3*K_true+G_true);

K_ratio_SC(j)=K_true./K;
E_ratio_SC(j)=E_true./E;
G_ratio_SC(j)=G_true./G;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Differential scheme %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Mori-Tanaka scheme  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_eff_MT=Cb+p_total(j).*((C1-Cb)*Omega(:,:,1))*((1-p_total(j)).*F+p_total(j).*Omega(:,:,1))^(-1);
S_eff_MT=inv(C_eff_MT);
E_eff_MT=1./S_eff_MT(3,3);    % uniaxial test - assume x3
K_eff_MT=(C_eff_MT(1,1)+C_eff_MT(1,2)+C_eff_MT(1,3)+C_eff_MT(2,1)+C_eff_MT(2,2)+C_eff_MT(2,3)+C_eff_MT(3,1)+C_eff_MT(3,2)+C_eff_MT(3,3))*(1/9);
K_ratio_MT(j)=K_eff_MT./K;
E_ratio_MT(j)=E_eff_MT./E;

% S_eff_MT2=Sb+p_total(j).*((S1-Sb)*Gamma(:,:,1))*((1-p_total(j)).*F+p_total(j).*Gamma(:,:,1))^(-1);
% C_eff_MT2=inv(S_eff_MT2);
% E_eff_MT2=1./S_eff_MT2(3,3);    % uniaxial test - assume x3
% K_eff_MT2=(C_eff_MT2(1,1)+C_eff_MT2(1,2)+C_eff_MT2(1,3)+C_eff_MT2(2,1)+C_eff_MT2(2,2)+C_eff_MT2(2,3)+C_eff_MT2(3,1)+C_eff_MT2(3,2)+C_eff_MT2(3,3))*(1/9);
% K_ratio_MT2(j)=K_eff_MT2./K;
% E_ratio_MT2(j)=E_eff_MT2./E;
end


%%
figure; 
plot(p_total,K_ratio_NIA1,'green'); hold on
plot(p_total,K_ratio_NIA2,'blue'); hold on
plot(p_total,K_ratio_R,'black'); hold on
plot(p_total,K_ratio_V,'black--'); hold on
plot(p_total,K_ratio_HSmin,'red--'); hold on
plot(p_total,K_ratio_HSmax,'red'); hold on
plot(p_total,K_ratio_SC,'cyan'); hold on
plot(p_total,K_ratio_MT,'yellow--'); hold on
legend('NIA compl','NIA stiff','Reuss','Voigt','HSmin','HSmax','SC stiff','M-T')
title('Effect of spheres viewed by various EMA')
%yticks([0 0.2 0.4 0.6 0.8 1])
% xticks([0 0.2 0.4 0.6 0.8 1])
ylabel('K/K0')
xlabel('Porosity')
%ylim([0 1])
% xlim([0 1])
grid on

figure;
plot(n_total,K_ratio_NIA1,'green'); hold on
plot(n_total,K_ratio_NIA2,'blue'); hold on
plot(n_total,K_ratio_R,'black'); hold on
plot(n_total,K_ratio_V,'black--'); hold on
plot(n_total,K_ratio_HSmin,'red--'); hold on
plot(n_total,K_ratio_HSmax,'red'); hold on
plot(n_total,K_ratio_SC,'cyan'); hold on
plot(n_total,K_ratio_MT,'yellow--'); hold on
legend('NIA compl','NIA stiff','Reuss','Voigt','HSmin','HSmax','SC stiff','M-T')
title('Effect of spheres viewed by various EMA')
%yticks([0 0.2 0.4 0.6 0.8 1])
% xticks([0 0.2 0.4 0.6 0.8 1])
ylabel('K/K0')
xlabel('Number of cracks')
%ylim([0 1])
% xlim([0 1])
grid on

figure;
plot(e_total,K_ratio_NIA1,'green'); hold on
plot(e_total,K_ratio_NIA2,'blue'); hold on
plot(e_total,K_ratio_R,'black'); hold on
plot(e_total,K_ratio_V,'black--'); hold on
plot(e_total,K_ratio_HSmin,'red--'); hold on
plot(e_total,K_ratio_HSmax,'red'); hold on
plot(e_total,K_ratio_SC,'cyan'); hold on
plot(e_total,K_ratio_MT,'yellow--'); hold on
legend('NIA compl','NIA stiff','Reuss','Voigt','HSmin','HSmax','SC stiff','M-T')
title('Effect of spheres viewed by various EMA')
%yticks([0 0.2 0.4 0.6 0.8 1])
% xticks([0 0.2 0.4 0.6 0.8 1])
ylabel('K/K0')
xlabel('Crack density')
%ylim([0 1])
% xlim([0 1])
grid on

figure;
plot(p_total,E_ratio_NIA1,'green'); hold on
plot(p_total,E_ratio_NIA2,'blue'); hold on
plot(p_total,E_ratio_R,'black'); hold on
plot(p_total,E_ratio_V,'black--'); hold on
plot(p_total,E_ratio_HSmin,'red--'); hold on
plot(p_total,E_ratio_HSmax,'red'); hold on
plot(p_total,E_ratio_SC,'cyan'); hold on
plot(p_total,E_ratio_MT,'yellow--'); hold on
legend('NIA compl','NIA stiff','Reuss','Voigt','HSmin','HSmax','SC stiff','M-T')
title('Effect of spheres viewed by various EMA')
yticks([0 0.2 0.4 0.6 0.8 1])
% xticks([0 0.2 0.4 0.6 0.8 1])
ylabel('E/E0')
xlabel('Porosity')
ylim([0 1])
% xlim([0 1])
grid on

figure;
plot(n_total,E_ratio_NIA1,'green'); hold on
plot(n_total,E_ratio_NIA2,'blue'); hold on
plot(n_total,E_ratio_R,'black'); hold on
plot(n_total,E_ratio_V,'black--'); hold on
plot(n_total,E_ratio_HSmin,'red--'); hold on
plot(n_total,E_ratio_HSmax,'red'); hold on
plot(n_total,E_ratio_SC,'cyan'); hold on
plot(n_total,E_ratio_MT,'yellow--'); hold on
legend('NIA compl','NIA stiff','Reuss','Voigt','HSmin','HSmax','SC stiff','M-T')
title('Effect of spheres viewed by various EMA')
yticks([0 0.2 0.4 0.6 0.8 1])
% xticks([0 0.2 0.4 0.6 0.8 1])
ylabel('E/E0')
xlabel('Number of cracks')
ylim([0 1])
% xlim([0 1])
grid on

figure;
plot(e_total,E_ratio_NIA1,'green'); hold on
plot(e_total,E_ratio_NIA2,'blue'); hold on
plot(e_total,E_ratio_R,'black'); hold on
plot(e_total,E_ratio_V,'black--'); hold on
plot(e_total,E_ratio_HSmin,'red--'); hold on
plot(e_total,E_ratio_HSmax,'red'); hold on
plot(e_total,E_ratio_SC,'cyan'); hold on
plot(e_total,E_ratio_MT,'yellow--'); hold on
legend('NIA compl','NIA stiff','Reuss','Voigt','HSmin','HSmax','SC stiff','M-T')
title('Effect of spheres viewed by various EMA')
yticks([0 0.2 0.4 0.6 0.8 1])
% xticks([0 0.2 0.4 0.6 0.8 1])
ylabel('E/E0')
xlabel('Crack density')
ylim([0 1])
% xlim([0 1])
grid on



 E_ratio_NIA1(50)
 E_ratio_NIA2(50)
 E_ratio_SC(50)
 E_ratio_MT(50)
 E_ratio_HSmax(50)
 E_ratio_V(50)

%%%%%% STABILITY CONDITIONS %%%%%%%
% det(C_eff_V)
% C_eff_V2=[C_eff_V(1,1) C_eff_V(1,2) C_eff_V(1,3);
%           C_eff_V(2,1) C_eff_V(2,2) C_eff_V(2,3);
%           C_eff_V(3,1) C_eff_V(3,2) C_eff_V(3,3);];
% det(C_eff_V2)
% C_eff_V(1,1)*C_eff_V(2,2)-C_eff_V(1,2)*C_eff_V(2,1)
% C_eff_V(1,1)
% C_eff_V(2,2)
% C_eff_V(3,3)
% 
% det(C_eff_R)
% C_eff_R2=[C_eff_R(1,1) C_eff_R(1,2) C_eff_R(1,3);
%           C_eff_R(2,1) C_eff_R(2,2) C_eff_R(2,3);
%           C_eff_R(3,1) C_eff_R(3,2) C_eff_R(3,3);];
% det(C_eff_R2)
% C_eff_R(1,1)*C_eff_R(2,2)-C_eff_R(1,2)*C_eff_R(2,1)
% C_eff_R(1,1)
% C_eff_R(2,2)
% C_eff_R(3,3)

% det(C_eff1)
% C_eff11=[C_eff1(1,1) C_eff1(1,2) C_eff1(1,3);
%           C_eff1(2,1) C_eff1(2,2) C_eff1(2,3);
%           C_eff1(3,1) C_eff1(3,2) C_eff1(3,3);];
% det(C_eff11)
% C_eff1(1,1)*C_eff1(2,2)-C_eff1(1,2)*C_eff1(2,1)
% C_eff1(1,1)
% C_eff1(2,2)
% C_eff1(3,3)
% C_eff1(4,4)
% C_eff1(5,5)
% C_eff1(6,6)
 
% det(S_eff1)
% S_eff11=[S_eff1(1,1) S_eff1(1,2) S_eff1(1,3);
%           S_eff1(2,1) S_eff1(2,2) S_eff1(2,3);
%           S_eff1(3,1) S_eff1(3,2) S_eff1(3,3);];
% det(S_eff11)
% S_eff1(1,1)*S_eff1(2,2)-S_eff1(1,2)*S_eff1(2,1)
% S_eff1(1,1)
% S_eff1(2,2)
% S_eff1(3,3)
% S_eff1(4,4)
% S_eff1(5,5)
% S_eff1(6,6)

% det(C_eff2)
% C_eff22=[C_eff2(1,1) C_eff2(1,2) C_eff2(1,3);
%           C_eff2(2,1) C_eff2(2,2) C_eff2(2,3);
%           C_eff2(3,1) C_eff2(3,2) C_eff2(3,3);];
% det(C_eff22)
% C_eff2(1,1).*C_eff2(2,2)-C_eff2(1,2).*C_eff2(2,1)
% C_eff2(1,1)
% C_eff2(2,2)
% C_eff2(3,3)
% C_eff2(4,4)
% C_eff2(5,5)
% C_eff2(6,6)