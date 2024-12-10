function FF = root2d(x,jj)
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
% for jj=1:100
% n(1)=1*jj;                     % number of defined inclusions
a(1)=3;                     % radius of circular inclusion (in mm)
y(1)=1;                    % aspect ratio of inclusion
or(1)=0;                    % orientation of inclusion 0=horizontal 1=vertical

% second inclusion type
% n(2)=0;                     % number of defined inclusions
% a(2)=2;                     % radius of circular inclusion (in mm)
% y(2)=2;                    % aspect ratio of inclusion
% or(2)=1;                    % orientation of inclusion 0=horizontal 1=vertical

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
% for i=1:2
% if y(i)>1   
%  s(i)=1./(1-y(i).^2)-y(i)./(1-y(i).^2).*(1/(sqrt(y(i).^2-1))).*acosh(y(i));
%  t=0.5*(1-s(i));
% psi1=(y(i).^2.*(4*t-1)-t)./(4*(1-y(i).^2));
% psi2=(y(i).^2.*(1-2*t)-t)./(4*(1-y(i).^2));
% psi4=(3.*t-1)./(2*(1-y(i).^2));
%   p1=(1./x(1)).*(t+psi1./(1-x(2)));
%   p2=psi2./(x(1).*(1-x(2)));
%   p3=p2;
%   p4=(1./x(1)).*(1-2*t+psi4./(1-x(2)));
%   p5=(1./x(1)).*(t+psi1./(2*(1-x(2))));
%   p6=(1./x(1)).*(0.5*(1-t)+2*psi2./(1-x(2)));
%  P_eff=p1*H1+p2*H2+p3*H3+p4*H4+p5*H5+p6*H6; 
% elseif y(i)<1
%  s(i)=1./(1-y(i).^2)-y(i)./(1-y(i).^2).*(1/(sqrt(1-y(i).^2))).*acos(y(i));
%  t=0.5*(1-s(i));
% psi1=(y(i).^2.*(4*t-1)-t)./(4*(1-y(i).^2));
% psi2=(y(i).^2.*(1-2*t)-t)./(4*(1-y(i).^2));
% psi4=(3.*t-1)./(2*(1-y(i).^2));
%   p1=(1./x(1)).*(t+psi1./(1-x(2)));
%   p2=psi2./(x(1).*(1-x(2)));
%   p3=p2;
%   p4=(1./x(1)).*(1-2*t+psi4./(1-x(2)));
%   p5=(1./x(1)).*(t+psi1./(2.*(1-x(2))));
%   p6=(1./x(1)).*(0.5*(1-t)+2.*psi2./(1-x(2)));
%  P_eff=p1*H1+p2*H2+p3*H3+p4*H4+p5*H5+p6*H6; 
% else
%   s(i)=1/3;
%   p1=(1-2.*x(2))./(6.*x(1).*(1-x(2)));
%   p2=(4-5.*x(2))./(15.*x(1).*(1-x(2)));
%   P_eff=p1*I1+p2*I2;
% end
% 
%  Prot_eff=MB'*P_eff*MB;
% 
% if or(i)==0
%  N_eff(:,:,i)=((C1-Cb)^(-1)+P_eff)^(-1);
%  else
%  N_eff(:,:,i)=((C1-Cb)^(-1)+Prot_eff)^(-1);
%  end
% 
%  Q_eff=Cb-F*Cb*(F*P_eff*Cb);
% 
% if or(i)==0
% H_eff(:,:,i)=inv((S1-Sb)^(-1)+Q_eff);
% else
% H_eff(:,:,i)=inv((S1-Sb)^(-1)+Cb-F*Cb*(Prot_eff)*F*Cb);
% end
% 
% p(i)=(4/3).*(pi.*y(i).*a(i).^3)./vol;  % volume fraction of a defined pore
% end
% 
% p_total(jj)=n(1).*p(1)+n(2).*p(2);
% e_total(jj)=(n(1)*a(1).^3./vol)+(n(2)*a(2).^3./vol);
% n_total(jj)=n(1)+n(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  RESULTS  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SC compliance %%%%
%%%%%%%%%%%%%%%%%%%%%%%%
S_eff=[1./(2.*x(1).*(1+x(2)))           -x(2)./(2.*x(1).*(1+x(2)))     -x(2)./(2.*x(1).*(1+x(2))) 0 0 0;
      -x(2)./(2.*x(1).*(1+x(2)))    1./(2.*x(1).*(1+x(2)))             -x(2)./(2.*x(1).*(1+x(2)))  0 0 0;
      -x(2)./(2.*x(1).*(1+x(2)))    -x(2)./(2.*x(1).*(1+x(2)))     1./(2.*x(1).*(1+x(2))) 0 0 0;
      0 0 0 1./x(1) 0 0;
      0 0 0 0 1./x(1) 0;
      0 0 0 0 0 1./x(1)];

p(1)=(4/3).*(pi.*y(1).*a(1).^3)./vol;
H_eff(:,:)=inv((S1-Sb)^(-1)+Cb-F*Cb*(F*(I1*(1-2.*x(2))./(6.*x(1).*(1-x(2)))+I2*(4-5.*x(2))./(15.*x(1).*(1-x(2))))*Cb));

FF(1)=-Sb(1,1)+S_eff(1,1)-jj*p(1)*(H_eff(1,1));
FF(2)=-Sb(2,2)+S_eff(2,2)-jj*p(1)*(H_eff(2,2));
% FF(3)=-Sb(3,3)+S_eff(3,3)-n(1)*p(1)*(H_eff(3,3));
% FF(4)=-Sb(1,2)+S_eff(1,2)-n(1)*p(1)*(H_eff(1,2));
% FF(5)=-Sb(1,3)+S_eff(1,3)-n(1)*p(1)*(H_eff(1,3));
% FF(6)=-Sb(2,3)+S_eff(2,3)-n(1)*p(1)*(H_eff(2,3));
% FF(7)=-Sb(4,4)+S_eff(4,4)-n(1)*p(1)*(H_eff(4,4));
% FF(8)=-Sb(5,5)+S_eff(5,5)-n(1)*p(1)*(H_eff(5,5));
% FF(9)=-Sb(6,6)+S_eff(6,6)-n(1)*p(1)*(H_eff(6,6));

