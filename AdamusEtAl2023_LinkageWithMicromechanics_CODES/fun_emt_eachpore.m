function [H_c, phi, Kf, K]=fun_emt_eachpore(m, y, theta_p1, theta_p2, theta_p3, Kf, E, v)

%  %  %  %  %  %  
% DESCRIPTION  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  % 
% This function computes dry micromechanical coefficients using EMT  %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  

x=3*m;                      % total number of pores
k=1;                        % extra dimensionality (keep k=1 here)
Kf=Kf*ones(x,k);            % fluid compressibility
E=E*ones(x,k);              % Young's modulus of a background solid
v=v*ones(x,k);              % Poisson's ratio of a background solid
l=E.*v./((1+v).*(1-2.*v));  % Lame' parameter
G=E./(2+2.*v);              % rigidity
K=E./(3.*(1-2.*v));         % bulk modulus
 

for j=1:k
for i=1:x
Cb(:,:,i,j)= [ l(i,j)+2.*G(i,j), l(i,j), l(i,j), 0, 0, 0;    % background elasticity matrix
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

%%%_shape_%%%

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

 % 6x6 matrix representation of the Eshelby tensor for horizontally oriented spheroids
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

%%%%%%%%%%%%%%
%%%_orient_%%%
%%%%%%%%%%%%%%  
  
A = [1,0,0;0,0,1;0,1,0;];   % get 3x3 rotation matrix (to first get the vertical orientation of pores)
A11=A(1,1); A12=A(1,2); A13=A(1,3); A21=A(2,1); A22=A(2,2); A23=A(2,3); A31=A(3,1); A32=A(3,2); A33=A(3,3); 

% get 6x6 rotation matrix to rotate the 6x6 matrix representation of the Eshelby tensor
MAA= [A11.*A11 A12.*A12 A13.*A13 A12.*A13 A11.*A13 A11.*A12;
     A21.*A21 A22.*A22 A23.*A23 A22.*A23 A21.*A23 A21.*A22;
     A31.*A31 A32.*A32 A33.*A33 A32.*A33 A31.*A33 A31.*A32;
     2.*A21.*A31 2.*A22.*A32 2.*A23.*A33 A22.*A33+A23.*A32 A21.*A33+A23.*A31 A21.*A32+A22.*A31;
     2.*A11.*A31 2.*A12.*A32 2.*A13.*A33 A12.*A33+A13.*A32 A11.*A33+A13.*A31 A11.*A32+A12.*A31;
     2.*A11.*A21 2.*A12.*A22 2.*A13.*A23 A12.*A23+A13.*A22 A11.*A23+A13.*A21 A11.*A22+A12.*A21;];

for j=1:k
for i=1:x 
Erot0(:,:,i,j)=MAA'*e(:,:,i,j)*MAA;  % 6x6 matrix representation of the Eshelby tensor describes vertical pores now (instead of horizontal)
 
x1=zeros(1,x,k); x2=zeros(1,x,k); x3=ones(1,x,k); x4=cat(1,x1,x2,x3);
x6=theta_p1*ones(m,k); x7=theta_p2*ones(m,k); x8=theta_p3*ones(m,k);
x5=cat(1,x6,x7,x8); w=x4; theta=x5;             % angular distribution of three vertical subsets m pores each

% Rotate the vertical pores to form three subsets distributed according to theta
bw(:,:,i,j) = [0, -w(3,i,j), w(2,i,j); w(3,i,j), 0, -w(1,i,j); -w(2,i,j), w(1,i,j), 0];
A(:,:,i,j) = eye(3) + sind(theta(i,j))*bw(:,:,i,j) + (1-cosd(theta(i,j)))*bw(:,:,i,j)*bw(:,:,i,j);

A11(i,j)=A(1,1,i,j); A12(i,j)=A(1,2,i,j); A13(i,j)=A(1,3,i,j); A21(i,j)=A(2,1,i,j); A22(i,j)=A(2,2,i,j); A23(i,j)=A(2,3,i,j); A31(i,j)=A(3,1,i,j); A32(i,j)=A(3,2,i,j); A33(i,j)=A(3,3,i,j);   
MA(:,:,i,j)= [A11(i,j).*A11(i,j) A12(i,j).*A12(i,j) A13(i,j).*A13(i,j) A12(i,j).*A13(i,j) A11(i,j).*A13(i,j) A11(i,j).*A12(i,j);
     A21(i,j).*A21(i,j) A22(i,j).*A22(i,j) A23(i,j).*A23(i,j) A22(i,j).*A23(i,j) A21(i,j).*A23(i,j) A21(i,j).*A22(i,j);
     A31(i,j).*A31(i,j) A32(i,j).*A32(i,j) A33(i,j).*A33(i,j) A32(i,j).*A33(i,j) A31(i,j).*A33(i,j) A31(i,j).*A32(i,j);
     2.*A21(i,j).*A31(i,j) 2.*A22(i,j).*A32(i,j) 2.*A23(i,j).*A33(i,j) A22(i,j).*A33(i,j)+A23(i,j).*A32(i,j) A21(i,j).*A33(i,j)+A23(i,j).*A31(i,j) A21(i,j).*A32(i,j)+A22(i,j).*A31(i,j);
     2.*A11(i,j).*A31(i,j) 2.*A12(i,j).*A32(i,j) 2.*A13(i,j).*A33(i,j) A12(i,j).*A33(i,j)+A13(i,j).*A32(i,j) A11(i,j).*A33(i,j)+A13(i,j).*A31(i,j) A11(i,j).*A32(i,j)+A12(i,j).*A31(i,j);
     2.*A11(i,j).*A21(i,j) 2.*A12(i,j).*A22(i,j) 2.*A13(i,j).*A23(i,j) A12(i,j).*A23(i,j)+A13(i,j).*A22(i,j) A11(i,j).*A23(i,j)+A13(i,j).*A21(i,j) A11(i,j).*A22(i,j)+A12(i,j).*A21(i,j);];

Erot(:,:,i,j)=MA(:,:,i,j)'*Erot0(:,:,i,j)*MA(:,:,i,j);     % matrix representation of the Eshelby tensor describing rotated vertical cracks
H_c(:,:,i,j)=inv(Cb(:,:,i,j)-F*Cb(:,:,i,j)*Erot(:,:,i,j)); % excess dry compliance of each pore
end
end

%%%_size_%%%

for i=1:x
d(i)=1./(9*m);                          % crack density parameter being a^3/V where 'a' is the radius of a cricular pore and V is the material's volume
phi(i)=(1./3).*4.*pi.*y.*d(i);          % volume fraction of a single pore         
end
