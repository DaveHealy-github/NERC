function [S_p1, S_p2, S_p3, B_p1, B_p2, B_p3] = fun_poroelast_isolatedsets(H_c, phi, Kf, K)

%  %  %  %  %  %  
%  DESCRIPTION %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  % 
%  This function computes poroealstic coefficients of three vertical    %
%  sets - having m pores each - that are isolated from each other       %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %


x=length(H_c);
for i=1:x
phi_H_c(:,:,i)=phi(i).*H_c(:,:,i);
end 

phi_p1=sum(phi(1:(x/3)));               % volume fraction of a pore set - first set
phi_p2=sum(phi(1+(x/3):(2*x/3)));       % volume fraction of a pore set - second set
phi_p3=sum(phi(1+(2*x/3):x));           % volume fraction of a pore set - third set

H_p1=sum(phi_H_c(:,:,1:(x/3)),3)./phi_p1;              % excess dry compliance of a pore set - first set
H_p2=sum(phi_H_c(:,:,1+(x/3):2*x/3),3)./phi_p2;        % excess dry compliance of a pore set - second set
H_p3=sum(phi_H_c(:,:,1+(2*x/3):x),3)./phi_p3;          % excess dry compliance of a pore set - third set

invKd_p1=H_p1(1,1)+H_p1(1,2)+H_p1(1,3)+H_p1(2,1)+H_p1(2,2)+H_p1(2,3)+H_p1(3,1)+H_p1(3,2)+H_p1(3,3);   % inverse of a bulk modulus of the dry pore set - first set
invKd_p2=H_p2(1,1)+H_p2(1,2)+H_p2(1,3)+H_p2(2,1)+H_p2(2,2)+H_p2(2,3)+H_p2(3,1)+H_p2(3,2)+H_p2(3,3);   % inverse of a bulk modulus of the dry pore set - second set
invKd_p3=H_p3(1,1)+H_p3(1,2)+H_p3(1,3)+H_p3(2,1)+H_p3(2,2)+H_p3(2,3)+H_p3(3,1)+H_p3(3,2)+H_p3(3,3);   % inverse of a bulk modulus of the dry pore set - third set

%%%_outputs_%%%

S_p1=phi_p1*(invKd_p1+(1/Kf(1))-(1/K(1)));                                 % Storage coefficient - first set
S_p2=phi_p2*(invKd_p2+(1/Kf(1))-(1/K(1)));                                 % Storage coefficient - second set
S_p3=phi_p3*(invKd_p3+(1/Kf(1))-(1/K(1)));                                 % Storage coefficient - third set

B11_p1=(3)*(phi_p1/S_p1)*(H_p1(1,1)+H_p1(1,2)+H_p1(1,3));
B22_p1=(3)*(phi_p1/S_p1)*(H_p1(2,1)+H_p1(2,2)+H_p1(2,3));
B33_p1=(3)*(phi_p1/S_p1)*(H_p1(3,1)+H_p1(3,2)+H_p1(3,3));
B12_p1=(3)*(phi_p1/S_p1)*(H_p1(6,1)+H_p1(6,2)+H_p1(6,3));
B13_p1=(3)*(phi_p1/S_p1)*(H_p1(5,1)+H_p1(5,2)+H_p1(5,3));
B23_p1=(3)*(phi_p1/S_p1)*(H_p1(4,1)+H_p1(4,2)+H_p1(4,3));

B11_p2=(3)*(phi_p2/S_p2)*(H_p2(1,1)+H_p2(1,2)+H_p2(1,3));
B22_p2=(3)*(phi_p2/S_p2)*(H_p2(2,1)+H_p2(2,2)+H_p2(2,3));
B33_p2=(3)*(phi_p2/S_p2)*(H_p2(3,1)+H_p2(3,2)+H_p2(3,3));
B12_p2=(3)*(phi_p2/S_p2)*(H_p2(6,1)+H_p2(6,2)+H_p2(6,3));
B13_p2=(3)*(phi_p2/S_p2)*(H_p2(5,1)+H_p2(5,2)+H_p2(5,3));
B23_p2=(3)*(phi_p2/S_p2)*(H_p2(4,1)+H_p2(4,2)+H_p2(4,3));

B11_p3=(3)*(phi_p3/S_p3)*(H_p3(1,1)+H_p3(1,2)+H_p3(1,3));
B22_p3=(3)*(phi_p3/S_p3)*(H_p3(2,1)+H_p3(2,2)+H_p3(2,3));
B33_p3=(3)*(phi_p3/S_p3)*(H_p3(3,1)+H_p3(3,2)+H_p3(3,3));
B12_p3=(3)*(phi_p3/S_p3)*(H_p3(6,1)+H_p3(6,2)+H_p3(6,3));
B13_p3=(3)*(phi_p3/S_p3)*(H_p3(5,1)+H_p3(5,2)+H_p3(5,3));
B23_p3=(3)*(phi_p3/S_p3)*(H_p3(4,1)+H_p3(4,2)+H_p3(4,3));

B_p1=[B11_p1 B12_p1 B13_p1; B12_p1 B22_p1 B23_p1; B13_p1 B23_p1 B33_p1];   % Skempton tensor - first set
B_p2=[B11_p2 B12_p2 B13_p2; B12_p2 B22_p2 B23_p2; B13_p2 B23_p2 B33_p2];   % Skempton tensor - second set
B_p3=[B11_p3 B12_p3 B13_p3; B12_p3 B22_p3 B23_p3; B13_p3 B23_p3 B33_p3];   % Skempton tensor - third set
