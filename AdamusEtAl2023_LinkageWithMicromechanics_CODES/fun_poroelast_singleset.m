function [S, B] = fun_poroelast_singleset(H_c, phi, Kf, K)

%  %  %  %  %  %  
%  DESCRIPTION %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  % 
%  This function computes poroealstic coefficients of three vertical    %
%  subsets - having m pores each - forming a single connected set       %
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %


for i=1:length(H_c)
phi_H_c(:,:,i)=phi(i).*H_c(:,:,i);
end 

phi_p=sum(phi);                                                                   % volume fraction of a single pore set (equivalent to the total volume fraction)
H_p=sum(phi_H_c,3)./phi_p;                                                        % excess dry compliance of a single pore set 
invKd_p=H_p(1,1)+H_p(1,2)+H_p(1,3)+H_p(2,1)+H_p(2,2)+H_p(2,3)+...                 % inverse of a bulk modulus of the dry pore set
                                                   H_p(3,1)+H_p(3,2)+H_p(3,3);    


%%%_outputs_%%%

S=phi_p*(invKd_p+(1/Kf(1))-(1/K(1)));               % storage coefficient for a single connected set
B11=(3)*(phi_p/S)*(H_p(1,1)+H_p(1,2)+H_p(1,3));   
B22=(3)*(phi_p/S)*(H_p(2,1)+H_p(2,2)+H_p(2,3));
B33=(3)*(phi_p/S)*(H_p(3,1)+H_p(3,2)+H_p(3,3));
B12=(3)*(phi_p/S)*(H_p(6,1)+H_p(6,2)+H_p(6,3));
B13=(3)*(phi_p/S)*(H_p(5,1)+H_p(5,2)+H_p(5,3));
B23=(3)*(phi_p/S)*(H_p(4,1)+H_p(4,2)+H_p(4,3));
B=[B11 B12 B13; B12 B22 B23; B13 B23 B33];          % Skempton tensor for a single connected set
end

