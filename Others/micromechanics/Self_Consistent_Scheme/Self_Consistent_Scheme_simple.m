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
Kf=1.0001/0.45; Gf=0.0001;     % min 0GPa (air) max 8GPa (ice at zero deg) 

% first inclusion type
for jj=1:40
n(1)=1*jj;                     % number of defined inclusions
a(1)=4;                     % radius of circular inclusion (in mm)
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

%%% LOOP - y %%%

for i=1:2
  
p(i)=(4/3).*(pi.*y(i).*a(i).^3)./vol;  % volume fraction of a defined pore

end

p_total(jj)=n(1).*p(1)+n(2).*p(2);
e_total(jj)=(n(1)*a(1).^3./vol)+(n(2)*a(2).^3./vol);
n_total(jj)=n(1)+n(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  RESULTS  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SC Kachanov-Sevostianov %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:300000
v_eff(j)=rand-0.5;
G_eff(j)=rand*G;

K_eff(j)=(2.*G_eff(j).*(1+v_eff(j)))./(3.*(1-2.*v_eff(j)));
alpha=3*K_eff(j)/(3*K_eff(j)+4*G_eff(j));
beta=(6/5)*(K_eff(j)+2.*G_eff(j))./(3.*K_eff(j)+2.*G_eff(j));

K_SC=K_eff(j)-p_total(jj).*K_eff(j).*(Kf-K)./(K_eff(j)+alpha.*(Kf-K_eff(j)));
G_SC=G_eff(j)-p_total(jj).*G_eff(j).*(Gf-G)./(G_eff(j)+beta.*(Gf-G_eff(j)));

err1(j)=(K_SC-K)^2;
err2(j)=(G_SC-G)^2;
err(j)=sqrt(err1(j)+err2(j));
end
[row1] = find (err==min(err));


K_true=K_eff(row1);
G_true=G_eff(row1);
E_true=9*K_true*G_true/(3*K_true+G_true);

K_ratio_SC1(jj)=K_true./K;
E_ratio_SC1(jj)=E_true./E;
G_ratio_SC1(jj)=G_true./G;
end
%%%%%%%%%%%%%%%%%%%%%%%%
min(err1)
min(err2)
min(err)
%%
figure; 
% plot(p_total,K_ratio_NIA1,'green'); hold on
% plot(p_total,K_ratio_NIA2,'blue'); hold on
 plot(p_total,K_ratio_SC1,'black'); hold on
% plot(p_total,K_ratio_SC2,'black--'); hold on
% legend('NIA compl','NIA stiff','SC1')
title('Effect of spheres viewed by dual NIA')
%yticks([0 0.2 0.4 0.6 0.8 1])
% xticks([0 0.2 0.4 0.6 0.8 1])
ylabel('K/K0')
xlabel('Porosity')
%ylim([0 1])
% xlim([0 1])
grid on

figure; 
% plot(p_total,K_ratio_NIA1,'green'); hold on
% plot(p_total,K_ratio_NIA2,'blue'); hold on
plot(p_total,E_ratio_SC1,'black'); hold on
% plot(p_total,E_ratio_SC2,'black--'); hold on
% legend('NIA compl','NIA stiff','SC1')
title('Effect of spheres viewed by dual NIA')
%yticks([0 0.2 0.4 0.6 0.8 1])
% xticks([0 0.2 0.4 0.6 0.8 1])
ylabel('E/E0')
xlabel('Porosity')
%ylim([0 1])
% xlim([0 1])
grid on

figure; 
% plot(p_total,K_ratio_NIA1,'green'); hold on
% plot(p_total,K_ratio_NIA2,'blue'); hold on
plot(p_total,G_ratio_SC1,'black'); hold on
% plot(p_total,G_ratio_SC2,'black--'); hold on
% legend('NIA compl','NIA stiff','SC1')
title('Effect of spheres viewed by dual NIA')
%yticks([0 0.2 0.4 0.6 0.8 1])
% xticks([0 0.2 0.4 0.6 0.8 1])
ylabel('G/G0')
xlabel('Porosity')
%ylim([0 1])
% xlim([0 1])
grid on