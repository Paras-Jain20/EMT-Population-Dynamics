% Author: Sukanta Sarkar et al
% Date: 07/11/2019

% Function Y_m 


function s = Ym(x2,bind_sites,mu0)

gammami=[0.04 0.2 1 1 1 1];
Ym1=0.0;
for i=1:bind_sites
Ym1=Ym1+gammami(i).*nchoosek(bind_sites,i).*(x2./mu0).^i;
end
s=(1./(1+x2./mu0).^bind_sites).*Ym1;