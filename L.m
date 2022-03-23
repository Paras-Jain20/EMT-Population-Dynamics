% Author: Sukanta Sarkar et al
% Date: 07/11/2019

% Function L  


function s = L(x2,bind_sites,mu0)

li=[0.6,0.3,0.1,0.05,0.05,0.05];
L1=1.0;
for i=1:bind_sites
L1=L1+li(i).*nchoosek(bind_sites,i).*(x2./mu0).^i;
end
s=(1./(1+x2./mu0).^bind_sites).*L1;