% Author: Sukanta Sarkar et al
% Date: 07/11/2019

% Function Y_mu  


function s = Ymu(x2,bind_sites,mu0)

gammamui=[0.005 0.05 0.5 0.5 0.5 0.5];
  Ymu1=0;
  for i=1:bind_sites
      Ymu1=Ymu1+gammamui(i).*i.*nchoosek(bind_sites,i).*(x2./mu0).^i;      
  end
s=(1./(1+x2./mu0).^bind_sites).*Ymu1;