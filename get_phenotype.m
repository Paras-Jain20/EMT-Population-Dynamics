%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                 Phenotype determination              %%%%%%%%
%%%(Based on cell's position in attractor basin of an phenotype %%%
%%%                 as per mZEB-SNAIL bifurcation plot)        %%%%
%%%  Method adopted from Tripathi et al. 2020 PLOS Comp Bio  %%%%%%
%%%%%           Authors: Jain et al. Date: 19/01/22        %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function phenotype = get_phenotype(m,I)


 x1 = 193.2; x2 = 208.7; y1 = 243.09; y2 = 90.93;
 a1 = y2 - y1; b1 = x1 - x2; c1 = x2*y1 - y2*x1;

 u1 = 185.12; u2 = 224.67; v1 = 698.395; v2 = 495.802;
 a2 = v2 - v1; b2 = u1 - u2; c2 = u2*v1 - v2*u1;

phenotype =-1;
x=I/10^3;
y=m;

    if(x < u1)

        phenotype = 0;

    elseif(x > u2)

        phenotype = 2;

    else
        fac1 = (a1*x + b1*y + c1) / b1;
        fac2 = (a2*x + b2*y + c2) / b2;
        if(fac2 >= 0)
            phenotype = 2;

        elseif(fac1 < 0)
            phenotype = 0;
        else
            phenotype = 1;
        end
    end
    
   if(phenotype==-1)
       fprintf("phenotype not properly assigned\n");
   end

end

