%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                     Cells' Steady State              %%%%%%%%
%%%%%           Authors: Jain et al. Date: 19/01/22        %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function y_steady=steady_points(y_init)
      
        while(true)
            [t,y]=ode45(@(t,y) Fx(t,y),[0 10],y_init);
            if(abs(y(end-4,2)-y(end,2))<0.0001)
                break;
            else
                y_init=y(end,:);
                
            end
        end
        y_steady=y(end,:);
       
end


function dxx=Fx(t,x)

        mu0=10000;Zmu0=220000;Zm0=25000;Smu0=180000;Sm0=180000;
        gmu=2100; gm=11; kmu=0.05; km=0.5;gZ=100;kZ=0.1; n_mu=6;
        nZ_mu=3;lambda_Z_mu=0.1;n2=2;lambda_Z_m=7.5;lambda_S_m=10.0;
        nZ_m=n2;nS_m=n2;nS_mu=n2;lambda_S_mu=lambda_Z_mu;

        mu=x(1);
        m=x(2);
        Z=x(3);
        S=x(4);
        
        hZ_mu=H(Z,nZ_mu,Zmu0,lambda_Z_mu);
        hS_mu=H(S,nS_mu,Smu0,lambda_S_mu);
        hZ_m=H(Z,nZ_m,Zm0,lambda_Z_m);
        hS_m=H(S,nS_m,Sm0,lambda_S_m);
       
        
        h1=hZ_mu.*hS_mu;
        h2=hZ_m.*hS_m;
        dxx=[h1*gmu-m*Ymu(mu,n_mu,mu0)-mu*kmu; h2*gm-m*Ym(mu,n_mu,mu0)-m*km; m*gZ*L(mu,n_mu,mu0)-Z*kZ;0];

end