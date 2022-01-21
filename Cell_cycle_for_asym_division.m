%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Cell cycles required for first asymmetric division %%%%%%%%
%%%%%                      (Noise in SANIL)                %%%%%%%%
%%%%%           Authors: Jain et al. Date: 19/01/22        %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I_dau=zeros(2,1); % to store SNAIL values of daughter cells

p_name=['E','H','M']; % phenotype names

I= [100000 150000 189900 206600 219500 250000 300000]; % specific SNAIL levels

eta1 = [0.0 0.2 0.20 0.2 0.1 0.4]; % scaling factor for stochastic duplication of parent cell's SNAIL content
eta2 = [0.1 0.0 0.05 0.2 0.1 0.1]; % scaling factor for stochastic partitioning of parent cell's SNAIL content

total_obs=16; % total number of independent runs for each SNAIL value and phenotype pair

phenotype=zeros(1,3); % to store the phenotype of parent and two daughter cells  

gen_count_eta_I_pair=cell(length(I),1); % to store minimum cell cycles for asym division for every phenotype possible at a given SNAIL value 

gen_count_eta_I_pair_initial_phenotype = cell(length(I),1); % to identify which phenotypes are stable at a given SNAIL value

flag=0; % raised when a asymmetric division is observed
y_steady_dau=zeros(4,2); % steady levels of miR200, mZEB, ZEB and SNAIL for the two daughter cells

%cd('/home/csb/CollaborationWork-Paras/Population_Dynamics/Final figures, data and codes/Mean First Passage time/Codes');


load 'equibrium_pts_and_eigenvalues_EMT_reduced.mat' % to initialize a cells' state/variables,
 %%% contains 'x' variable storing steady levels of miR200, mZEB, ZEB and SNAIL and 'f' variable storing their corresponding 
 %%% eigen values for various SNAIL levels

for eta_id  = 1:length(eta1)
        
        cd('/home/csb/CollaborationWork-Paras/Population_Dynamics/Final figures, data and codes/Mean First Passage time/Codes');
        
        
        for indx=1:length(I) 
                %%% to initialize a cell state/variables for all possible
                %%% phentype (one at a time) at a given SNAIL value
                first_time=1;
                y_steady=[];
                
                for g =1:size(x,2)
                    x1 = find(f(:,g) > 0);
                    if(abs(x(4,g)-I(indx))<=100 && isempty(x1))
                        if(first_time==1)
                            diff=abs(x(4,g)-I(indx));
                            min_indx=g;
                            first_time=2;
                        else
                            if(abs(x(4,g)-I(indx))<diff)
                                min_indx=g;
                            end
                        end
                    elseif(first_time==2)
                        y_steady=[y_steady,x(:,min_indx)];
                        first_time=1;
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%% for each stable phenotype at a given SNAIL level run the analysis 
                for n_states = 1:size(y_steady,2)
                    
                    gen_count_eta_I_pair_initial_phenotype{indx}(n_states)=get_phenotype(y_steady(2,n_states),y_steady(end,n_states));
                    
                    obs_count=1;
                    
                    gen_count=ones(total_obs,1);
                    
                    %%% to perform multiple runs of same analysis
                    while(obs_count<=total_obs)
                        fprintf(['\n I=' num2str(I(indx)) ' n_state = ' num2str(n_states) ' eta1_=' num2str(eta1(eta_id)) ' eta2_=' num2str(eta2(eta_id)) ' obser. no.=' num2str(obs_count) '\n'])
                        dau_indx_prev=1;
                        y_steady_par = zeros(4,5000);
                        y_steady_par(:,dau_indx_prev)=y_steady(:,n_states);
                        
                        while(true)
                            dau_indx=0; %% dau_indx account for cells in the next generation
                            for z=1:dau_indx_prev %% dau_indx_prev account for cells in current generation
                                phenotype(1)=get_phenotype(y_steady_par(2,z),y_steady_par(4,z));
                                rand_num = randn(1,2);
                                for i=1:2
                                    
                                    if(i == 1)
                                        I_dau(i)=y_steady(4,n_states)+y_steady(4,n_states)*eta1(eta_id)*rand_num(1,1)/2 + ...
                                            eta2(eta_id)*(2*y_steady(4,n_states)+y_steady(4,n_states)*eta1(eta_id)*rand_num(1,1))*rand_num(1,2);
                                    else
                                        I_dau(i)=y_steady(4,n_states)+y_steady(4,n_states)*eta1(eta_id)*rand_num(1,1)/2 - ...
                                            eta2(eta_id)*(2*y_steady(4,n_states)+y_steady(4,n_states)*eta1(eta_id)*rand_num(1,1))*rand_num(1,2);
                                    end
                                    
                                    if(I_dau(i)<0)
                                        I_dau(i)=0;
                                    end
                                    
                                    %%% finding the new steady state of the
                                    %%% daugther cellular content
                                    y_steady_dau(:,i)=steady_points([y_steady_par(1:3,z);I_dau(i)]);
                                    
                                    %%% assigning phenotype to daughter cells
                                    %%% based upon their cellular content
                                    phenotype(i+1)=get_phenotype(y_steady_dau(2,i),I_dau(i));
                                end
                                
                                %%% checking for asymmetric division, if
                                %%% not happened for a generation of cells then simulate the cell division of next
                                %%% generation and check again, otherwise
                                %%% stop current run, save cell generation count and go to next iteration. 
                                if(phenotype(2)==phenotype(1) && phenotype(3)==phenotype(1))
                                    y_steady_par(:,z)=y_steady_dau(:,1);
                                    y_steady_par(:,dau_indx_prev+z)=y_steady_dau(:,2);
                                    dau_indx=dau_indx+2;
                                else
                                    flag=1;
                                    break
                                end
                            end
                            
                            if(flag==1)
                                break;
                            else
                                gen_count(obs_count)=gen_count(obs_count)+1;
                                if(gen_count(obs_count) > 12) %%% restricting our serach upto 12 cell cycles/generations
                                    break;
                                end
                                dau_indx_prev=dau_indx;
                            end
                        end
                        flag=0;
                        fprintf(['gen count=' num2str(gen_count(obs_count)) '\n'])
                        obs_count=obs_count+1;
                    end
                    gen_count_eta_I_pair{indx}(n_states,:)=gen_count;
                end
        end
        
        %cd('/home/csb/CollaborationWork-Paras/Population_Dynamics/Final figures, data and codes/Dependent_formalism/First asymmetric division')
        save(['gen_count_eta_I_pair_eta1_' num2str(eta1(eta_id)) '_eta2_' num2str(eta2(eta_id)) '.mat'],'gen_count_eta_I_pair','gen_count_eta_I_pair_initial_phenotype');
        
   
end