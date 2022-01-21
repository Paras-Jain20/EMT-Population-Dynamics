%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Probability of undergoing Sym, Asym and Div division %%%%%%%%
%%%%%                  (Noise in SNAIL or ALL)             %%%%%%%%
%%%%%           Authors: Jain et al. Date: 19/01/22        %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I_dau=zeros(2,1); % to store SNAIL values of daughter cells

p_name=['E','H','M']; % phenotype names

I= [100000 150000 189900 206600 219500 250000 300000];
total_obs=100; % total iterations (cell divisions) in a run

phenotype=zeros(total_obs,3); % to store the phenotype of parent and two daughter cells in every iteration 

num_runs = 10; % total independent runs of 100 iterations (cell divisions) each

eta1 = [0.2 0.0 0.2 0.20 0.2 0.1 0.4]; % scaling factor for stochastic duplication of parent cell's SNAIL content
eta2 = [0.1 0.1 0.0 0.05 0.2 0.1 0.1]; % scaling factor for stochastic partitioning of parent cell's SNAIL content

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noise_pert = 2; %% 1 : only in snail, 2: in all nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


y_perturbed_dau = zeros(4,2); %  levels of miR200, mZEB, ZEB and SNAIL for the just born daughter cells

%cd('/home/csb/CollaborationWork-Paras/Population_Dynamics/Final figures, data and codes/Mean First Passage time/Codes');
load 'equibrium_pts_and_eigenvalues_EMT_reduced.mat' %%% contains 'x' variable storing steady levels of miR200, mZEB, ZEB and SNAIL and 'f' variable storing their corresponding 
                                                     %%% eigen values for various SNAIL levels

for eta_id  = 1:length(eta2)

        %cd('/home/csb/CollaborationWork-Paras/Population_Dynamics/Final figures, data and codes/Dependent_formalism/Division type probabilities')
        
        % to check the presence of output file and account of data already
        % present in it. Accordingly, setting the required number of runs
        if(isfile(['Multiple_runs_Division_type_prob_data_eta1_' num2str(eta1(eta_id)) '_eta2_' num2str(eta2(eta_id)) '.csv']))
            data = readtable(['Multiple_runs_Division_type_prob_data_eta1_' num2str(eta1(eta_id)) '_eta2_' num2str(eta2(eta_id)) '.csv']);
            num_runs_req = num_runs - height(data)/(length(I)*3*3);
           
            probability = zeros(length(I)*3*3*num_runs_req,6);
        else
            num_runs_req = num_runs;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %cd('/home/csb/CollaborationWork-Paras/Population_Dynamics/Final figures, data and codes/Mean First Passage time/Codes');
        
        for sim_num = 1:num_runs_req
            
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
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%% run the division type analysis for each stable phenotype at a given SNAIL level
                for n_states = 1:size(y_steady,2)
                    
                    obs_count=1;
                    sym_count = 0;
                    asym_count = 0;
                    div_count = 0;
                    
                    if(noise_pert == 1)
                        rand_num = randn(total_obs,2);
                    else
                        rand_num = randn(total_obs,2,4);
                    end
                    
                    %%% iterating over same SNAIL level
                    while(obs_count<=total_obs)
%                         disp(obs_count)
                        phenotype(obs_count,1) = get_phenotype(y_steady(2,n_states),y_steady(4,n_states));
                        for i=1:2
                            
                            %%% assigning SNAIl levels to daughter cells by
                            %%% acounting noise/fluctuations in duplication
                            %%% and partitioning of parent's SNAIL content
                            if(noise_pert == 2)
                                
                                if(i == 1)
                                    y_perturbed_dau(:,i)=y_steady(:,n_states) + y_steady(:,n_states) .* eta1(eta_id) .* reshape(rand_num(obs_count,1,:),4,1)/2 +...
                                        eta2(eta_id) * (2*y_steady(:,n_states) + y_steady(:,n_states) .* eta1(eta_id) .* reshape(rand_num(obs_count,1,:),4,1)) .*reshape(rand_num(obs_count,2,:),4,1);
                                else
                                    y_perturbed_dau(:,i)=y_steady(:,n_states)+y_steady(:,n_states).*eta1(eta_id).*reshape(rand_num(obs_count,1,:),4,1)/2 -...
                                        eta2(eta_id) * (2*y_steady(:,n_states)+y_steady(:,n_states).*eta1(eta_id).*reshape(rand_num(obs_count,1,:),4,1)).*reshape(rand_num(obs_count,2,:),4,1);
                                end
                                
                                if(~isempty(find(y_perturbed_dau(:,i)<0,1)))
                                    y_perturbed_dau(y_perturbed_dau(:,i)<0,i) = 0;
                                end
                                
                                %%% finding the new steady state of the
                                %%% daugther cellular content
                                y_steady_dau=steady_points(y_perturbed_dau(:,i));
                                
                                %%% assigning phenotype to daughter cells
                                %%% based upon their cellular content
                                phenotype(obs_count,i+1)=get_phenotype(y_steady_dau(end,2),y_steady_dau(end,4));
                                
                            elseif(noise_pert == 1)
                                if(i == 1)
                                    I_dau(i)=y_steady(4,n_states)+y_steady(4,n_states)*eta1(eta_id)*rand_num(obs_count,1)/2 + ...
                                        eta2(eta_id)*(2*y_steady(4,n_states)+y_steady(4,n_states)*eta1(eta_id)*rand_num(obs_count,1))*rand_num(obs_count,2);
                                else
                                    I_dau(i)=y_steady(4,n_states)+y_steady(4,n_states)*eta1(eta_id)*rand_num(obs_count,1)/2 - ...
                                        eta2(eta_id)*(2*y_steady(4,n_states)+y_steady(4,n_states)*eta1(eta_id)*rand_num(obs_count,1))*rand_num(obs_count,2);
                                end
                                
                                if(I_dau(i)<0)
                                    I_dau(i)=0;
                                end
                                
                                %%% finding the new steady state of the
                                %%% daugther cellular content
                                y_steady_dau(:,i)=steady_points([y_steady(1:3,n_states);I_dau(i)]);
                                
                                %%% assigning phenotype to daughter cells
                                %%% based upon their cellular content
                                phenotype(obs_count,i+1)=get_phenotype(y_steady_dau(2,i),I_dau(i));
                            end
                            
                        end
                        
                        %%% checking which type of division has taken place
                        if(phenotype(obs_count,2)==phenotype(obs_count,1) && phenotype(obs_count,3)==phenotype(obs_count,1))
                            sym_count=sym_count+1;
                        elseif(~isempty(find(phenotype(obs_count,2:3)==phenotype(obs_count,1),1)))
                            asym_count=asym_count+1;
                        elseif(isempty(find(phenotype(obs_count,2:3)==phenotype(obs_count,1),1)))
                            div_count=div_count+1;
                        end
                        obs_count=obs_count+1;
                    end
                    
                    %%% determining fraction of each division type over all
                    %%% iterations
                    sym_count=sym_count/total_obs;
                    asym_count=asym_count/total_obs;
                    div_count=div_count/total_obs;
                    
                    
                    %%% storing data frame format so the plot in the paper can be generated  
                    if(get_phenotype(y_steady(2,n_states),y_steady(4,n_states)) == 0)
                        probability((indx-1)*9+1+(sim_num-1)*length(I)*9,1) = sym_count;
                        probability((indx-1)*9+2+(sim_num-1)*length(I)*9,1) = asym_count;
                        probability((indx-1)*9+3+(sim_num-1)*length(I)*9,1) = div_count;
                    elseif(get_phenotype(y_steady(2,n_states),y_steady(4,n_states)) == 1)
                        probability((indx-1)*9+4+(sim_num-1)*length(I)*9,1) = sym_count;
                        probability((indx-1)*9+5+(sim_num-1)*length(I)*9,1) = asym_count;
                        probability((indx-1)*9+6+(sim_num-1)*length(I)*9,1) = div_count;
                    else
                        probability((indx-1)*9+7+(sim_num-1)*length(I)*9,1) = sym_count;
                        probability((indx-1)*9+8+(sim_num-1)*length(I)*9,1) = asym_count;
                        probability((indx-1)*9+9+(sim_num-1)*length(I)*9,1) = div_count;
                    end
                    
                    probability((indx-1)*9+1+(sim_num-1)*length(I)*9,2) = 0;
                    probability((indx-1)*9+2+(sim_num-1)*length(I)*9,2) = 1;
                    probability((indx-1)*9+3+(sim_num-1)*length(I)*9,2) = 2;
                    
                    probability((indx-1)*9+4+(sim_num-1)*length(I)*9,2) = 0;
                    probability((indx-1)*9+5+(sim_num-1)*length(I)*9,2) = 1;
                    probability((indx-1)*9+6+(sim_num-1)*length(I)*9,2) = 2;
                    
                    probability((indx-1)*9+7+(sim_num-1)*length(I)*9,2) = 0;
                    probability((indx-1)*9+8+(sim_num-1)*length(I)*9,2) = 1;
                    probability((indx-1)*9+9+(sim_num-1)*length(I)*9,2) = 2;
                    
                    probability((indx-1)*9+1+(sim_num-1)*length(I)*9,3) = 0;
                    probability((indx-1)*9+2+(sim_num-1)*length(I)*9,3) = 0;
                    probability((indx-1)*9+3+(sim_num-1)*length(I)*9,3) = 0;
                    
                    probability((indx-1)*9+4+(sim_num-1)*length(I)*9,3) = 1;
                    probability((indx-1)*9+5+(sim_num-1)*length(I)*9,3) = 1;
                    probability((indx-1)*9+6+(sim_num-1)*length(I)*9,3) = 1;
                    
                    probability((indx-1)*9+7+(sim_num-1)*length(I)*9,3) = 2;
                    probability((indx-1)*9+8+(sim_num-1)*length(I)*9,3) = 2;
                    probability((indx-1)*9+9+(sim_num-1)*length(I)*9,3) = 2;
                    
                    probability((indx-1)*9+1+(sim_num-1)*length(I)*9,4) = I(indx);
                    probability((indx-1)*9+2+(sim_num-1)*length(I)*9,4) = I(indx);
                    probability((indx-1)*9+3+(sim_num-1)*length(I)*9,4) = I(indx);
                    
                    probability((indx-1)*9+4+(sim_num-1)*length(I)*9,4) = I(indx);
                    probability((indx-1)*9+5+(sim_num-1)*length(I)*9,4) = I(indx);
                    probability((indx-1)*9+6+(sim_num-1)*length(I)*9,4) = I(indx);
                    
                    probability((indx-1)*9+7+(sim_num-1)*length(I)*9,4) = I(indx);
                    probability((indx-1)*9+8+(sim_num-1)*length(I)*9,4) = I(indx);
                    probability((indx-1)*9+9+(sim_num-1)*length(I)*9,4) = I(indx);
                    
                    probability((indx-1)*9+1+(sim_num-1)*length(I)*9,5) = sim_num;
                    probability((indx-1)*9+2+(sim_num-1)*length(I)*9,5) = sim_num;
                    probability((indx-1)*9+3+(sim_num-1)*length(I)*9,5) = sim_num;
                    
                    probability((indx-1)*9+4+(sim_num-1)*length(I)*9,5) = sim_num;
                    probability((indx-1)*9+5+(sim_num-1)*length(I)*9,5) = sim_num;
                    probability((indx-1)*9+6+(sim_num-1)*length(I)*9,5) = sim_num;
                    
                    probability((indx-1)*9+7+(sim_num-1)*length(I)*9,5) = sim_num;
                    probability((indx-1)*9+8+(sim_num-1)*length(I)*9,5) = sim_num;
                    probability((indx-1)*9+9+(sim_num-1)*length(I)*9,5) = sim_num;
                    
                    probability((indx-1)*9+1+(sim_num-1)*length(I)*9,6) = noise_pert;
                    probability((indx-1)*9+2+(sim_num-1)*length(I)*9,6) = noise_pert;
                    probability((indx-1)*9+3+(sim_num-1)*length(I)*9,6) = noise_pert;
                    
                    probability((indx-1)*9+4+(sim_num-1)*length(I)*9,6) = noise_pert;
                    probability((indx-1)*9+5+(sim_num-1)*length(I)*9,6) = noise_pert;
                    probability((indx-1)*9+6+(sim_num-1)*length(I)*9,6) = noise_pert;
                    
                    probability((indx-1)*9+7+(sim_num-1)*length(I)*9,6) = noise_pert;
                    probability((indx-1)*9+8+(sim_num-1)*length(I)*9,6) = noise_pert;
                    probability((indx-1)*9+9+(sim_num-1)*length(I)*9,6) = noise_pert;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
            end
        end
        
        div_cat = categorical(probability(:,2),[0 1 2],{'Sym','Asym','Div'});
        phen_cat = categorical(probability(:,3),[0 1 2],{'E','H','M'});
        noise_cat = categorical(probability(:,6),[1 2],{'SNAIL','ALL'});
        
        varnames = ["Probability";"Division_type"; "Phenotype";"SNAIL_level";"Sim_run";"Noise"];
        prob_data = table(probability(:,1),div_cat,phen_cat,probability(:,4),probability(:,5),noise_cat,'VariableNames',varnames);
        
        
        %%% to update the output file if it was already present
        cd('/home/csb/CollaborationWork-Paras/Population_Dynamics/Final figures, data and codes/Dependent_formalism/Division type probabilities')
        if(isfile(['Multiple_runs_Division_type_prob_data_eta1_' num2str(eta1(eta_id)) '_eta2_' num2str(eta2(eta_id)) '.csv']))
            data = readtable(['Multiple_runs_Division_type_prob_data_eta1_' num2str(eta1(eta_id)) '_eta2_' num2str(eta2(eta_id)) '.csv']);
            prob_data.Sim_run =  prob_data.Sim_run + height(data)/(length(I)*3*3);
            prob_data = [data;prob_data];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         cd('/home/csb/CollaborationWork-Paras/Population_Dynamics/Final figures, data and codes/Dependent_formalism/Division type probabilities')
        if(noise_pert == 1)
            writetable(prob_data,['Multiple_runs_Division_type_prob_data_eta1_' num2str(eta1(eta_id)) '_eta2_' num2str(eta2(eta_id)) '.csv']);
        elseif(noise_pert == 2)
            writetable(prob_data,['Noise_in_all_Division_type_prob_data_eta1_' num2str(eta1(eta_id)) '_eta2_' num2str(eta2(eta_id)) '.csv']);
        end
end