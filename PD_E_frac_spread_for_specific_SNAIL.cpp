////////////////////////////////////////////////////////////////////////////////////////////////
////                     Population Dynamics single cell simulations                       /////
////                          (For predtermined SNAIL levels)				               /////
////                Codes adopted from Tripathi et. al. 2020 PLOS Comp Bio                 /////
////                Modification done by Jain et al., Date: 19th Jan 2022                  /////
//////////////////////////////////////////////////////////////////////////////////////////////// 

#include <iostream>
#include <cmath>
#include <ctime>
#include <random>
#include <string>
#include <fstream>
#include <boost/array.hpp>
#include <boost/unordered_map.hpp>
#include <boost/numeric/odeint.hpp>

#define NUMNODES 4

typedef boost::array<double, NUMNODES> cell_state;
typedef boost::unordered_map<int, boost::array<double, NUMNODES>> population;

std::mt19937 generator; // for random number generation

double temp_var[13];  // temp_var hold the value for initialization parameters provided at start 

//////////////////////// k combination out of n //////////////////////////////
double nchoosek(int n, int k) 
{
	double result = 1.0;
	double result0 = 1.0;

	for(int i = 0; i < k; i++)
	{
		result *= (n - i);
	}
	for(int i = 1; i <= k; i++)
	{
		result0 *= i;
	}

	return(result / result0);
}
///////////////////////////////////////////////////////////////////////////////////

// function to break a string and separate numeral character based on the delimiter
void tokenize(std::string s, std::string del = " ")
{
	int start = 0;
	int end = s.find(del);
	int i = 0;
	while (end != -1) {
		//std::cout << s.substr(start, end - start) << std::endl;
		temp_var[i] = std::stod(s.substr(start, end - start));
		//std::cout << temp_var[i] << std::endl;
		start = end + del.size();
		end = s.find(del, start);
		i++;
	}
	//std::cout << s.substr(start, end - start);
	temp_var[i] = std::stod(s.substr(start, end - start));
	//std::cout << temp_var[i];
}
///////////////////////////////////////////////////////////////////////////////////

/////////////////// EMT core circuit ODEs (updates a cells' state) ////////////////
void EMT_system(const cell_state &x, cell_state &dxdt, double t)
{
	double ku200 = 0.05, kmz = 0.5, kz = 0.1;
	 
	// Transcription rate:
	double gu200 = 2100, gmz = 11, gz = 100;

	// Hills function threshold :
	double z0u200 = 220000, z0mz = 25000, s0u200 = 180000, s0mz = 180000, u2000 = 10000;

	// Cooperativity:
	double nzu200 = 3, nsu200 = 2, nzmz = 2, nsmz = 2, nu200 = 6;

	// fold change
	double lamdazu200 = 0.1, lamdasu200 = 0.1, lamdazmz = 7.5, lamdasmz = 10;



	double Mu0=1/std::pow((1+x[0]/u2000),nu200);
	double Mu1=std::pow((x[0]/u2000),1)/std::pow((1+x[0]/u2000),nu200);
	double Mu2=std::pow((x[0]/u2000),2)/std::pow((1+x[0]/u2000),nu200);
	double Mu3=std::pow((x[0]/u2000),3)/std::pow((1+x[0]/u2000),nu200);
	double Mu4=std::pow((x[0]/u2000),4)/std::pow((1+x[0]/u2000),nu200);
	double Mu5=std::pow((x[0]/u2000),5)/std::pow((1+x[0]/u2000),nu200);
	double Mu6=std::pow((x[0]/u2000),6)/std::pow((1+x[0]/u2000),nu200);


		
	double Hillszu200=(1+lamdazu200*std::pow((x[2]/z0u200),nzu200))/(1+std::pow((x[2]/z0u200),nzu200));
	double Hillssu200=(1+lamdasu200*std::pow((x[3]/s0u200),nsu200))/(1+std::pow((x[3]/s0u200),nsu200));
	double Hillszmz=(1+lamdazmz*std::pow((x[2]/z0mz),nzmz))/(1+std::pow((x[2]/z0mz),nzmz));
	double Hillssmz=(1+lamdasmz*std::pow((x[3]/s0mz),nsmz))/(1+std::pow((x[3]/s0mz),nsmz));
	 

	dxdt[0] = gu200*Hillszu200*Hillssu200-x[1]*(0.005*6*Mu1+2*0.05*15*Mu2+3*0.5*20*Mu3+4*0.5*15*Mu4+5*0.5*6*Mu5+6*0.5*Mu6)-ku200*x[0];
	dxdt[1] = gmz*Hillszmz*Hillssmz-x[1]*(0.04*6*Mu1+0.2*15*Mu2+20*Mu3+15*Mu4+6*Mu5+Mu6)-kmz*x[1];
	dxdt[2] = gz*x[1]*(Mu0+0.6*6*Mu1+0.3*15*Mu2+0.1*20*Mu3+0.05*15*Mu4+0.05*6*Mu5+0.05*Mu6)-kz*x[2];
	dxdt[3] = 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////// to assign phenotypes to cells /////////////////////////
void get_phenotypes(population &P, boost::unordered_map <int, int> &phenotype)
{
	double x1 = 193.2, x2 = 208.7, y1 = 243.09, y2 = 90.93;
	double a1 = y2 - y1, b1 = x1 - x2, c1 = x2*y1 - y2*x1;

	double u1 = 185.12, u2 = 224.67, v1 = 698.395, v2 = 495.802;
	double a2 = v2 - v1, b2 = u1 - u2, c2 = u2*v1 - v2*u1;

	double x, y, fac1, fac2;
	int state = -1;

	for(int i = 0; i < P.size(); i++)
	{
		x = P[i][3] / 1e3;
		y = P[i][1];
		state = -1;

		if(x < u1)
		{
			state = 0;
		}
		else if(x > u2)
		{
			state = 2;
		}

		else
		{
			fac1 = (a1*x + b1*y + c1) / b1;
			fac2 = (a2*x + b2*y + c2) / b2;
			if(fac2 >= 0)
			{
				state = 2;
			}
			else if(fac1 < 0)
			{
				state = 0;
			}
			else
				state = 1;
			
		}

		if(state == -1)
		{
			std::cout << "Error in phenotype assignment." << "\n";
		}
		else
		{
			phenotype[i] = state;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////

//// To initialize cellular states of the initial population
void initialize_Signal_lognormal(population& P, double* pop_fraction, const int pop_size, const double SNAIL)
{

	//std::normal_distribution <> dist{ 0.0, 1.0 };

	

	//std::vector<int> V[3];

	int flag = 0;
	//const unsigned int N = 100000;
	double I = SNAIL;
	
	/*double CV2 = 1.0;
	double SD = std::sqrt(std::log(CV2 * CV2 + 1.0));
	double M = 200e3;
	*/


	boost::unordered_map <int, int> phenotype;
	boost::array<int, 4> count = { 0,0,0,0 };

	population temp_P, temp_P_vector;
	std::ifstream steady_states, eigen_values;

	steady_states.open("steady_states.csv");
	eigen_values.open("eigen_values.csv");
	std::string line;

// load steady states matrix (save steady_states.csv in working directory)////

	int num_states = 0;

	while(getline(steady_states,line)){
		
		//std::cout<<line;
		tokenize(line,",");

		for(int i = 0; i < NUMNODES; i++){
			temp_P[num_states][i] = temp_var[i]; 
		}
			num_states++;
	}


// to load eigen values matrix (save steady_states.csv in working directory) ////

	num_states = 0;

	while(getline(eigen_values,line)){
		
		//std::cout<<line;
		tokenize(line,",");

		for(int i = 0; i < NUMNODES; i++){
			temp_P_vector[num_states][i] = temp_var[i]; 
		}
			num_states++;
	}
	///////////////////////////////////////////////////////////////

	steady_states.close();
	eigen_values.close();

	//////// to initialize cells based on random generated SNAIL levels and population fraction.///
	///// The logic is to search for nearest SNAIL value, compared to sampled I, in steady states and eigen values matrices 
	/// and initialize cells with miR200, mZEB, and ZEB variable vector at that SNAIL value  

	int first_time = 1, min_indx = 0, selected_P_indx = 0 ;
	double abs_diff = 0, diff = 0;

	while (count[3] < (pop_size))
	{

		if(I <= temp_P[temp_P.size()-1][3]){
			first_time = 1;
			min_indx = 0;
			selected_P_indx = 0 ;
			abs_diff = 0;
			diff = 0;

			population selected_P;

			for(int g = 0; g < temp_P.size(); g++){
				
				if(temp_P[g][3] - I >= 0){
					abs_diff = temp_P[g][3] - I;
				}
				else
				{
					abs_diff = -(temp_P[g][3] - I);	
				}

				if(abs_diff<= 100 && temp_P_vector[g][0] < 0 && temp_P_vector[g][1] < 0 && temp_P_vector[g][2] < 0){

					if(first_time==1){
						diff = abs_diff;
						min_indx = g;
						first_time = 2;
					}
					else{
						if(abs_diff < diff)
							min_indx = g;
					}

				}
				else if(first_time == 2)
				{

				for(int i = 0; i < NUMNODES; i++){
					selected_P[selected_P_indx][i] = temp_P[g][i];
					
					/*std::cout<< selected_P[selected_P_indx][i] << " ";

					if(i ==  NUMNODES-1)
						std::cout<<"\n";*/				

				}

				selected_P_indx++;
				first_time = 1;

				}
			}


			get_phenotypes(selected_P,phenotype);
			/*std::cout << selected_P.size() << "\n";

			for(int i = 0; i < phenotype.size(); i++){
				std::cout<<phenotype[i] << "\n";
			}*/ 

			for(int i = 0; i < selected_P.size();i++){

				if (phenotype[i] == 0 && count[0] < pop_size * pop_fraction[0]) {	
					for(int j = 0; j < NUMNODES; j++){
						P[count[3]][j] = selected_P[i][j];
						
						/*std::cout<< selected_P[i][j] << " ";

						if(j ==  NUMNODES-1)
							std::cout<<"\n";*/
					}
					count[0]++;
					

					
				}
				else if (phenotype[i] == 1 && count[1] < pop_size * pop_fraction[1]){
					for(int j = 0; j < NUMNODES; j++){
						P[count[3]][j] = selected_P[i][j];

						/*std::cout<< selected_P[i][j] << " ";

						if(j ==  NUMNODES-1)
							std::cout<<"\n";*/

					}
					count[1]++;
					
					
				}
				else if (phenotype[i] == 2 && count[2] < pop_size * pop_fraction[2]){
					for(int j = 0; j < NUMNODES; j++){
						P[count[3]][j] = selected_P[i][j];
						/*std::cout<< selected_P[i][j] << " ";

						if(j ==  NUMNODES-1)
							std::cout<<"\n";*/
					}
					count[2]++;
					
					
				}

				count[3] = count[0] + count[1] + count[2];

			}
			

			selected_P.clear();
			phenotype.clear();
		}
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	
}


///// uniform sampling a cell from a subpopulation////////////////
int get_sample(boost::unordered_map <int, int> &phenotype, int pid)
{
	std::vector<int> V;
	int count = 0;
	for(int i = 0; i < phenotype.size(); i++)
	{
		if(phenotype[i] == pid)
		{
			V.push_back(i);
			count += 1;
		}
	}

	std::uniform_int_distribution<int> distribution(0, count - 1);

	return(V[distribution(generator)]);
}
//////////////////////////////////////////////////////////////////////


////////////////// adding duplication and partitioning noise/fluctuation on cell division ////////// 
void add_sym_noise_normal(population &P, int rep_id, int new_id, double * eta, int eta_id)
{
	std::normal_distribution <> distribution{0, 1};

	std::uniform_real_distribution <> distribution_unifrom{0.2, 0.4};
	
	double temp_snail[2];

	double rand_num[2];

	rand_num[0] =  distribution(generator);
	rand_num[1] =  distribution(generator);

	//do{
	if(eta_id == 0){
	        temp_snail[0] = P[rep_id][3] + distribution(generator)*distribution_unifrom(generator)*P[rep_id][3];
	        temp_snail[1] = P[rep_id][3] + distribution(generator)*distribution_unifrom(generator)*P[rep_id][3];
        }
        else{
        	temp_snail[0] = P[rep_id][3] + rand_num[0]*eta[0]*P[rep_id][3]/2 + rand_num[1]*eta[1]*(2*P[rep_id][3] + rand_num[0]*eta[0]*P[rep_id][3]) ;
	        temp_snail[1] = P[rep_id][3] + rand_num[0]*eta[0]*P[rep_id][3]/2 - rand_num[1]*eta[1]*(2*P[rep_id][3] + rand_num[0]*eta[0]*P[rep_id][3]) ;
        }

        //std::cout << temp_snail[0] << " " << temp_snail[1];

	//}while(temp_snail[0] <= 0.0 || temp_snail[1] <=0.0);

	P[rep_id][3] = temp_snail[0];
	P[new_id][3] = temp_snail[1];

        if(P[rep_id][3] < 0.0)
        {
                P[rep_id][3] = 0.0;
        }
        if(P[new_id][3] < 0.0)
        {
                P[new_id][3] = 0.0;
        }
        
}

/////////////////////////////////////////////////////////////////////////////////


////////// sampling a fraction of cells from a subpopulation when overall population size increase 80% of carrying capacity K
int FACS(population &P, int N, int type_start, std::vector <int> &cell_index)
{
	boost::unordered_map <int, int> phenotype;
	get_phenotypes(P, phenotype);
	std::vector<int> usefulIndex;
	int count = 0, flag = 1;

	for(int i = 0; i < phenotype.size(); i++)
	{
		if(phenotype[i] == type_start)
		{
			usefulIndex.push_back(i);
			count += 1;
		}
	}
	if(count < N)
	{
		std::cout << "Not enough cells of type" << std::to_string(type_start) << " present in the culture." << "\n";
		flag = 0;
		return(flag);
	}

	std::uniform_int_distribution<int> distribution(0, usefulIndex.size() - 1);

	int size = 0;
	int index = 0;
	while(size < N)
	{
		//index = usefulIndex[distribution(generator)];
		
		index = usefulIndex[size];
		cell_index.push_back(index);

		/*sorted_P[size+cell_count] = {0.0, 0.0, 0.0, 0.0};
		for(int i = 0; i < NUMNODES; i++)
		{
			sorted_P[size+cell_count][i] = P[index][i];
		}*/
		size += 1;
	}

	return(flag);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////// population dynamics ////////////////////////////////////////////////////////////
int simulate_normal(population& P, double end_time, double * eta, int eta_id,int frac_id, int sim_type, boost::array <double, 3> GR, int GR_ratio_id, int pop_size, int sim_num, int SNAIL_id)
{
	boost::unordered_map <int, double> update_time;
	boost::unordered_map <int, int> phenotype;
	boost::array<double, 3> count;
	//boost::array<int, 3> cell_count;
	std::vector <int> cell_index;

	double temp_count[3];
	double day_indx = 1;

	for(int i = 0; i < P.size(); i++)
	{
		update_time[i] = 0.0;
	}

	int last_index = P.size();

	get_phenotypes(P, phenotype);

	//std::cout << GR[0] << " " << GR[1] << " " << GR[2] << "\n"; 

	boost::array <double, 3> growth;
	growth[0] = std::log(2) / GR[0];
	growth[1] = std::log(2) / GR[1];
	growth[2] = std::log(2) / GR[2];
	boost::array<double, 3> r0 = {growth[0], growth[1], growth[2]};
	boost::array<double, 3> d0 = {growth[0] / 10.0, growth[1] / 10.0, growth[2] / 10.0};
	boost::array<double, 6> rates = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};


	int K = 5000;

	if(sim_type == 1)
	{
		K = 2000;
	}

	for(int i = 0; i < P.size(); i++)
	{
		if(phenotype[i] == 0)
		{
			if(P.size() < K)
			{
				rates[0] += r0[0]*(1 - float(P.size()) / K);
			}
			rates[3] += d0[0];
		}
		else if(phenotype[i] == 1)
		{
			if(P.size() < K)
			{
				rates[1] += r0[1]*(1 - float(P.size()) / K);
			}
			rates[4] += d0[1];
		}
		else
		{
			if(P.size() < K)
			{
				rates[2] += r0[2]*(1 - float(P.size()) / K);
			}
			rates[5] += d0[2];
		}
	}

	std::uniform_real_distribution<double> distribution(0.0, 1.0);

	double t = 0.0, a0, dt, p0, p1;
	double sum_prev, sum_new;
	double last_updated = 0.0;
	int event_id, rep_id, death_id;


	std::ofstream end_of_day, SNAIL_distribution, mZEB_distribution;

	if(sim_type == 1){

		count = {0.0, 0.0, 0.0};
		for(int i = 0; i < P.size(); i++)
		{
			count[phenotype[i]] += 1.0;
		}
		for(int i = 0; i < 3; i++)
		{
			count[i] = count[i] / P.size();
		}

		
		// Below, 'E_frac_spread_Data_20hrs_DT_4_weeks' is the folder name in which output data files will be stored.
		// Create it in current working directory before running the simulation. 
		end_of_day.open("E_frac_spread_Data_20hrs_DT_4_weeks/end_of_day_counts_" + std::to_string(eta_id) + "_" + std::to_string(frac_id) + "_" + std::to_string(GR_ratio_id) + "_" + std::to_string(SNAIL_id) + ".csv", std::ios::app);
		
		//SNAIL_distribution.open("E_frac_spread_SNAIL_dist_20hrs_DT_4_weeks/snail_dist_" + std::to_string(eta_id) + "_" + std::to_string(frac_id) + "_" + std::to_string(GR_ratio_id) + "_" + std::to_string(SNAIL_id) + "_" + std::to_string(sim_num) + ".csv");

		//mZEB_distribution.open("E_frac_spread_ZEB_dist_20hrs_DT_4_weeks/mZEB_dist_" + std::to_string(eta_id) + "_" + std::to_string(frac_id) + "_" + std::to_string(GR_ratio_id) + "_" + std::to_string(SNAIL_id) + "_" + std::to_string(sim_num) + ".csv");
		
		/*for(int i = 0; i < P.size(); i++){
			if(i < P.size() -1){
				SNAIL_distribution << P[i][3] << ",";
				mZEB_distribution << P[i][1] << ",";
			}
			else{
				SNAIL_distribution << P[i][3] << "\n";
				mZEB_distribution << P[i][1] << "\n";
			}	
		}*/
		//phenotype_switch.open("phenotype_switch/end_of_day_counts_" + std::to_string(eta_id) + "_" + std::to_string(frac_id) + "_" + std::to_string(GR_ratio_id) + "_" + std::to_string(sim_num) + ".csv", std::ios::app);
	}
	

	while(t <= end_time)
	{
		a0 = 0.0;

		for(int i = 0; i < 6; i++)
		{
			a0 += rates[i];
		}

		p0 = distribution(generator);
		dt = (1.0 / a0)*std::log(1 / p0);
		t += dt;

		p1 = distribution(generator);
		sum_prev = 0.0;
		sum_new = 0.0;
		event_id = -1;

		for(int i = 0; i < 6; i++)
		{
			sum_new = sum_prev + rates[i] / a0;
			if(p1 >= sum_prev && p1 < sum_new)
			{
				event_id = i;
				break;
			}
			sum_prev = sum_new;
		}

		if(event_id == 0)
		{
			rep_id = get_sample(phenotype, 0);
		}
		else if(event_id == 1)
		{
			rep_id = get_sample(phenotype, 1);
		}
		else if(event_id == 2)
		{
			rep_id = get_sample(phenotype, 2);
		}
		else if(event_id == 3)
		{
			death_id = get_sample(phenotype, 0);
		}
		else if(event_id == 4)
		{
			death_id = get_sample(phenotype, 1);
		}
		else if(event_id == 5)
		{
			death_id = get_sample(phenotype, 2);
		}

		if(event_id < 3)
		{	
			
			boost::numeric::odeint::integrate(EMT_system, P[rep_id], 0.0, t - update_time[rep_id], 0.1);
			
			P[last_index] = {0.0, 0.0, 0.0, 0.0};  // last index holds the length of Population Array vector 
			for(int j = 0; j < NUMNODES; j++)
			{
				P[last_index][j] = P[rep_id][j];
			}

			add_sym_noise_normal(P, rep_id, last_index, eta, eta_id);

			update_time[rep_id] = t;
			update_time[last_index] = t;
			last_index += 1;
		}
		else
		{
			for(int j = 0; j < NUMNODES; j++)
			{
				P[death_id][j] = P[last_index - 1][j];
			}
			update_time[death_id] = update_time[last_index - 1];
			phenotype[death_id] = phenotype[last_index - 1];
			P.erase(last_index - 1);
			update_time.erase(last_index - 1);
			phenotype.erase(last_index - 1);
			last_index -= 1;
		}

		if(P.size()==0)
		{
			if(sim_type == 1)
			{
				end_of_day.close();
				//SNAIL_distribution.close();
				//mZEB_distribution.close();
			}
		
			return(1);
		}

		if(true)
		{
			for(int i = 0; i < P.size(); i++)
			{
				if(t - update_time[i] > 0.0)
				{
					boost::numeric::odeint::integrate(EMT_system, P[i], 0.0, t - update_time[i], 0.1);
				}
				update_time[i] = t;
			}
		}

		get_phenotypes(P, phenotype);


		/// reducing population size by sampling different subpopulations as per their last distribution

		if(P.size() > 0.8 * K && sim_type == 1){

			count = {0.0, 0.0, 0.0};

			for(int i = 0; i < P.size(); i++)
			{
				count[phenotype[i]] += 1.0;
			}

			for(int i = 0; i < 3; i++)
			{
				count[i] = count[i] / P.size();

				double b = 0; // to truncate the double to two decimal places

				for (double k = 0; k<=100; k++)
				{
					if(k > 100*count[i])
					{
						b = 100*count[i] - (k-1);
						count[i] = count[i] - b/100; 
						break;
					}
				}

				//std::cout << count[i] << "\n";

				int flag = FACS(P,pop_size*count[i],i,cell_index);
				if (flag == 0){
					std::cout <<" stopped simulation due to non-availability of cells\n";
					return(-1);
				}
			}
			
			std::sort(cell_index.begin(),cell_index.end());

			//for (int i = 0; i< cell_index.size(); i++) std::cout << cell_index[i] << "\n"; 
			//std::cout<< cell_index.size()<<"\n";

			for (int i = 0; i < cell_index.size(); i++){
					if(i != cell_index[i]){
						for(int j = 0; j < NUMNODES; j++){
							P[i][j] = P[cell_index[i]][j];
							//std::cout << P[i][j] << " ";
						}
						phenotype[i] = phenotype[cell_index[i]];
					}
					update_time[i] = t; 

					//std::cout<<"\n";
			}
			for (int i = P.size()-1; i >= cell_index.size(); i--)
			{
				P.erase(i);
				update_time.erase(i);
				phenotype.erase(i);
			}
			last_index = P.size();

			/*std::cout << "\nsize of the sampled pop is " << last_index <<"\n";
			
			for (int i = 0; i < P.size(); i++)
			{
				std::cout << P[i][1] << " " << P[i][3] << "\n";
			}*/
			

			cell_index.clear();

		}
		///////////////////////////////////////////////////////////////////////////////////////

		
		count = {0.0, 0.0, 0.0};

		for(int i = 0; i < P.size(); i++)
		{
			count[phenotype[i]] += 1.0;
		}

		/*for(int i = 0; i < 3; i++)
		{
			count[i] = count[i] / P.size();
		}*/

		
		rates = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

		for(int i = 0; i < P.size(); i++)
		{
			if(phenotype[i] == 0)
			{
				if(P.size() < K)
				{
					rates[0] += r0[0]*(1 - float(P.size()) / K);
				}
				rates[3] += d0[0];
			}
			else if(phenotype[i] == 1)
			{
				if(P.size() < K)
				{
					rates[1] += r0[1]*(1 - float(P.size()) / K);
				}
				rates[4] += d0[1];
			}
			else
			{
				if(P.size() < K)
				{
					rates[2] += r0[2]*(1 - float(P.size()) / K);
				}
				rates[5] += d0[2];
			}
		}

		if(sim_type == 1)
		{
			

			if (t < 24* day_indx) {
				temp_count[0] = count[0];
				temp_count[1] = count[1];
				temp_count[2] = count[2];
				
			}
			else {

				//std::cout << day_indx << "," << temp_count[0] << "," << temp_count[1] << "," << temp_count[2] << "," << P.size() << "\n";
			
				//// the phenotypic distribution is stored each day in absolute cell counts ////

				end_of_day << day_indx << "," << temp_count[0] << "," << temp_count[1] << "," << temp_count[2] << "," << (P.size()-1) << "\n";

				/*if((int(day_indx)%14) == 0){ // save distribution data every four weeks
					for(int i = 0; i < P.size(); i++){
						if(i < P.size() -1){
							SNAIL_distribution << P[i][3] << ",";
							mZEB_distribution << P[i][1] << ",";
						}
						else{
							SNAIL_distribution << P[i][3] << "\n";
							mZEB_distribution << P[i][1] << "\n";
						}	
					}
				}*/
				

				day_indx++;

				temp_count[0] = count[0];
				temp_count[1] = count[1];
				temp_count[2] = count[2];

			}
		}
	}

	if(sim_type == 1)
	{
		end_of_day.close();
	//	SNAIL_distribution.close();
	//	mZEB_distribution.close();
	}

	return(0);
}





int main(int argc, char* argv[])
{

	tokenize(argv[1], "_"); /// to separate parameters from input string 

	boost::unordered_map <int, int> phenotype;
	boost::array <double, 3> GR; // doubling time of the three population
	std::vector <int> cell_index; // used with FACS
	int aborted;
	
	int count[3];

	// order the input parameter values delimited by '_' as follows:
// e.g. 1_0.2_0.1_0_1_0_0_0_20_20_20_0_50000
	int eta_id = temp_var[0];
	double eta[2] = { temp_var[1], temp_var[2] };
	int frac_id = temp_var[3];
	double pop_fraction[3] = { temp_var[4], temp_var[5], temp_var[6] };
	int GR_ratio_id = temp_var[7];
	GR[0] = temp_var[8];
	GR[1] = temp_var[9];
	GR[2] = temp_var[10];
	int SNAIL_id = temp_var[11];
	double SNAIL = temp_var[12];

	int sim_num = 1;
	int pop_size = 200;
	int total_inde_runs = 50;
	int time_in_days = 4 * 7 ; // 4 respresent 4 weeks
	int file_rows = 0;

	std::ifstream csv_file;

	std::string line;

	while(true){
		
		//// check how many runs of data is stored in output file
		csv_file.open("E_frac_spread_Data_20hrs_DT_32_weeks/end_of_day_counts_" + std::to_string(eta_id) + "_" + std::to_string(frac_id) + "_" + 
			std::to_string(GR_ratio_id) + "_"  + std::to_string(SNAIL_id) + ".csv");

		/*if(csv_file.is_open())
		{
			std::cout<<"file is correctly opened\n";
		}*/
		
		file_rows = 0;		
		while(getline(csv_file,line)){
			file_rows++;
		}

		csv_file.close();
		//std::cout << file_rows << "\n";
		
		/////////// don't proceed if file has data of required number of runs or incomplete previous simulation data  
		if(file_rows >= (total_inde_runs * time_in_days) || (file_rows%time_in_days) != 0){
			
			if(file_rows%time_in_days != 0)
				std::cout << "stopped bacause of already missing data \n";
			else if(file_rows > (total_inde_runs * time_in_days))
				std::cout << total_inde_runs * time_in_days << " enough data already present\n";
			break;
		} 
		//////////////////////////////////////////////////////////////////////////////////////
		
		generator = std::mt19937(std::time(NULL) + file_rows);

		population P;

		initialize_Signal_lognormal(P, pop_fraction, 1, SNAIL); // initial pop size as 1 for single cell simualtion
		
		/*for(int i = 0; i < P.size(); i++){

			for(int j = 0; j < NUMNODES; j++){
				std::cout << P[i][j] << " ";
			}
			std::cout<<"\n";
		}*/

		
		std::cout << "initialization complete for eta id " << eta_id << " frac_id " << frac_id << " GR ratio id " << GR_ratio_id << " SNAIL_id " << SNAIL_id <<" and sim num is " << sim_num << std::endl;

		aborted = simulate_normal(P, 24 * time_in_days, eta, eta_id, frac_id, 1, GR, GR_ratio_id, pop_size, sim_num, SNAIL_id);

		if(aborted == 1)
		{
			std::cout<<"the simulation is aborted due to extinction of population\n";
		}else{

			std::cout << "simulation complete for eta id " << eta_id << " frac_id " << frac_id << " GR ratio id " << GR_ratio_id << " SNAIL_id " << SNAIL_id <<" and sim num is " << sim_num << std::endl;
			sim_num++;
		}

		P.clear(); 
	}
	


	return(0);
	
}