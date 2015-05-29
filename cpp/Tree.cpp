#include "Tree.h"
#include <vector>
#include <random>
#include <cmath>

using std::vector;

/* 	RNG. Seeded from platform's random device.
	These have to be labelled differently than those in Sim.cpp
	to keep g++ happy. Not sure why.	  */
std::random_device rd2;		
std::mt19937 generator2(rd2());	
std::normal_distribution<> distribution2(0,1);

void Tree::setValues(int &x, int &y, int &nt, int start[], int en[], double len[], double tip[])
{
	edge_count 	= x;
	tip_count 	= y;
	Ntraits 	= nt;
	root		= start[0];
	st.assign	(x, 0);
	end.assign	(x, 0);
	length.assign	(x, 0);

	std::vector<double> vec_start;
	vec_start.assign (x, 0);	

	vals.assign (Ntraits, vec_start);	

	for(int i = 0; i < x; i++) 
	{
		st[i]		= start[i];
		end[i]		= en[i];
		length[i]	= len[i];
	}
}

void Tree::simulation(double &a, double &sigma, double &sigma2, double &dt, double &lim, int &kernel, int &ratecut)
{
	double s 	= sigma * sqrt(dt);
	double s2 	= sigma2 * sqrt(dt);

	vector<vector<double> > run_vals;	// values of running lines
	vector<double> 	single_value (1,0);
	run_vals.assign(Ntraits, single_value);			

	vector<int> 	node;				// nodes due to begin running

	vector<vector<double> > nodeVal;
	vector<double> single_nodeVal;
	nodeVal.assign(Ntraits, single_nodeVal);

	// NF simulation
	if(kernel==3){
		simNF(a, s, dt, lim, kernel, node, nodeVal, run_vals);		
	} else {
		simSeg(a, s, s2, dt, lim, kernel, ratecut, node, nodeVal, run_vals);
	} 
}

void Tree::simNF(double &a, double &s, double &dt, double &lim, int &kernel, vector<int> &node, vector<vector<double> > &nodeVal, vector<vector<double> > &run_vals)
{
	/* variance of niche distribution. The parameter sigma controls this for the NF model. */
	double nicheVar = s/sqrt(dt);		
 
	bool run[edge_count];
	//int x=0;
	for(int i = 0; i < edge_count; i++) 
	{
		if(st[i] == root){
			run[i] = true;
			// new vals get random value from niche distribution (1 trait only here)
			vals[0][i] += nicheVar * distribution2(generator2);
		} else {
			run[i] = false;			
		} 	
	}

	// While there are branches to simulate
	while(run_vals[0].empty() == false){

		// run_vals is the working copy of species being simulated
		// clear it ready for a new simulation.
		for (int i = 0; i < Ntraits; ++i)
		{
			run_vals[i].clear();
		}

		// 'vals' are the trait values on the branches.
		// Copy the vals of running branches to run_vals
		// & Set 'time' to the time until the next node. 
		double time = 1e308;
		for(int i=0; i<edge_count; i++){
			if (run[i] == true) {
				for (int j = 0; j < Ntraits; ++j)
				{
					run_vals[j].push_back( vals[j][i] );
				}
				if (length[i] < time)	time = length[i];
			}
		}

		// Copy the run_vals back to the corresponding vals.
		// Reduce the running branch lengths by 'time'
		int z = 0;
		for(int i=0; i < edge_count; i++){
			if (run[i] == true)
			{
				for (int j = 0; j < Ntraits; ++j)
				{
					vals[j][i] = run_vals[j][z];
				}
				length[i] -= time;
				z++;
			}
		}

		// If a branch has ended, add its end-node to vector 'node'
		// & add its vals value to nodeVal
		for (int i = 0; i < edge_count; i++)
		{
			if (length[i] <= 0 && run[i] == true) {
				node.push_back( end[i] );
				for (int j = 0; j < Ntraits; ++j)
				{
					nodeVal[j].push_back( vals[j][i] );
				}
			}
		}

		// Branches that start at a node in 'node' are due to begin simulation.
		// Copy nodeVals to the corresponding new branches.
		// And set those branches running.
		// And make 1 new branch jump to a random niche, which is closest to itself.
		int x = 0;
		for(int i=0; i < edge_count; i++)
		{
			for(unsigned int j=0; j < node.size(); j++)
			{
				if(st[i] == node[j]) 
				{
					for (int k = 0; k < Ntraits; ++k)
					{
						vals[k][i] = nodeVal[k][j];
						// alternate new vals get random value from niche distribution
						if(x%2==0)
						{
							// Generate niches until niche is closest to vals[k][i] out of all vals[k]
							double newNiche;
							int use = 1;
							while(use==1){
								use = 0;
								newNiche	= nicheVar * distribution2(generator2);
								double difference = abs(newNiche - vals[k][i]);
								for(int m=0; m < edge_count; m++)
								{
									if(abs(newNiche - vals[k][i]) < (difference-0.00001))
									{
										use = 1;
									}
								}
							}
							vals[k][i] = newNiche;
						}	
						x = x+1;	
					}
					run[i] = true;
				}
			}
		}

		// Stop branches running when they reach zero length.
		for (int i = 0; i < edge_count; i++)
		{
			if (length[i] <= 0){
				run[i] = false;
			}
		}

		// Clear the 'working' list of evolving nodes and values.
		node.clear(); 
		for (int i = 0; i < Ntraits; ++i)
		{
			nodeVal[i].clear();
		}
	}
	
}


void Tree::simSeg(double &a, double &s, double &s2, double &dt, double &lim, int &kernel, int &ratecut, vector<int> &node, vector<vector<double> > &nodeVal, vector<vector<double> > &run_vals)
{

	bool run[edge_count];
	for(int i = 0; i < edge_count; i++) 
	{
		if(st[i] == root)	run[i] = true;
 		else 				run[i] = false;
	}

	int node_counter = 0;

	while(run_vals[0].empty() == false){

		for (int i = 0; i < Ntraits; ++i)
		{
			run_vals[i].clear();
		}

		// set run_val and time
		double time = 1e308;
		for(int i=0; i<edge_count; i++){
			if (run[i] == true) {
				for (int j = 0; j < Ntraits; ++j)
				{
					run_vals[j].push_back( vals[j][i] );
				}
				if (length[i] < time)	time = length[i];
			}
		}

		// // runSim on lines where run=TRUE
		int l = run_vals[0].size();
		int count = ( time / dt );

		// Use BM rate s up to the ratecut node, then use the rate s2.
		double rate = 0;
		if(node_counter < ratecut){
			rate = s;
		} else if(node_counter >= ratecut) {
			rate = s2;
		}

		Sim segment;
		if(kernel==0 || kernel==4){
			segment.BMsim(run_vals, &Ntraits, &l, &a, &rate, &count, &dt);
		} else if(kernel==1) {
			segment.CEsim(run_vals, &Ntraits, &l, &a, &rate, &count, &dt);
		} else if(kernel==2) {
			segment.LIMsim(run_vals, &Ntraits, &l, &a, &rate, &count, &dt, &lim);
		}
		node_counter++;

		int z = 0;
		for(int i=0; i < edge_count; i++){
			if (run[i] == true)
			{
				for (int j = 0; j < Ntraits; ++j)
				{
					vals[j][i] = run_vals[j][z];
				}
				length[i] -= time;
				z++;
			}
		}

		for (int i = 0; i < edge_count; i++)
		{
			if (length[i] <= 0 && run[i] == true) {
				node.push_back( end[i] );
				for (int j = 0; j < Ntraits; ++j)
				{
					nodeVal[j].push_back( vals[j][i] );
				}
			}
		}

		for(int i=0; i < edge_count; i++)
		{
			for(unsigned int j=0; j < node.size(); j++){
				if(st[i] == node[j]) {
					for (int k = 0; k < Ntraits; ++k)
					{
						vals[k][i] = nodeVal[k][j];
					}
					run[i] = true;
				}
			}
		}

		for (int i = 0; i < edge_count; i++)
		{
			if (length[i] <= 0){
				run[i] = false;
			}
		}

		node.clear(); 
		for (int i = 0; i < Ntraits; ++i)
		{
			nodeVal[i].clear();
		}
	}
}
