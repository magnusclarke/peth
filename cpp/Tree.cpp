#include "Tree.h"
using std::vector;

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

void Tree::simulation(double &a, double &sigma, double &dt, int &kernel)
{
	double s = sigma * sqrt(dt);

	vector<vector<double> > run_vals;	// values of running lines
	vector<double> 	single_value (1,0);
	run_vals.assign(Ntraits, single_value);			

	vector<int> 	node;				// nodes due to begin running

	vector<vector<double> > nodeVal;
	vector<double> single_nodeVal;
	nodeVal.assign(Ntraits, single_nodeVal);

	simSeg(a, s, dt, kernel, node, nodeVal, run_vals);
}

void Tree::simSeg(double &a, double &s, double &dt, int &kernel, vector<int> &node, vector<vector<double> > &nodeVal, vector<vector<double> > &run_vals)
{

	bool run[edge_count];
	for(int i = 0; i < edge_count; i++) 
	{
		if(st[i] == root)	run[i] = true;
 		else 				run[i] = false;
	}

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

		Sim segment;
		if(kernel==0){
			segment.BMsim(run_vals, &Ntraits, &l, &a, &s, &count, &dt);
		} else if(kernel==1) {
			segment.CEsim(run_vals, &Ntraits, &l, &a, &s, &count, &dt);
		}

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
