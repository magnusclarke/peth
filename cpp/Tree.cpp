#include "Tree.h"

void Tree::setValues(int x,int y, int start[], int en[], double len[], double tip[], double tip2[])
{
	edge_count 	= x;
	tip_count 	= y;
	root		= start[0];
	value.assign	(x, 0);
	value2.assign	(x, 0);
	st.assign	(x, 0);
	end.assign	(x, 0);
	length.assign	(x, 0);

	for(int i = 0; i < x; i++) 
	{
		value[i]	= tip[i];
		value2[i]	= tip2[i];
		st[i]		= start[i];
		end[i]		= en[i];
		length[i]	= len[i];
	}
}

void Tree::simulation(double a, double sigma, double dt)
{
	double s = sigma * sqrt(dt);
	std::vector<double> 	run_val (1,0);		// values of running lines
	std::vector<double> 	run_val2 (1,0);		// values of running lines
	std::vector<int> 	node;			// nodes due to begin running
	std::vector<double>	val;			// value at each of 'node'
	std::vector<double>	val2;			// value at each of 'node'

	simSeg(a, s, dt, node, val, val2, run_val, run_val2);
}

void Tree::simSeg(double a, double s, double dt, std::vector<int> node, std::vector<double> val, std::vector<double> val2, std::vector<double> run_val, std::vector<double> run_val2)
{

	bool run[edge_count];
	for(int i = 0; i < edge_count; i++) 
	{
		if(st[i] == root)	run[i] = true;
 		else 			run[i] = false;
	}

	while(run_val.empty() == false){

		run_val.clear();
		run_val2.clear();

		// set run_val and time
		double time = 1e308;
		for(int i=0; i<edge_count; i++){
			if (run[i] == true) {
				run_val.push_back( value[i] );
				run_val2.push_back( value2[i] );
				if (length[i] < time)	time = length[i];
			}
		}

		// runSim on lines where run=TRUE
		int l = run_val.size();
		int count = ( time / dt );
		Sim segment;
		segment.setSeg(run_val.data(), run_val2.data(), &l);
		segment.runSim( &a, &s, &count, &dt);
		segment.getSeg(run_val.data(), run_val2.data(), &l);

		int z = 0;
		for(int i=0; i < edge_count; i++){
			if (run[i] == true)
			{
				value[i] = run_val[z];
				value2[i] = run_val2[z];
				length[i] -= time;
				z++;
			}
		}

		for (int i = 0; i < edge_count; i++)
		{
			if (length[i] <= 0 && run[i] == true) {
				node.push_back( end[i] );
				val.push_back( value[i] );
				val2.push_back( value2[i] );
			}
		}

		for(int i=0; i < edge_count; i++)
		{
			for(unsigned int j=0; j < node.size(); j++){
				if(st[i] == node[j]) {
					value[i] = val[j];
					value2[i] = val2[j];
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
		val.clear();
		val2.clear();
	}
}
