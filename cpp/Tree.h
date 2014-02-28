#ifndef __TREE_H__
#define __TREE_H__

#include "boost/random.hpp"	
#include "Sim.h"

using std::vector;

class Tree 
{
	int			edge_count, tip_count, root;
	std::vector<int>	st, end;
	std::vector<double>	length;

	void	simSeg(double&, double&, double&, vector<int>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&);

public:
	std::vector<double>	value;
	std::vector<double>	value2;
	
	void 	setValues (int&, int&, int[], int[], double[], double[], double[]);
	void 	simulation (double&, double&, double&);
};

#endif