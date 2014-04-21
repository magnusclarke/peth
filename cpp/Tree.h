#ifndef __TREE_H__
#define __TREE_H__

#include "Sim.h"

using std::vector;

class Tree 
{
	int		edge_count, tip_count, Ntraits, root;
	std::vector<int>	st, end;
	std::vector<double>	length;

	void	simSeg(double&, double&, double&, vector<int>&, vector<vector<double> >&, vector<vector<double> >&);

public:
	std::vector<std::vector<double> > vals;

	void 	setValues (int&, int&, int&, int[], int[], double[], double[]);
	void 	simulation (double&, double&, double&);
};

#endif
