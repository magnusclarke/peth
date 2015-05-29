#ifndef __TREE_H__
#define __TREE_H__

#include "Sim.h"

using std::vector;

class Tree 
{
	int		edge_count, tip_count, Ntraits, root;
	std::vector<int>	st, end;
	std::vector<double>	length;

	void simSeg(double &a, double &s, double &s2, double &dt, double &lim, int &kernel, int &ratecut, vector<int> &node, vector<vector<double> > &nodeVal, vector<vector<double> > &run_vals);
	void simNF(double&, double&, double&, double&, int&, vector<int>&, vector<vector<double> >&, vector<vector<double> >&);

public:
	std::vector<std::vector<double> > vals;

	void 	setValues (int&, int&, int&, int[], int[], double[], double[]);
	void 	simulation(double &a, double &sigma, double &sigma2, double &dt, double &lim, int &kernel, int &ratecut);

};

#endif