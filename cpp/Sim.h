#ifndef __RUNSIM_H__
#define __RUNSIM_H__

#include <vector>
#include "boost/random/mersenne_twister.hpp"
#include "boost/random/normal_distribution.hpp"

using std::vector;

class Sim
{			
	int len, Ntraits;
	void denSim(vector<vector<double> > run_vals, double *a, double *s, double *dt);
public:
	void CEsim(vector<vector<double> > &run_vals, int *nt, int *l, double *a, double *s, int *count, double *dt);
	void BMsim(vector<vector<double> > &run_vals, int *nt, int *l, double *a, double *s, int *count, double *dt);
};

#endif