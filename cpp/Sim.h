#ifndef __RUNSIM_H__
#define __RUNSIM_H__

#include <vector>
#include "boost/random/mersenne_twister.hpp"
#include "boost/random/normal_distribution.hpp"

class Sim
{			
	int len;
	std::vector<double> seg;
	std::vector<double> seg2;

	void denSim(double *a, double *s, double *dt);
	
public:
	void setSeg(double *dat, double *dat2, int *l);
	void getSeg(double *dat, double *dat2, int *l);
	void runSim(double *a, double *s, int *count, double *dt);
};

#endif
