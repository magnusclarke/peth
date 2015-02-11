#include "Sim.h"	
#include "Tree.h"

extern "C" void genTree(int *ec, int *nc, int *nt, int *kernel, int *ratecut, double *a, int *start, int *end, double *len, double *sigma, double *sigma2, double *dt, double *lim, double *tip, double *cov)
{
	Tree phy;
	
	phy.setValues (*ec, *nc, *nt, start, end, len, tip, cov);

	phy.simulation(*a, *sigma, *sigma2, *dt, *lim, *kernel, *ratecut);

	for (int j = 0; j < *nt; ++j)
	{
		for(int i = 0; i < *ec; ++i)
		{
			tip[i + j*(*ec)] = phy.vals[j][i];
		} 
	}
}
