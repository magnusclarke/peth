#include "Sim.h"	
#include "Tree.h"

extern "C" void genTree(int *ec, int *nc, int *nt, int *kernel, double *a, int *start, int *end, double *len, double *sigma, double *dt, double *tip)
{
	Tree phy;
	
	phy.setValues (*ec, *nc, *nt, start, end, len, tip);

	phy.simulation(*a, *sigma, *dt, *kernel);

	for (int j = 0; j < *nt; ++j)
	{
		for(int i = 0; i < *ec; ++i)
		{
			tip[i + j*(*ec)] = phy.vals[j][i];
		} 
	}
}