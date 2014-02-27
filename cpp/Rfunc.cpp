#include "Sim.h"	
#include "Tree.h"

extern "C" void genTree(const int *ec,  const int *nc, double *a, int *start, int *end, double *len, double *sigma, double *dt, double *tip, double *tip2)
{
	Tree phy;
	
	phy.setValues (*ec, *nc, start, end, len, tip, tip2);
	
	phy.simulation(*a, *sigma, *dt);

	for(int i = 0; i < *ec; ++i)
	{
		tip[i] = phy.value[i];
		tip2[i]= phy.value2[i];
	} 

}
