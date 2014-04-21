#include "Sim.h"	

/* RNG Seed: number cycles since processor on (assumed unique per core ) */
unsigned long long rdtsc()
{
	unsigned int lo,hi;
	__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
        return ((unsigned long long)hi << 32) | lo;
}

void Sim::denSim(vector<vector<double> > run_vals, double *a, double *s, double *dt)
{
	for (int i=0; i < len; ++i) 
	{
		for (int j = 0; j < len; ++j)
		{
			// Distance between trait values, scaled by 1/sigma
			std::vector<double> dists (Ntraits, 0);
			std::vector<double> sqDists (Ntraits, 0);
			double sumDist = 0;
			double sumSqDist = 0;
			for (int k = 0; k < Ntraits; ++k)
			{
				dists[k]	= (run_vals[k][i] - run_vals[k][j]);
				sqDists[k]	= dists[k]*dists[k];
				sumDist		+= std::abs(dists[k]);
				sumSqDist	+= sqDists[k];
			}

			double q	= (sqrt(*dt) / (*s) ) * 0.5 * sqrt(sumSqDist);

			// simple approximation to pnorm(q, 0, 1, FALSE, FALSE); good to 2dp
			double pn;
			if(q <= 2.2)					{ pn = 0.1 * q * (4.4 - q); }
			else if (q > 2.2 && q < 2.6)	{ pn = 0.49; }
			else if (q > 2.6)				{ pn = 0.50; }
			else							{ pn = 0.50; }

			double g = (*dt) * (*a) * (0.5-pn);
			
			// Change in trait values, scaled by sigma
			if(sumDist != 0)
			{
				for (int k = 0; k < Ntraits; ++k)
				{
					run_vals[k][i] += ((*s)/sqrt(*dt))*g*(dists[i]/sumDist);
					run_vals[k][j] -= ((*s)/sqrt(*dt))*g*(dists[i]/sumDist);
				}
			}
		}
	}
}

void Sim::CEsim(std::vector<std::vector<double> > &run_vals, int *nt, int *l, double *a, double *s, int *count, double *dt)
{
	len = *l;

	Ntraits = *nt;		

	// random numbers and random normal distribution
	boost::mt19937_64 generator; 			
	generator.seed(rdtsc());  
	boost::normal_distribution<double> distribution;

	// Time loop
	for(int j = 0; j < *count; j++) 
	{
		// BM evolution, here uncorrelated between traits for simplicity
		// Loop over number of species
		for(int i=0; i < len; i++)
		{
			for (int k = 0; k < Ntraits; ++k)
			{
				run_vals[k][i] += (*s) * distribution(generator);
			}
		} 		
		// Density-dependent evolution	
		if(*a != 0)		denSim(run_vals, a, s, dt);
	}
}

void Sim::BMsim(std::vector<std::vector<double> > &run_vals, int *nt, int *l, double *a, double *s, int *count, double *dt)
{

}