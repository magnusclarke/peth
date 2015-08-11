#include "Sim.h"	

/* RNG Seed: number cycles since processor on (assumed unique per core ) */
unsigned long long rdtsc()
{
	unsigned int lo,hi;
	__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
        return ((unsigned long long)hi << 32) | lo;
}

void Sim::BMsim(std::vector<std::vector<double> > &run_vals, int *nt, int *l, double *a, double *s, int *count, double *dt)
{
	len = *l;
	Ntraits = *nt;	
	double s2 = (*s) * (*s);
	double s2_time = *count * s2;

	// random numbers and random normal distribution
	boost::mt19937_64 generator; 			
	generator.seed(rdtsc());  
	boost::normal_distribution<double> distribution;

	for(int i=0; i < len; ++i)
	{
		for (int j = 0; j < Ntraits; ++j)
		{
			run_vals[j][i] += s2_time * distribution(generator);
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
		// Adding this ( for a!=0 ) roughly doubles runtime
		if(*a != 0)		denSim(run_vals, a, s, dt);
	}
}

void Sim::LIMsim(std::vector<std::vector<double> > &run_vals, int *nt, int *l, double *a, double *s, int *count, double *dt, double *lim)
{
	len = *l;
	Ntraits = *nt;	

	//double lim = 2;				// max size of deviation from zero.

	boost::mt19937_64 generator; 			
	generator.seed(rdtsc());  
	boost::normal_distribution<double> distribution;

	for(int j = 0; j < *count; j++) 
	{
		for(int i=0; i < len; i++)
		{
			for (int k = 0; k < Ntraits; ++k)
			{
				run_vals[k][i] += (*s) * distribution(generator);
				if(run_vals[k][i] > *lim)	run_vals[k][i] = *lim;
				if(run_vals[k][i] < -*lim)	run_vals[k][i] = -*lim;
			}
		} 		
		if(*a != 0)		denSim(run_vals, a, s, dt);
	}
}

void Sim::denSim(vector<vector<double> > &run_vals, double *a, double *s, double *dt)
{
	std::vector<double> dists (Ntraits, 0);
	std::vector<double> sqDists (Ntraits, 0);
	double sumDist 		= 0;
	double sumSqDist 	= 0;
	double sqrtDT_s 	= 0.5 * sqrt(*dt) / (*s);
	double s_sqrtDT 	= (*s) / sqrt(*dt);
	double dta 			= (*dt) * (*a);

	for (int i=0; i < len; ++i) 
	{
		for (int j = 0; j < len; ++j)
		{
			sumDist 	= 0;
			sumSqDist 	= 0;

			for (int k = 0; k < Ntraits; ++k)
			{
				dists[k]	= (run_vals[k][i] - run_vals[k][j]);
				sqDists[k]	= dists[k]*dists[k];
				sumDist		+= std::abs(dists[k]);
				sumSqDist	+= sqDists[k];
			}

			// The first factor scales (INVERSELY) the width of the resource-use curves; 2nd is distance
			double q	= sqrtDT_s * sqrt(sumSqDist);

			// simple approximation to pnorm(q, 0, 1, FALSE, FALSE); good to 2dp
			double pn;
			if(q <= 2.2)					{ pn = 0.1 * q * (4.4 - q); }
			else if (q > 2.2 && q < 2.6)	{ pn = 0.49; }
			else if (q > 2.6)				{ pn = 0.50; }
			else							{ pn = 0.50; }

            //double g = 0.5 * dta * (0.5-pn);
			double g = 1 * dta * (0.5-pn);
			
			// Change in trait values, scaled by sigma
			if(sumDist != 0)
			{
				for (int k = 0; k < Ntraits; ++k)
				{
					run_vals[k][i] += s_sqrtDT*g*(dists[k]/sumDist);
					run_vals[k][j] -= s_sqrtDT*g*(dists[k]/sumDist);
				}
			}
		}
	}
}
