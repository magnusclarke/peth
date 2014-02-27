#include "Sim.h"	

/* RNG Seed: number cycles since processor on (assumed unique per core ) */
unsigned long long rdtsc()
{
	unsigned int lo,hi;
	__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
        return ((unsigned long long)hi << 32) | lo;
}

void Sim::setSeg(double *dat, double *dat2, int *l)
{
	seg.assign(dat, dat+*l);
	seg2.assign(dat2, dat2+*l);
	len = *l;
}

void Sim::getSeg(double *dat, double *dat2, int *l)
{
	for (int i = 0; i < *l; i++) 
	{
		dat[i] = seg[i];
		dat2[i] = seg2[i];
	}
}

void Sim::denSim(double *a, double *s, double *dt)
{
	for (int i=0; i < len; ++i) 
	{
		for (int j = 0; j < len; ++j)
		{
			// Distance between trait values, scaled by 1/sigma
			double q1	= seg[i] - seg[j];
			double q2	= seg2[i] - seg2[j];
			double q12	= std::abs(q1) + std::abs(q2);
			double q 	= (sqrt(*dt) / (*s) ) * 0.5 * sqrt( (q1*q1) + (q2*q2) );

			// simple approximation to pnorm(q, 0, 1, FALSE, FALSE); good to 2dp
			double pn;
			if(q <= 2.2)				{ pn = 0.1 * q * (4.4 - q); }
			else if (q > 2.2 && q < 2.6)		{ pn = 0.49; }
			else if (q > 2.6)			{ pn = 0.50; }
			else					{ pn = 0.50; }

			double g = (*dt) * (*a) * (0.5-pn);
			
			// Change in trait values, scaled by sigma
			if(q12 != 0)
			{
				seg[i]	+= ( (*s) / sqrt(*dt) ) * g * (q1 / q12);
				seg2[i]	+= ( (*s) / sqrt(*dt) ) * g * (q2 / q12);
				seg[j]	-= ( (*s) / sqrt(*dt) ) * g * (q1 / q12);
				seg2[j]	-= ( (*s) / sqrt(*dt) ) * g * (q2 / q12);
			}
		}
	}
}

void Sim::runSim(double *a, double *s, int *count, double *dt) 
{
	// random numbers and random normal distribution
	boost::mt19937_64 generator; 			
	generator.seed(rdtsc());  
	boost::normal_distribution<double> distribution;

	// Time loop
	for(int j = 0; j < *count; j++) 
	{
		// BM evolution, here uncorrelated between traits for simplicity
		for(int i=0; i < len; i++)
		{
			seg[i]	+= (*s) * distribution(generator);  
			seg2[i]	+= (*s) * distribution(generator);  
		} 		
		// Density-dependent evolution
		if(*a != 0)		denSim(a, s, dt);
	}
}
