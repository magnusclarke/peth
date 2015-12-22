#include <vector>
#include <cmath>
#include <random>
#include "sim.h"
#include "tree.h"

using std::vector;

/* 	RNG. Declaring this inside the looped function would cost the earth in time. */
std::random_device rd;
std::mt19937 gen(rd());
std::normal_distribution<> d(0,1);	// concievably faster to gen all rnums at outset, then use.
std::uniform_real_distribution<> ud(0, 1);	// random number from uniform dist between 0 and 1.

/* 	Modify trait values for segment between speciation
	events. Needs number of time steps within that segment. 	
	This is where all the runtime is spent!	 */
void Sim::evolve_segment(int &nsteps)
{
	for (int step = 0; step < nsteps; ++step)
	{
		step_segment();
	}
}

/* 	Update trait values for one time step. 	*/
void Sim::step_segment()
{
	num_species = tval.size();

	for (int species = 0; species < num_species; ++species)
	{
		step_species(species);
	}
}

/* 	Do one evolutionary step on one species  */
void Sim::step_species(int &species)
{
	int Ntraits = 1;	// one trait for now

	tval[species] += d(gen) * rate;	// BM

	std::vector<double> dists (Ntraits, 0);
	std::vector<double> sqDists (Ntraits, 0);

	// Loop over species i that are interacting with this species.
	for (int i = 0; i < num_species; ++i)
	{	

		sumDist 	= 0;
		sumSqDist 	= 0;
		
		// Loop over traits. N is just 1 for now
		for (int k = 0; k < Ntraits; ++k)		
		{
			dists[k]	= (tval[species] - tval[i]);
			sqDists[k]	= dists[k]*dists[k];
			sumDist		+= std::abs(dists[k]);
			sumSqDist	+= sqDists[k];
		}

		// The first factor scales (INVERSELY) the width of the resource-use curves; 2nd is distance
		double q	= 0.5 * sqrt(sumSqDist);

		// simple approximation to pnorm(q, 0, 1, FALSE, FALSE); good to 2dp
		double pn;
		if(q <= 2.2)					{ pn = 0.1 * q * (4.4 - q); }
		else if (q > 2.2 && q < 2.6)	{ pn = 0.49; }
		else if (q > 2.6)				{ pn = 0.50; }
		else							{ pn = 0.50; }

		double g = 2 * dt * a * (0.5-pn);	// factor of 2, so overlap varies from 0 to 1
		
		// Change in trait values, scaled by rate
		if(sumDist != 0)
		{
			for (int k = 0; k < Ntraits; ++k)
			{
				tval[species] += rate*g*(dists[k]/sumDist);	// Would also have k index if Ntraits > 1.
				tval[i] -= rate*g*(dists[k]/sumDist);
			}
		}
	}
}

void Sim::path()
{
	for (int i = 0; i < num_segment; ++i)
	{
		int x = tree.speciators[i];
		tval.push_back(tval[x]);
		evolve_segment(segment_steps[i]);
	}
}

// Doesn't use up runtime!!
void Sim::set_values(double &r_dt, double &r_rate,double &r_a, double r_intervals[], Tree &t)
{
	tree = t;
	num_segment = tree.num_tips-1;

	dt = r_dt;
	rate = r_rate * sqrt(dt);			// Peth adds a sqrt(dt) factor to sigma. Not sure why.
	a = r_a;
	total_time_steps = tree.total_time / dt;

	// Compute number of time steps for each segment
	segment_steps.assign(num_segment, 0);
	for (int i = 0; i < num_segment; ++i)
	{
		segment_steps[i] = r_intervals[i] / dt;
	}

	// Root trait value is in middle of fitness surface; set all to this
	double root = 0;
	vector<double> root_species;
	root_species.assign(1, root);		// one trait
	tval.assign(1, root);		// initially one species

}
