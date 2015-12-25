#include <vector>
#include <cmath>
#include <random>
#include "sim.h"
#include "tree.h"

using std::vector;

/* 	Random Number Generator  */
std::random_device rd;
std::mt19937 gen(rd());
std::normal_distribution<> d(0,1);

/*  Simulate dataset  */
void Sim::path()
{
	for (int i = 0; i < num_segment; ++i)
	{
		int x = tree.speciators[i];
		tval.push_back(tval[x]);
		evolve_segment(segment_steps[i]);
	}
}

/* 	Modify trait values for one segment (i.e. between speciation
	events). Needs parameter: number of time steps within that segment. */
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

/* 	Do one evolutionary step on one species.  */
void Sim::step_species(int &species)
{
	tval[species] += d(gen) * rate;		// BM
	// Loop over species i that are interacting with this species.
	for (int i = 0; i < num_species; ++i)
	{	
		interaction(species, i);
	}
}

/* 	Interaction between two species; updates trait values.  */
void Sim::interaction(int s1, int s2)
{
	update_distance(s1, s2);
	
	/* Compute resource distribution overlap g (0 to 1) between two species. */
	double q = 0.5 * sqrt(sumSqDist);
	double pn = pnorm(q);
	double g = 2 * dt * a * (0.5-pn); 
	
	/* Update trait values. */
	if (sumDist != 0)
	{
		for (int k = 0; k < num_traits; ++k)
		{
			tval[s1] += rate*g*(dists[k]/sumDist);	// Would also have k index if Ntraits > 1.
			tval[s2] -= rate*g*(dists[k]/sumDist);
		}
	}
}

void Sim::update_distance(int s1, int s2)
{
	sumDist = 0;
	sumSqDist = 0;

	/* Loop over traits and get trait distances between species s1 and s2. */
	for (int trait = 0; trait < num_traits; ++trait)		
	{
		dists[trait]	= tval[s1] - tval[s2];
		sqDists[trait]	= dists[trait]*dists[trait];
		sumDist		+= std::abs(dists[trait]);
		sumSqDist	+= sqDists[trait];
	}
}

/*  Simple approximation to R's pnorm(q, 0, 1, FALSE, FALSE); good to 2dp.  */
double Sim::pnorm(double q)
{
	double pn;
	if(q <= 2.2)
	{ 
		pn = 0.1 * q * (4.4 - q); 
	} else if (q > 2.2 && q < 2.6){ 
		pn = 0.49; 
	} else if (q > 2.6){ 
		pn = 0.50; 
	}
	else{ 
		pn = 0.50; 
	}
	return(pn);
}

/*  Copy all the parameters from R  */
void Sim::set_values(double &r_dt, double &r_rate,double &r_a, double r_intervals[], Tree &t)
{
	tree = t;
	num_traits = 1;						// This needs to become an int parameter.
	num_segment = tree.num_tips-1;

	dt = r_dt;
	a = r_a;
	rate = r_rate * sqrt(dt);	
	total_time_steps = tree.total_time / dt;

	/* Compute number of time steps for each segment */
	segment_steps.assign(num_segment, 0);
	for (int i = 0; i < num_segment; ++i)
	{
		segment_steps[i] = r_intervals[i] / dt;
	}

	/* Root trait value is zero. We start with a single trait vector for the root. */
	double root = 0;
	vector<double> root_species;
	root_species.assign(num_traits, root);		// one trait
	tval.assign(1, root);		
	
	dists.assign (num_traits, 0);
	sqDists.assign (num_traits, 0);
}
