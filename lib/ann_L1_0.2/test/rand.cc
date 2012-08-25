//----------------------------------------------------------------------
//	File:		rand.cc
//	Programmer:	Sunil Arya and David Mount
//	Last modified:	03/04/98 (Release 0.1)
//	Description:	Routines for random point generation
//----------------------------------------------------------------------
// Copyright (c) 1997-1998 University of Maryland and Sunil Arya and David
// Mount.  All Rights Reserved.
// 
// This software and related documentation is part of the 
// Approximate Nearest Neighbor Library (ANN).
// 
// Permission to use, copy, and distribute this software and its 
// documentation is hereby granted free of charge, provided that 
// (1) it is not a component of a commercial product, and 
// (2) this notice appears in all copies of the software and
//     related documentation. 
// 
// The University of Maryland (U.M.) and the authors make no representations
// about the suitability or fitness of this software for any purpose.  It is
// provided "as is" without express or implied warranty.
//----------------------------------------------------------------------

#include "rand.h"			// random generator declarations

//----------------------------------------------------------------------
//  Globals
//----------------------------------------------------------------------
int	idum = 0;			// used for random number generation

//------------------------------------------------------------------------
//	ran0 - (safer) uniform random number generator
//
//	The code given here is taken from "Numerical Recipes in C" by
//	William Press, Brian Flannery, Saul Teukolsky, and William
//	Vetterling. The task of the code is to do an additional randomizing
//	shuffle on the system-supplied random number generator to make it
//	safer to use. 
//
//	Returns a uniform deviate between 0.0 and 1.0 using the
//	system-supplied routine "random()". Set idum to any negative value
//	to initialise or reinitialise the sequence.
//------------------------------------------------------------------------

           				// "standard" random number generators
#ifdef USING_SGI			// SGI versions
extern "C" long random(void);
extern "C" void srandom(unsigned);
#endif
#ifdef USING_SUNOS5			// SUNOS5 versions
extern "C" int random(void);
extern "C" int srandom(unsigned);
#endif

double ran0()
{
    int j;

    static double y, maxran, v[98];	// The exact number 98 is unimportant
    static int iff = 0;

    // As a precaution against misuse, we will always initialize on the first
    // call, even if "idum" is not set negative. Determine "maxran", the next
    // integer after the largest representable value of type int. We assume 
    // this is a factor of 2 smaller than the corresponding value of type
    // unsigned int. 

    if (idum < 0 || iff == 0) {	// initialize
	unsigned i, k;
	iff = 1;
	i = 2;
	do {
	    k = i;
	    i <<= 1;
	} while (i);
	maxran = (double) k;
   
	srandom(idum);
	idum = 1;

	for (j = 1; j <= 97; j++)	// exercise the system routine
	    random();			// (value intentionally ignored)

	for (j = 1; j <= 97; j++)	// Then save 97 values and a 98th
	    v[j] = random();
	y = random();
     }

    // This is where we start if not initializing. Use the previously saved
    // random number y to get an index j between 1 and 97. Then use the
    // corresponding v[j] for both the next j and as the output number. */

    j = 1 + (int) (97.0 * (y / maxran));
    y = v[j];
    v[j] = random();			// Finally, refill the table entry
					// with the next random number from
					// "random()" 
    return(y / maxran);
}

//------------------------------------------------------------------------
//  ran_gauss - Gaussian random number generator
//	Returns a normally distributed deviate with zero mean and unit
//	variance, using ran0() as the source of uniform deviates.
//------------------------------------------------------------------------

double ran_gauss()
{
    static int iset=0;
    static double gset;

    if (iset == 0) {			// we don't have a deviate handy
	double v1, v2;
	double r = 2.0;
	while (r >= 1.0) {
	    //------------------------------------------------------------
	    // Pick two uniform numbers in the square extending from -1 to
	    // +1 in each direction, see if they are in the circle of radius
	    // 1.  If not, try again 
	    //------------------------------------------------------------
	    v1 = 2.0 * ran0() - 1.0; 
	    v2 = 2.0 * ran0() - 1.0;
	    r = v1 * v1 + v2 * v2;
	}
        double fac = sqrt(-2.0 * log(r) / r);
	//-----------------------------------------------------------------
	// Now make the Box-Muller transformation to get two normal
	// deviates.  Return one and save the other for next time.
	//-----------------------------------------------------------------
	gset = v1 * fac;
	iset = 1;		    	// set flag
	return v2 * fac;
    }
    else {				// we have an extra deviate handy
	iset = 0;			// so unset the flag
	return gset;			// and return it
    }
}

//------------------------------------------------------------------------
//  ran_laplace - Laplacian random number generator
//	Returns a laplacian distributed deviate with zero mean and
//	unit variance, using ran0() as the source of uniform deviates. 
//
//		prob(x) = b/2 * exp(-b * |x|).
//
//	b is chosen to be sqrt(2.0) so that the variance of the laplacian
//	distribution [2/(b^2)] becomes 1. 
//------------------------------------------------------------------------

double ran_laplace() 
{
    const double b = 1.4142136;

    double laprand = -log(ran0()) / b;
    double sign = ran0();
    if (sign < 0.5) laprand = -laprand;
    return(laprand);
}

//----------------------------------------------------------------------
//  uniform_pts - Generate uniformly distributed points
//	A uniform distribution over [-1,1].
//----------------------------------------------------------------------

void uniform_pts(		// uniform distribution
	ANNpointArray	pa,		// point array (modified)
	int		n,		// number of points
	int		dim)		// dimension
{
    for (int i = 0; i < n; i++) {
	for (int d = 0; d < dim; d++) {
	    pa[i][d] = (ANNcoord) (2*ran0() - 1);
	}
    }
}

//----------------------------------------------------------------------
//  gauss_pts - Generate Gaussian distributed points
//	A Gaussian distribution with zero mean and the given standard
//	deviation.
//----------------------------------------------------------------------

void gauss_pts(			// Gaussian distribution
	ANNpointArray	pa,		// point array (modified)
	int		n,		// number of points
	int		dim,		// dimension
	double		std_dev)	// standard deviation
{
    for (int i = 0; i < n; i++) {
	for (int d = 0; d < dim; d++) {
	    pa[i][d] = (ANNcoord) (ran_gauss() * std_dev);
	}
    }
}

//----------------------------------------------------------------------
//  laplace_pts - Generate Laplacian distributed points
//	Generates a Laplacian distribution (zero mean and unit variance).
//----------------------------------------------------------------------

void laplace_pts(		// Laplacian distribution
	ANNpointArray	pa,		// point array (modified)
	int		n,		// number of points
	int		dim)		// dimension
{
    for (int i = 0; i < n; i++) {
	for (int d = 0; d < dim; d++) {
            pa[i][d] = (ANNcoord) ran_laplace();
	}
    }
}

//----------------------------------------------------------------------
//  co_gauss_pts - Generate correlated Gaussian distributed points
//	Generates a Gauss-Markov distribution of zero mean and unit
//	variance.
//----------------------------------------------------------------------

void co_gauss_pts(		// correlated-Gaussian distribution
	ANNpointArray	pa,		// point array (modified)
	int		n,		// number of points
	int		dim,		// dimension
	double		correlation)	// correlation
{
    double std_dev_w = sqrt(1.0 - correlation * correlation);
    for (int i = 0; i < n; i++) {
	double previous = ran_gauss();
	pa[i][0] = (ANNcoord) previous;
	for (int d = 1; d < dim; d++) {
	    previous = correlation*previous + std_dev_w*ran_gauss();
	    pa[i][d] = (ANNcoord) previous;
	} 
    }
}

//----------------------------------------------------------------------
//  co_laplace_pts - Generate correlated Laplacian distributed points
//	Generates a Laplacian-Markov distribution of zero mean and unit
//	variance.
//----------------------------------------------------------------------

void co_laplace_pts(		// correlated-Laplacian distribution
	ANNpointArray	pa,		// point array (modified)
	int		n,		// number of points
	int		dim,		// dimension
	double		correlation)	// correlation
{
    double wn;
    double corr_sq = correlation * correlation;

    for (int i = 0; i < n; i++) {
	double previous = ran_laplace();
	pa[i][0] = (ANNcoord) previous;
	for (int d = 1; d < dim; d++) {
	    double temp = ran0();
	    if (temp < corr_sq)
		wn = 0.0;
	    else
		wn = ran_laplace();
	    previous = correlation * previous + wn;
	    pa[i][d] = (ANNcoord) previous;
        } 
    }
}

//----------------------------------------------------------------------
//  clus_gauss_pts - Generate clusters of Gaussian distributed points
//	Cluster centers are uniformly distributed over [0,1], and the
//	standard deviation within each cluster is fixed.
//
//	Note: Once cluster centers have been set, they are not changed,
//	unless new_clust = true.  This is so that subsequent calls generate
//	points from the same distribution.  It follows, of course, that any
//	attempt to change the dimension or number of clusters without
//	generating new clusters is asking for trouble.
//
//	Note: Cluster centers are not generated by a call to uniform_pts().
//	Although this could be done, it has been omitted for
//	compatibility with clus_gauss_pts() in the colored version,
//	rand_c.cc.
//----------------------------------------------------------------------

void clus_gauss_pts(		// clustered-Gaussian distribution
	ANNpointArray	pa,		// point array (modified)
	int		n,		// number of points
	int		dim,		// dimension
	int		n_col,		// number of colors
	ANNbool		new_clust,	// generate new clusters.
	double		std_dev)	// standard deviation within clusters
{
    static ANNpointArray clusters = NULL;// cluster storage

    if (clusters == NULL || new_clust) {// need new cluster centers
	if (clusters != NULL)		// clusters already exist
	    annDeallocPts(clusters);	// get rid of them
	clusters = annAllocPts(n_col, dim);
					// generate cluster center coords
	for (int i = 0; i < n_col; i++) {
	    for (int d = 0; d < dim; d++) {
		clusters[i][d] = (ANNcoord) ran0();
	    }
	}
    }

    for (int i = 0; i < n; i++) {
	int c = (int) (ran0() * n_col);// generate cluster index
	for (int d = 0; d < dim; d++) {
          pa[i][d] = (ANNcoord) (std_dev*ran_gauss() + clusters[c][d]);
	}
    }
}

//----------------------------------------------------------------------
//  par_lines_pts - points uniform over 2 parallel lines
//	The points are generated along two parallel line segments.
//	In all cases the first coordinate is uniformly distributed
//	over the interval [-1,1].  Half of the remaining coordinates
//	are very close to 1 and the the other half are very close to
//	-1.
//
//	This distribution is nasty for a number of splitting rules,
//	and the optimized kd-tree in particular.  There is essentially
//	only one choice for the orientation of the cutting plane, and
//	so a large number of parallel cuts are generated, resulting in
//	the limit in a series of "wafer-like" cells.
//
//	To make this a bad scenario at query time, it is also important
//	that query points be selected from a different distribution from
//	the data points.  A good choice is a uniform or gaussian 
//	distribution.
//
//----------------------------------------------------------------------

void par_lines_pts(			// uniform on 2 parallel line segments
	ANNpointArray	pa,		// point array (modified)
	int		n,		// number of points
	int		dim)		// dimension
{
    int i;
    const double EPS = 1.0e-6;		// a very small value

    int mid = n/2;

    for (i = 0; i < mid; i++) {		// generate first half
					// first coord is uniform [-1,1]
	pa[i][0] = (ANNcoord) (2*ran0() - 1);
	for (int d = 1; d < dim; d++) {	// remainder close to -1
	    pa[i][d] = (ANNcoord) (-1 + ran_gauss()*EPS);
	}
    }
    for (i = mid; i < n; i++) {		// generate latter half
					// first coord is uniform [-1,1]
	pa[i][0] = (ANNcoord) (2*ran0() - 1);
	for (int d = 1; d < dim; d++) {	// remainder close to 1
	    pa[i][d] = (ANNcoord) (1 + ran_gauss()*EPS);
	}
    }
}
