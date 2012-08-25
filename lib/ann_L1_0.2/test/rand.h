//----------------------------------------------------------------------
//	File:		rand.h
//	Programmer:	Sunil Arya and David Mount
//	Last modified:	03/04/98 (Release 0.1)
//	Description:	Basic include file for random point generators
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

#ifndef rand_H
#define rand_H

//----------------------------------------------------------------------
//  Basic includes
//----------------------------------------------------------------------
#include <math.h>			// math routines
#include <ANN/ANN.h>			// basic ANN includes

//----------------------------------------------------------------------
//  Globals
//----------------------------------------------------------------------
extern	int	idum;			// used for random number generation

//----------------------------------------------------------------------
//  External entry points
//----------------------------------------------------------------------

double ran0();			// uniform random number generator

void uniform_pts(		// uniform distribution
	ANNpointArray	pa,		// point array (modified)
	int		n,		// number of points
	int		dim);		// dimension

void gauss_pts(			// Gaussian distribution
	ANNpointArray	pa,		// point array (modified)
	int		n,		// number of points
	int		dim,		// dimension
	double		std_dev);	// standard deviation

void co_gauss_pts(		// correlated-Gaussian distribution
	ANNpointArray	pa,		// point array (modified)
	int		n,		// number of points
	int		dim,		// dimension
	double		correlation);	// correlation

void laplace_pts(		// Laplacian distribution
	ANNpointArray	pa,		// point array (modified)
	int		n,		// number of points
	int		dim);		// dimension

void co_laplace_pts(		// correlated-Laplacian distribution
	ANNpointArray	pa,		// point array (modified)
	int		n,		// number of points
	int		dim,		// dimension
	double		correlation);	// correlation

void clus_gauss_pts(		// clustered-Gaussian distribution
	ANNpointArray	pa,		// point array (modified)
	int		n,		// number of points
	int		dim,		// dimension
	int		n_col,		// number of colors (clusters)
	ANNbool		new_clust,	// generate new cluster centers
	double		std_dev);	// standard deviation within clusters

void par_lines_pts(		// uniform on two parallel lines
	ANNpointArray	pa,		// point array (modified)
	int		n,		// number of points
	int		dim);		// dimension

#endif
