//----------------------------------------------------------------------
//	File:		brute.cc
//	Programmer:	Sunil Arya and David Mount
//	Last modified:	03/04/98 (Release 0.1)
//	Description:	Brute-force nearest neighbors
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

#include <ANN/ANNx.h>			// all ANN includes
#include "pr_queue_k.h"			// k element priority queue

//----------------------------------------------------------------------
//  	Brute-force search simply stores a pointer to the list of
//	data points and searches linearly for the nearest neighbor.
//	The k nearest neighbors are stored in a k-element priority
//	queue (which is implemented in a pretty dumb way as well).
//
//	If ANN_ALLOW_SELF_MATCH is ANNfalse then data points at distance
//	zero are not considered.
//----------------------------------------------------------------------

ANNbruteForce::ANNbruteForce(		// constructor from point array
    ANNpointArray	pa,		// point array
    int			n,		// number of points
    int			dd)		// dimension
{
    dim = dd;  n_pts = n;  pts = pa;
}

ANNbruteForce::~ANNbruteForce() { }	// destructor (empty)

void ANNbruteForce::annkSearch(		// approx k near neighbor search
    ANNpoint		q,		// query point
    int			k,		// number of near neighbors to return
    ANNidxArray		nn_idx,		// nearest neighbor indices (returned)
    ANNdistArray	dd,		// dist to near neighbors (returned)
    double		eps)		// error bound (ignored)
{
    ANNmin_k mk(k);			// construct a k-limited priority queue
    ANNdist  dist;			// distance to point
    int i;

    if (k > n_pts) {			// too many near neighbors?
	annError("Requesting more near neighbors than data points", ANNabort);
    }
					// run every point through queue
    for (i = 0; i < n_pts; i++) {
	dist = annDist(dim, pts[i], q);	// compute distance to point
	if (ANN_ALLOW_SELF_MATCH || dist != 0)
	    mk.insert(dist, i);
    }
    for (i = 0; i < k; i++) {		// extract the k closest points
	dd[i] = mk.ith_smallest_key(i);
	nn_idx[i] = mk.ith_smallest_info(i);
    }
}
