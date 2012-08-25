//----------------------------------------------------------------------
//	File:		ann_sample.cc
//	Programmer:	Sunil Arya and David Mount
//	Last modified:	03/04/98 (Release 0.1)
//	Description:	Sample program for ANN
//----------------------------------------------------------------------
// Copyright (c) 1997-1998 University of Maryland and Sunil Arya and David
// Mount All Rights Reserved.
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

#include <stdio.h>			// C I/O
#include <fstream>			// file I/O
#include <string.h>			// string manipulation
#include <math.h>			// math routines
#include <ANN/ANN.h>			// ANN declarations

//----------------------------------------------------------------------
// ann_sample
//
// This is a simple sample program for the ANN library.  After compiling,
// it can be run as follows.
// 
// 	ann_sample [-d dim] [-max mpts] [-nn k] [-e eps] [-df data] [-qf query]
//
// where
//	dim		is the dimension of the space (default = 2)
//	mpts		maximum number of data points (default = 1000)
//	k		number of nearest neighbors per query (default 1)
//	eps		is the error bound (default = 0.0)
//	data		file containing data points
//	query		file containing query points
//
// Results are sent to the standard output.
//----------------------------------------------------------------------

//----------------------------------------------------------------------
//  Parameters that are set in getArgs()
//----------------------------------------------------------------------
void getArgs(int argc, char **argv);		// get command-line arguments

int		k		= 1;		// number of nearest neighbors
int		dim		= 2;		// dimension
double		eps		= 0;		// error bound
int		m_pts		= 1000;		// maximum number of data points

istream		*data_in	= NULL;		// input for data points
istream		*query_in	= NULL;		// input for query points

ANNbool readPt(istream &in, ANNpoint p)		// read point (false on EOF)
{
    for (int i = 0; i < dim; i++) {
	if(!(in >> p[i])) return ANNfalse;
    }
    return ANNtrue;
}

void printPt(ostream &out, ANNpoint p)		// print point
{
    out << "(" << p[0];
    for (int i = 1; i < dim; i++) {
	out << ", " << p[i];
    }
    out << ")\n";
}

main(int argc, char **argv)
{
    int			n_pts;			// actual number of data points
    ANNpointArray	data_pts;		// data points
    ANNpoint		query_pt;		// query point
    ANNidxArray		nn_idx;			// near neighbor indices
    ANNdistArray	dists;			// near neighbor distances
    ANNkd_tree		*the_tree;		// search structure

    getArgs(argc, argv);			// read command-line arguments

    query_pt = annAllocPt(dim);			// allocate query point
    data_pts = annAllocPts(m_pts, dim);		// allocate data points
    nn_idx = new ANNidx[k];			// allocate near neigh indices
    dists = new ANNdist[k];			// allocate near neighbor dists

    n_pts = 0;					// read data points

    cout << "Data Points:\n";
    while (n_pts < m_pts && readPt(*data_in, data_pts[n_pts])) {
	printPt(cout, data_pts[n_pts]);
	n_pts++;
    }

    the_tree = new ANNkd_tree(			// build search structure
		    data_pts,			// the data points
		    n_pts,			// number of points
		    dim);			// dimension of space

    while (readPt(*query_in, query_pt)) {	// read query points
	cout << "Query point: ";		// echo query point
	printPt(cout, query_pt);

	the_tree->annkSearch(			// search
		query_pt,			// query point
		k,				// number of near neighbors
		nn_idx,				// nearest neighbors (returned)
     		dists,				// distance (returned)
		eps);				// error bound

	cout << "\tNN:\tIndex\tDistance\n";
	for (int i = 0; i < k; i++) {		// print summary
	    dists[i] = sqrt(dists[i]);		// unsquare distance
	    cout << "\t" << i << "\t" << nn_idx[i] << "\t" << dists[i] << "\n";
	}
    }
}

//----------------------------------------------------------------------
//  getArgs - get command line arguments
//----------------------------------------------------------------------

void getArgs(int argc, char **argv)
{
    static ifstream dataStream;			// data file stream
    static ifstream queryStream;		// query file stream

    if (argc <= 1) {				// no arguments
  	cerr << "Usage:\n\n"
	<< "  ann_sample [-d dim] [-max m] [-nn k] [-e eps] [-df data]"
	   " [-qf query]\n\n"
	<< "  where:\n"
	<< "    dim             dimension of the space (default = 2)\n"
	<< "    m               maximum number of data points (default = 1000)\n"
	<< "    k               number of nearest neighbors per query (default 1)\n"
	<< "    eps             the error bound (default = 0.0)\n"
	<< "    data            name of file containing data points\n"
	<< "    query           name of file containing query points\n\n"
	<< " Results are sent to the standard output.\n"
	<< "\n"
	<< " To run this demo use:\n"
	<< "    ann_sample -df data_pts -qf query_pts\n";
	exit(0);
    }
    int i = 1;
    while (i < argc) {				// read arguments
	if (!strcmp(argv[i], "-d")) {		// -d option
	    dim = atoi(argv[++i]);		// get dimension to dump
	}
	else if (!strcmp(argv[i], "-max")) {	// -max option
	    m_pts = atoi(argv[++i]);		// get max number of points
	}
	else if (!strcmp(argv[i], "-nn")) {	// -nn option
	    k = atoi(argv[++i]);		// get number of near neighbors
	}
	else if (!strcmp(argv[i], "-e")) {	// -e option
	    sscanf(argv[++i], "%lf", &eps);	// get error bound
	}
	else if (!strcmp(argv[i], "-df")) {	// -df option
	    dataStream.open(argv[++i], ios::in);// open data file
	    if (!dataStream) {
		cerr << "Cannot open data file\n";
		exit(1);
	    }
	    data_in = &dataStream;		// make this the data stream
	}
	else if (!strcmp(argv[i], "-qf")) {	// -qf option
	    queryStream.open(argv[++i], ios::in);// open query file
	    if (!queryStream) {
		cerr << "Cannot open query file\n";
		exit(1);
	    }
	    query_in = &queryStream;		// make this query stream
	}
	else {					// illegal syntax
	    cerr << "Unrecognized option.\n";
	    exit(1);
	}
	i++;
    }
    if (data_in == NULL || query_in == NULL) {
	cerr << "-df and -qf options must be specified\n";
	exit(1);
    }
}
