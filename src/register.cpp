/* register.cpp */
/* Compute relationships between images */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#include "point.h"
#include "register.h"

#include "affine.h"
#include "dmap.h"
#include "homography.h"
#include "horn.h"
#include "matrix.h"
#include "tps.h"
#include "util.h"
#include "vector.h"

// #define LAMBDA 1.0e4 // 1.0e3 // 5.0e1 // 0.0 // 1.0e1
// #define NUM_BASIS_PTS 128 // 512

static int LeastSquaresFit(std::vector<point_t> p1, std::vector<point_t> p2, 
			   std::vector<point_match_t> matches, MotionModel mm,
			   const std::vector<int> &inliers, 
			   MotionParams *params);

static double Ufn(double rsq) 
{
    if (rsq == 0.0)
	return 0.0;
    else
	return rsq * log(rsq);
}

static void EvaluateTPS(std::vector<point_t> &basis_pts, 
		       double *x_affine, double *x_weights,
		       double *y_affine, double *y_weights,
		       double *p, double *q)
{
    int num_basis_pts = (int) basis_pts.size();
    
    double Tx = 
	p[0] * x_affine[0] + p[1] * x_affine[1] + p[2] * x_affine[2];
    double Ty = 
	p[0] * y_affine[0] + p[1] * y_affine[1] + p[2] * y_affine[2];

    for (int j = 0; j < num_basis_pts; j++) {
	double dx = basis_pts[j].x - p[0];
	double dy = basis_pts[j].y - p[1];
	
	Tx += x_weights[j] * Ufn(dx * dx + dy * dy);
	Ty += y_weights[j] * Ufn(dx * dx + dy * dy);
    }

    q[0] = Tx;
    q[1] = Ty;
    q[2] = 1.0;
}

void VectorizeDMAP(img_dmap_t *dmap1, img_dmap_t *dmap2, 
		   std::vector<point_t> &p1, std::vector<point_t> &p2,
		   std::vector<point_match_t> &matches)
{
    p1.clear();
    p2.clear();
    matches.clear();
    
    int idx = 0;
    int count = 0;
    for (int y = 0; y < dmap1->h; y++) {
	for (int x = 0; x < dmap1->w; x++, idx++) {
	    if (dmap1->dists[idx] == DBL_MAX) {
		continue;
	    }

	    point_t pt1, pt2;
	    
	    pt1.x = (double) x;
	    pt1.y = (double) y;
	    
	    pt2.x = (double) Vx(dmap1->nns[idx]);
	    pt2.y = (double) Vy(dmap1->nns[idx]);

	    p1.push_back(pt1);
	    p2.push_back(pt2);

	    point_match_t m;
	    m.idx1 = count;
	    m.idx2 = count;

	    matches.push_back(m);

	    count++;
	}
    }

    printf("num_matches = %d\n", count);
}

void PruneDMAP(img_dmap_t *dmap1, img_dmap_t *dmap2, 
	       std::vector<point_t> p1, std::vector<point_t> p2, 
	       std::vector<point_match_t> matches, std::vector<int> inliers,
	       img_dmap_t *dmap1_out, img_dmap_t *dmap2_out)
{
    int idx = 0;
    for (int y = 0; y < dmap1_out->h; y++) {
	for (int x = 0; x < dmap1_out->w; x++, idx++) {
	    dmap1_out->dists[idx] = DBL_MAX;
	}
    }

    idx = 0;
    for (int y = 0; y < dmap2_out->h; y++) {
	for (int x = 0; x < dmap2_out->w; x++, idx++) {
	    dmap2_out->dists[idx] = DBL_MAX;
	}
    }

    for (int i = 0; i < (int) inliers.size(); i++) {
	int midx = inliers[i];
	int idx1 = matches[midx].idx1;
	int idx2 = matches[midx].idx2;
	    
	int x1 = iround(p1[idx1].x);
	int y1 = iround(p1[idx1].y);
	
	int x2 = iround(p2[idx2].x);
	int y2 = iround(p2[idx2].y);

	int didx1 = y1 * dmap1_out->w + x1;
	int didx2 = y2 * dmap2_out->w + x2;

	dmap1_out->dists[didx1] = dmap1->dists[didx1];
	dmap1_out->nns[didx1]   = dmap1->nns[didx1];
	    
	dmap2_out->dists[didx2] = dmap2->dists[didx2];
	dmap2_out->nns[didx2]   = dmap2->nns[didx2];
    }    
}

    
/* Wrapper for dmaps */
std::vector<int> EstimateTransform(img_dmap_t *dmap1, img_dmap_t *dmap2, 
				   MotionModel mm,
				   int nRANSAC, double RANSACthresh,
				   MotionParams *params,
				   img_dmap_t *dmap1_out, 
				   img_dmap_t *dmap2_out, bool output_full)
{
    std::vector<point_t> p1, p2;
    std::vector<point_match_t> matches;
    
    int idx = 0;
    int count = 0;
    for (int y = 0; y < dmap1->h; y++) {
	for (int x = 0; x < dmap1->w; x++, idx++) {
	    if (dmap1->dists[idx] == DBL_MAX) {
		continue;
	    }

	    point_t pt1, pt2;
	    
	    pt1.x = (double) x;
	    pt1.y = (double) y;
	    
	    pt2.x = (double) Vx(dmap1->nns[idx]);
	    pt2.y = (double) Vy(dmap1->nns[idx]);

	    p1.push_back(pt1);
	    p2.push_back(pt2);

	    point_match_t m;
	    m.idx1 = count;
	    m.idx2 = count;

	    matches.push_back(m);

	    count++;
	}
    }

    printf("num_matches = %d\n", count);

    std::vector<int> inliers;
    inliers = EstimateTransform(p1, p2, matches, mm, 
				nRANSAC, RANSACthresh, params);

    /* Create the dmaps */
    if (!output_full) {
	idx = 0;
	for (int y = 0; y < dmap1_out->h; y++) {
	    for (int x = 0; x < dmap1_out->w; x++, idx++) {
		dmap1_out->dists[idx] = DBL_MAX;
	    }
	}

	idx = 0;
	for (int y = 0; y < dmap2_out->h; y++) {
	    for (int x = 0; x < dmap2_out->w; x++, idx++) {
		dmap2_out->dists[idx] = DBL_MAX;
	    }
	}

	for (int i = 0; i < (int) inliers.size(); i++) {
	    int midx = inliers[i];
	    int idx1 = matches[midx].idx1;
	    int idx2 = matches[midx].idx2;
	    
	    int x1 = iround(p1[idx1].x);
	    int y1 = iround(p1[idx1].y);

	    int x2 = iround(p2[idx2].x);
	    int y2 = iround(p2[idx2].y);

	    int didx1 = y1 * dmap1_out->w + x1;
	    int didx2 = y2 * dmap2_out->w + x2;

	    dmap1_out->dists[didx1] = dmap1->dists[didx1];
	    dmap1_out->nns[didx1]   = dmap1->nns[didx1];

	    dmap2_out->dists[didx2] = dmap2->dists[didx2];
	    dmap2_out->nns[didx2]   = dmap2->nns[didx2];
	}
    } else {
	idx = 0;
	for (int y = 0; y < dmap1_out->h; y++) {
	    for (int x = 0; x < dmap1_out->w; x++, idx++) {
		dmap1_out->dists[idx] = DBL_MAX;
	    }
	}

	idx = 0;
	for (int y = 0; y < dmap2_out->h; y++) {
	    for (int x = 0; x < dmap2_out->w; x++, idx++) {
		dmap2_out->dists[idx] = DBL_MAX;
	    }
	}

	for (int y = 0; y < dmap1_out->h; y++) {
	    for (int x = 0; x < dmap1_out->w; x++, idx++) {
		double p[3] = { (double) x, (double) y, 1.0 };
		double q[3];
		
		EvaluateTPS(params->basis_pts, 
			    params->x_affine, params->x_weights,
			    params->y_affine, params->y_weights,
			    p, q);

		int didx1 = y * dmap1_out->w + x;
		
		dmap1_out->dists[didx1] = 0.0;
		dmap1_out->nns[didx1] = v2_new(q[0], q[1]);

		if (q[0] > 0.0 && q[0] < dmap2_out->w - 1 &&
		    q[1] > 0.0 && q[1] < dmap2_out->h - 1) {
		    
		    int didx2 = iround(q[1]) * dmap2_out->w + iround(q[0]);
		    dmap2_out->dists[didx2] = 0.0;
		    dmap2_out->nns[didx2] = v2_new(p[0], p[1]);
		}
	    }
	}

#if 1
	/* Invert the mapping */
	v3_t *r_pts = new v3_t[params->num_basis_pts];
	v3_t *l_pts = new v3_t[params->num_basis_pts];

	for (int i = 0; i < params->num_basis_pts; i++) {
	    l_pts[i] = v3_new(params->mapped_pts[i].x, 
			      params->mapped_pts[i].y, 1.0);
	    r_pts[i] = v3_new(params->basis_pts[i].x, 
			      params->basis_pts[i].y, 1.0);
	}

	double x_affine[3], y_affine[3];
	double *x_weights = new double[params->num_basis_pts];
	double *y_weights = new double[params->num_basis_pts];

	align_tps(params->num_basis_pts, r_pts, l_pts, params->lambda,
		  x_affine, x_weights, y_affine, y_weights);
	
	for (int y = 0; y < dmap2_out->h; y++) {
	    for (int x = 0; x < dmap2_out->w; x++, idx++) {
		double p[3] = { (double) x, (double) y, 1.0 };
		double q[3];
		
		EvaluateTPS(params->mapped_pts, 
			    x_affine, x_weights,
			    y_affine, y_weights,
			    p, q);

		int didx2 = y * dmap2_out->w + x;
		
		dmap2_out->dists[didx2] = 0.0;
		dmap2_out->nns[didx2] = v2_new(q[0], q[1]);
	    }
	}

	delete [] r_pts;
	delete [] l_pts;
	delete [] x_weights;
	delete [] y_weights;
#endif
    }
    
    return inliers;
}


/* Estimate a transform between two sets of keypoints */
std::vector<int> EstimateTransform(std::vector<point_t> p1, 
				   std::vector<point_t> p2, 
				   std::vector<point_match_t> matches, 
				   MotionModel mm,
				   int nRANSAC, double RANSACthresh, 
				   MotionParams *params) 
{
    params->model = mm;

    int min_matches = -1;
    switch (mm) {
	case MotionRigid:
	    min_matches = 3;
	    break;
	case MotionAffine:
	    min_matches = 3;
	    break;
	case MotionHomography:
	    min_matches = 4;
	    break;
	case MotionThinPlateSpline:
	    min_matches = params->num_basis_pts; // NUM_BASIS_PTS;
	    break;
    }

    int *match_idxs = new int[min_matches];

    int num_matches = (int) matches.size();
    int max_inliers = 0;
    MotionParams Mbest;

    Mbest.model = mm;
    Mbest.x_weights = NULL;
    Mbest.y_weights = NULL;
    
    if (num_matches < min_matches) {
	std::vector<int> empty;
	printf("Cannot estimate rigid transform\n");
	return empty;
    }

    v3_t *r_pts = new v3_t[min_matches];
    v3_t *l_pts = new v3_t[min_matches];
    double *weight = new double[min_matches];

    for (int round = 0; round < nRANSAC; round++) {
	for (int i = 0; i < min_matches; i++) {
	    bool found;
	    int idx;
	    
	    do {
		found = true;
		idx = rand() % num_matches;
		
		for (int j = 0; j < i; j++) {
		    if (match_idxs[j] == idx) {
			found = false;
			break;
		    }
		}
	    } while (!found);

	    match_idxs[i] = idx;
	}

	/* Solve for the motion */
	std::vector<point_t> basis_pts;
	std::vector<point_t> mapped_pts;

	for (int i = 0; i < min_matches; i++) {
	    int idx1 = matches[match_idxs[i]].idx1;
	    int idx2 = matches[match_idxs[i]].idx2;
	    
	    Vx(l_pts[i]) = p1[idx1].x;
	    Vy(l_pts[i]) = p1[idx1].y;
	    Vz(l_pts[i]) = 1.0;
		    
	    Vx(r_pts[i]) = p2[idx2].x;
	    Vy(r_pts[i]) = p2[idx2].y;
	    Vz(r_pts[i]) = 1.0;

	    weight[i] = 1.0;

	    point_t pt;
	    pt.x = Vx(l_pts[i]);
	    pt.y = Vy(l_pts[i]);
	    
	    basis_pts.push_back(pt);

	    pt.x = Vx(r_pts[i]);
	    pt.y = Vy(r_pts[i]);
	    mapped_pts.push_back(pt);
	}

	double Mcurr[9];
	double *x_weights = new double[min_matches], x_affine[3];
	double *y_weights = new double[min_matches], y_affine[3];

	switch (mm) {
	    case MotionRigid: {
		double R[9], T[9], Tout[9], scale;
		align_horn(min_matches, r_pts, l_pts, R, T, Tout, &scale, weight);
		memcpy(Mcurr, Tout, 9 * sizeof(double));
		break;
	    }

	    case MotionAffine: {
		align_affine(min_matches, r_pts, l_pts, Mcurr);
		break;
	    }

	    case MotionThinPlateSpline: {
		align_tps(min_matches, r_pts, l_pts, params->lambda,
			  x_affine, x_weights, y_affine, y_weights);

		printf("[%d] x_affine: %0.3e, %0.3e, %0.3e\n", 
		       round, x_affine[0], x_affine[1], x_affine[2]);

		double error = 0.0;
		for (int i = 0; i < min_matches; i++) {
		    double p[3];

		    p[0] = Vx(l_pts[i]);
		    p[1] = Vy(l_pts[i]);
		    p[2] = 1.0;
		    
		    double q[3];

		    EvaluateTPS(basis_pts, 
				x_affine, x_weights, y_affine, y_weights,
				p, q);

		    double dx = q[0] - Vx(r_pts[i]);
		    double dy = q[1] - Vy(r_pts[i]);
		    
		    double resid = dx * dx + dy * dy;
		    error += resid;
		}

		printf("error = %0.3e\n", error / params->num_basis_pts);

		// matrix_print(1, min_matches, x_weights);

		break;
	    }
		
	    case MotionHomography: {
		align_homography(min_matches, r_pts, l_pts, Mcurr, 0);
		break;
	    }
	}
		

	std::vector<int> inliers;
	
	int num_inliers;

	if (mm != MotionThinPlateSpline) {
	    num_inliers = CountInliers(p1, p2, matches, Mcurr, 
				       RANSACthresh, inliers);
	} else {
	    num_inliers = CountInliers(p1, p2, matches, 
				       basis_pts,
				       x_affine, x_weights, 
				       y_affine, y_weights,
				       RANSACthresh, inliers);
	    printf("[%d] num_inliers = %d\n", round, num_inliers);
	    fflush(stdout);
	}
	
	if (num_inliers > max_inliers) {
	    max_inliers = num_inliers;

	    // memcpy(Mbest, Mcurr, 9 * sizeof(double));
	    memcpy(Mbest.M, Mcurr, 9 * sizeof(double));
	    memcpy(Mbest.x_affine, x_affine, 3 * sizeof(double));
	    memcpy(Mbest.y_affine, y_affine, 3 * sizeof(double));

	    Mbest.basis_pts = basis_pts;
	    Mbest.mapped_pts = mapped_pts;

	    if (Mbest.x_weights != NULL)
		delete [] Mbest.x_weights;
	    if (Mbest.y_weights != NULL)
		delete [] Mbest.y_weights;
	    
	    Mbest.x_weights = new double[basis_pts.size()];
	    Mbest.y_weights = new double[basis_pts.size()];
	    
	    memcpy(Mbest.x_weights, x_weights, 
		   basis_pts.size() *sizeof(double));
	    memcpy(Mbest.y_weights, y_weights, 
		   basis_pts.size() *sizeof(double));	    
	}

	delete [] x_weights;
	delete [] y_weights;
    }

    std::vector<int> inliers;
    if (mm != MotionThinPlateSpline) {
	CountInliers(p1, p2, matches, Mbest.M, RANSACthresh, inliers);
    } else {
	CountInliers(p1, p2, matches, 
		     Mbest.basis_pts, 
		     Mbest.x_affine, Mbest.x_weights, 
		     Mbest.y_affine, Mbest.y_weights,
		     RANSACthresh, inliers);

	printf("[final] x_affine = %0.3e, %0.3e, %0.3e\n", 
	       Mbest.x_affine[0], Mbest.x_affine[1], Mbest.x_affine[2]);

	printf("[final] y_affine = %0.3e, %0.3e, %0.3e\n", 
	       Mbest.y_affine[0], Mbest.y_affine[1], Mbest.y_affine[2]);

	printf("[final] x_weights: ");
	matrix_print(1, Mbest.basis_pts.size(), Mbest.x_weights);

	printf("[final] y_weights: ");
	matrix_print(1, Mbest.basis_pts.size(), Mbest.y_weights);

	printf("[final] num_inliers = %d\n", (int) inliers.size());
    }

    // CountInliers(p1, p2, matches, Mbest, RANSACthresh, inliers);

    if (mm != MotionThinPlateSpline) {
	LeastSquaresFit(p1, p2, matches, mm, inliers, params);
    } else {
	memcpy(params->x_affine, Mbest.x_affine, 3 * sizeof(double));
	memcpy(params->y_affine, Mbest.y_affine, 3 * sizeof(double));

	params->basis_pts = Mbest.basis_pts;
	params->mapped_pts = Mbest.mapped_pts;

	params->x_weights = new double[Mbest.basis_pts.size()];
	params->y_weights = new double[Mbest.basis_pts.size()];
	    
	memcpy(params->x_weights, Mbest.x_weights, 
	       Mbest.basis_pts.size() *sizeof(double));
	memcpy(params->y_weights, Mbest.y_weights, 
	       Mbest.basis_pts.size() *sizeof(double));
    }
    
    delete [] match_idxs;
    delete [] r_pts;
    delete [] l_pts;
    delete [] weight;

    return inliers;
}

int CountInliers(std::vector<point_t> p1, std::vector<point_t> p2, 
		 std::vector<point_match_t> matches,
		 double *M, double thresh, std::vector<int> &inliers)
{
    inliers.clear();
    int count = 0;

    for (unsigned int i = 0; i < matches.size(); i++) {
	/* Determine if the ith feature in f1, when transformed by M,
	 * is within RANSACthresh of its match in f2 (if one exists)
	 *
	 * if so, increment count and append i to inliers */

	double p[3];

	p[0] = p1[matches[i].idx1].x;
	p[1] = p1[matches[i].idx1].y;
	p[2] = 1.0;

	double q[3];
	matrix_product(3, 3, 3, 1, M, p, q);

	double qx = q[0] / q[2];
	double qy = q[1] / q[2];

	double dx = qx - p2[matches[i].idx2].x;
	double dy = qy - p2[matches[i].idx2].y;
	
	double dist = sqrt(dx * dx + dy * dy);
	
	if (dist <= thresh) {
	    count++;
	    inliers.push_back(i);
	}
    }

    return count;
}

int CountInliers(std::vector<point_t> p1, std::vector<point_t> p2, 
		 std::vector<point_match_t> matches,
		 std::vector<point_t> basis_pts,
		 double *x_affine, double *x_weights,
		 double *y_affine, double *y_weights,
		 double thresh, std::vector<int> &inliers)
{
    inliers.clear();
    int count = 0;

    int num_basis_pts = (int) basis_pts.size();

    for (unsigned int i = 0; i < matches.size(); i++) {
	/* Determine if the ith feature in f1, when transformed by M,
	 * is within RANSACthresh of its match in f2 (if one exists)
	 *
	 * if so, increment count and append i to inliers */

	double p[3];

	p[0] = p1[matches[i].idx1].x;
	p[1] = p1[matches[i].idx1].y;
	p[2] = 1.0;

	double Tx = 
	    p[0] * x_affine[0] + p[1] * x_affine[1] + p[2] * x_affine[2];
	double Ty = 
	    p[0] * y_affine[0] + p[1] * y_affine[1] + p[2] * y_affine[2];

	for (int j = 0; j < num_basis_pts; j++) {
	    double dx = basis_pts[j].x - p[0];
	    double dy = basis_pts[j].y - p[1];

	    Tx += x_weights[j] * Ufn(dx * dx + dy * dy);
	    Ty += y_weights[j] * Ufn(dx * dx + dy * dy);
	}

	double dx = Tx - p2[matches[i].idx2].x;
	double dy = Ty - p2[matches[i].idx2].y;
	
	double dist = sqrt(dx * dx + dy * dy);
	
	if (dist <= thresh) {
	    count++;
	    inliers.push_back(i);
	}
    }

    return count;
}


static int LeastSquaresFit(std::vector<point_t> p1, std::vector<point_t> p2, 
			   std::vector<point_match_t> matches, MotionModel mm,
			   const std::vector<int> &inliers, 
			   MotionParams *params)
{
    v3_t *r_pts = new v3_t[inliers.size()];
    v3_t *l_pts = new v3_t[inliers.size()];
    double *weight = new double[inliers.size()];
    std::vector<point_t> basis_pts;

    for (int i=0; i < (int) inliers.size(); i++) {
	int idx1 = matches[inliers[i]].idx1;
	int idx2 = matches[inliers[i]].idx2;
	
	Vx(l_pts[i]) = p1[idx1].x;
	Vy(l_pts[i]) = p1[idx1].y;
	Vz(l_pts[i]) = 1.0;
		    
	Vx(r_pts[i]) = p2[idx2].x;
	Vy(r_pts[i]) = p2[idx2].y;
	Vz(r_pts[i]) = 1.0;

	weight[i] = 1.0;

	point_t pt;
	pt.x = Vx(l_pts[i]);
	pt.y = Vy(l_pts[i]);
	    
	basis_pts.push_back(pt);
    }
    
    switch (mm) {
	case MotionRigid: {
	    double R[9], T[9], Tout[9], scale;
	    align_horn(inliers.size(), r_pts, l_pts, R, T, Tout, &scale, weight);
	    memcpy(params->M, Tout, 9 * sizeof(double));
	    break;
	}

	case MotionAffine: {
	    align_affine(inliers.size(), r_pts, l_pts, params->M);
	    break;
	}

	case MotionThinPlateSpline: {
	    double *x_weights = new double[inliers.size()], x_affine[3];
	    double *y_weights = new double[inliers.size()], y_affine[3];

	    align_tps(inliers.size(), r_pts, l_pts, params->lambda,
		      x_affine, x_weights, y_affine, y_weights);

	    memcpy(params->x_affine, x_affine, 3 * sizeof(double));
	    memcpy(params->y_affine, y_affine, 3 * sizeof(double));

	    params->basis_pts = basis_pts;
	    
	    params->x_weights = new double[basis_pts.size()];
	    params->y_weights = new double[basis_pts.size()];
	    
	    memcpy(params->x_weights, x_weights, 
		   basis_pts.size() *sizeof(double));
	    memcpy(params->y_weights, y_weights, 
		   basis_pts.size() *sizeof(double));	    

	    delete [] x_weights;
	    delete [] y_weights;

	    break;
	}
	    
	case MotionHomography: {
	    align_homography(inliers.size(), r_pts, l_pts, params->M, 1);
	    break;
	}
    }

    delete [] r_pts;
    delete [] l_pts;
    delete [] weight;

    return 0;
}

