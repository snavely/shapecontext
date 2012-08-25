/* register.h */
/* Compute relationships between images */

#ifndef __register_h__
#define __register_h__

#include <vector>

#include "point.h"

#include "dmap.h"

enum MotionModel {
    MotionRigid,
    MotionAffine,
    MotionThinPlateSpline,
    MotionHomography,
};

typedef struct {
    MotionModel model;
    double M[9];  /* For linear models */

    std::vector<point_t> basis_pts;
    std::vector<point_t> mapped_pts;

    double x_affine[3], y_affine[3];  /* For TPS */
    double *x_weights, *y_weights;
    double lambda;
    int num_basis_pts;
    
} MotionParams;


void VectorizeDMAP(img_dmap_t *dmap1, img_dmap_t *dmap2, 
		   std::vector<point_t> &p1, std::vector<point_t> &p2,
		   std::vector<point_match_t> &matches);

void PruneDMAP(img_dmap_t *dmap1, img_dmap_t *dmap2, 
	       std::vector<point_t> p1, std::vector<point_t> p2, 
	       std::vector<point_match_t> matches, std::vector<int> inliers,
	       img_dmap_t *dmap1_out, img_dmap_t *dmap2_out);

/* Wrapper for dmaps */
std::vector<int> EstimateTransform(img_dmap_t *dmap1, img_dmap_t *dmap2, 
				   MotionModel mm,
				   int nRANSAC, double RANSACthresh,
				   MotionParams *params,
				   img_dmap_t *dmap1_out, 
				   img_dmap_t *dmap2_out, bool output_full);

/* Estimate a transform between two sets of keypoints */
std::vector<int> EstimateTransform(std::vector<point_t> k1, 
                                   std::vector<point_t> k2, 
                                   std::vector<point_match_t> matches, 
                                   MotionModel mm,
                                   int nRANSAC, double RANSACthresh, 
                                   MotionParams *params);

int CountInliers(std::vector<point_t> p1, std::vector<point_t> p2, 
		 std::vector<point_match_t> matches,
		 double *M, double thresh, std::vector<int> &inliers);

int CountInliers(std::vector<point_t> p1, std::vector<point_t> p2, 
		 std::vector<point_match_t> matches,
		 std::vector<point_t> basis_pts,
		 double *x_affine, double *x_weights,
		 double *y_affine, double *y_weights,
		 double thresh, std::vector<int> &inliers);

#endif /* __register_h__ */
