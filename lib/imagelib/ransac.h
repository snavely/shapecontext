/* ransac.h */
/* Perform the RANSAC (RANdom SAmple Consensus) algorithm to align a
 * pair of images */

#ifndef __ransac_h__
#define __ransac_h__

#include "image.h"
#include "pyramid.h"
#include "transform.h"

/* Compute the transform that aligns image B to image A using the
 * RANSAC method */
trans2D_t *align_image_ransac(img_t *a, img_t *b, int use_icp, double *error_out);

/* This is the same function as above but takes a distance map for a
 * and b as input */
trans2D_t *align_image_ransac_map(img_t *a, img_t *b, int diameter, 
				  img_dmap_t *amap, img_dmap_t *bmap, int use_icp, double *error_out);

/* Align image B to image B given an initial estimate of the transform */
trans2D_t *align_Timage_ransac(img_t *a, img_t *b, trans2D_t *Tin, int use_icp, double *error_out);

/* This is the same function as above but takes a distance map for a
 * and b as input */
trans2D_t *align_Timage_ransac_map(img_t *a, img_t *b, trans2D_t *Tin, int diameter, 
				   img_dmap_t *amap, img_dmap_t *bmap, int use_icp, double *error_out);

/* Set the maximum number of trials that RANSAC will use */
void set_ransac_max_trials(int n);

#endif /* __ransac_h__ */
