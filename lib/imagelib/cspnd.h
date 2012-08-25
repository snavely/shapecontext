/* cspnd.h */
/* Find correspondences between two images */

#ifndef __cspnd_h__
#define __cspnd_h__

#ifdef __cplusplus
extern "C" {
#endif

#include "dmap.h"
#include "image.h"

/* Set the number of F-matrix estimation RANSAC trials */
void set_num_fmatrix_trials(int n);

/* Set the similarity threshold for finding symmetric matches */
void cspnd_set_similarity_threshold(double t);

/* Find correspondences from image B to image A */
img_dmap_t *img_find_correspondence(img_t *a, img_t *b, int diameter, img_dmap_t *amap, img_dmap_t *bmap);

/* Estimate correspondences from image B to image A by detecting
 * symmetric matches.  The correspondences are returned in the form of
 * a distance map */
img_dmap_t *img_estimate_correspondence(img_t *a, img_t *b, img_dmap_t *amap, img_dmap_t *bmap);

/* Estimate the fundamental matrix for the images `a' and `b' given an
 * initial set of correspondences `bmap', then filter the set to only
 * contain matches that satisfy the epipolar constraint */
void img_estimate_correspondence_epipolar(img_t *a, img_t *b, img_dmap_t *amap, img_dmap_t *bmap, 
					  double *F, double *e1, double *e2, img_dmap_t **b2a, img_dmap_t **a2b);

#ifdef __cplusplus
}
#endif

#endif /* __cspnd_h__ */
