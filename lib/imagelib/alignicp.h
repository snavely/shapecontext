/* alignicp.h */
/* Routines for aligning pairs of images using the Iterative Closest Point
 * algorithm */

#ifndef __alignicp_h__
#define __alignicp_h__

#include "error_metric.h"
#include "transform.h"
#include "transform3D.h"
#include "transform-opt.h"

/* Set the error metric functions used by icp */
int icp_set_error_metric(error_metric_t em);

/* Find a linear transformation that aligns volume B to volume A */
trans2D_t *align_bin_volume_ICP(bvol_t *a, bvol_t *b);
trans2D_t *align_bin_volume_ICP2(bvol_t *a, bvol_t *b);
trans2D_t *align_bin_volume_ICP3(bvol_t *a, bvol_t *b);

/* Find a linear transformation that aligns image B to image A.
 * Uses ICP and resamples the image after each iteration.  The
 * transform that describes the initial guess is given. */
trans2D_t *align_Timage_ICP(img_t *a, img_t *b, trans2D_t *Tin, transform_class_t tclass);

/* Find a linear transformation that aligns image B to image A.
 * Uses ICP and resamples the image after each iteration.  The
 * transform that describes the initial guess is given. */
trans2D_t *align_Timage_ICP_symmetric(img_t *a, img_t *b, trans2D_t *Tin);

/* Find a linear transformation that aligns image B to image A.
 * Uses ICP and resamples the image after each iteration. */
trans2D_t *align_image_ICP(img_t *a, img_t *b, transform_class_t tclass);

/* Set the percent of worst matches that will be rejected during ICP */
void set_icp_reject_percent(double p);

/* Set the epsilon used in the nearest neighbor calculations */
void set_icp_ann_epsilon(double eps);

/* Set the percent of pairs which are randomly sampled during ICP */
void set_icp_pairs_sample_percent(double p);

/* Set the maximum number of iterations for ICP to execute */
void set_icp_max_rounds(int max_rounds);

/* Set whether we just want to sample the intersection of two images */
void set_icp_sample_intersection(int val);

/* Turn off/on recursive ICP */
void icp_recurse_off();
void icp_recurse_on();

#endif /* __alignicp_h__ */
