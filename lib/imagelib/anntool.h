#ifndef __anntool_h__
#define __anntool_h__

#ifdef __cplusplus
extern "C" {
#endif

#include "anniface.h"
#include "binvolume.h"

/* Create a kd tree of the points stored in a binary volume */
ANNkd_tree_t *bvol_create_kd_tree(bvol_t *b);

/* Create a 3D kd-tree from a (grayscale) image */
ANNkd_tree_t *img_create_kd_tree_grayscale(img_t *img);

/* Create a 3D kd-tree from a (grayscale) image, filling in some
 * interior points using linear interpolation */
ANNkd_tree_t *img_create_kd_tree_grayscale_fill(img_t *img, int fill_ratio);

/* Create a 4D kd-tree from a grayscale image.  The fourth dimension
 * is the gradient */
ANNkd_tree_t *img_create_kd_tree_gradient(img_t *img);

/* Create a 2+n^2 dimensional kd-tree from a grayscale image, 
 * with dimensions for each pixel value in the neighborhood */
ANNkd_tree_t *img_create_kd_tree_nhood(img_t *img, int n);

/* Create a 2+3*n^2 dimensional kd-tree from a rgb image, 
 * with dimensions for each pixel value in the neighborhood */
ANNkd_tree_t *img_create_kd_tree_nhood_rgb(img_t *img, int n, int fill_ratio,
					   int xmin, int xmax, int ymin, int ymax);

/* Create a 5D kd-tree from a (rgb) image, with separate dimensions
 * for each of the tree color channels */
ANNkd_tree_t *img_create_kd_tree_rgb(img_t *img, int fill_ratio);

/* Return the approximate nearest neighbor to q out of the points in
 * the (weighted) tree */
v3_t query_weighted_ann_tree_3D_pt(ANNkd_tree_t *tree, v3_t q, double eps);

#ifdef __cplusplus
}
#endif

#endif /* __anntool_h__ */
