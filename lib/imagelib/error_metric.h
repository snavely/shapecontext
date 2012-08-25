/* error.h */

#ifndef __error_h__
#define __error_h__

#include "anniface.h"
#include "image.h"
#include "transform.h"
#include "vector.h"

typedef enum {
    CS_RGB,
    CS_LUV
} color_space_t;

typedef enum {
    EM_GRAYSCALE,      /* (x, y, intensity) */
    EM_GRAYSCALE_FILL, /* (x, y, intensity), interpolated for
			* in-between points */
    EM_RGB,            /* (x, y, r, g, b) */
    EM_GRADIENT,       /* (x, y, intensity, gradient) */
    EM_NHOOD_GS,       /* (x, y, i1, i2, ... in) */
    EM_NHOOD_RGB,      /* (x, y, r1, g1, b1, r2, g2, b2, ...) */
    EM_NHOOD_RGB_NORM, /* Normalized version of the above */
    EM_HISTOGRAM_GS,   /* (x, y, b1, b2, ..., bk) */
    EM_HISTOGRAM_RGB,  /* (x, y, br1, ..., brk, bg1, ..., bgk,
			*  bb1, ..., bbk) */
    EM_RING_HISTOGRAM_GS,  /* (x, y, bi1, ..., bik) */
    EM_RING_HISTOGRAM_RGB,  /* (x, y, br1, ..., brk, bg1, ..., bgk,
			     *  bb1, ..., bbk, br1, ..., brk, bg1, ..., bgk,
			     *  bb1, ..., bbk, ...) */
    EM_UNKNOWN
} error_metric_t;

/* Set the values of certain error metric parameters */
void em_set_nhood_size(int size);
void em_set_nhood_dims(int w, int h);
void em_set_hgram_nhood_radius(double radius);
void em_set_num_hgram_buckets(int num_buckets);
void em_set_grid_fill_ratio(int ratio);
void em_set_num_ring_hgrams(int n);

void em_use_RGB();
void em_use_LUV();
color_space_t em_get_color_space();

/* Create a query point for the nearest-neighbor search for the
 * current error metric */
vec_t em_create_query_pt(img_t *img, error_metric_t em, double x, double y);

/* Create a query point for the nearest-neighbor search for the
 * current error metric */
vec_t em_create_query_pt_downsize(img_t *img, error_metric_t em, double x, double y, int downsize);

/* Create a tree for an image based on the given error metric */
ANNkd_tree_t *img_create_ann_tree(img_t *img, error_metric_t em, int *nhood_radius);
ANNkd_tree_t *img_create_ann_tree_params(img_t *img, error_metric_t em, int *nhood_radius, int rectified, int downsize);
ANNkd_tree_t *img_create_ann_tree_rectified(img_t *img, error_metric_t em, int *nhood_radius);

/* Create a tree for an image based on the given error metric, and
 * make an index map along with it */
ANNkd_tree_t *img_create_ann_tree_idx_map(img_t *img, error_metric_t em, int *nhood_radius, int *idx_map);
ANNkd_tree_t *img_create_ann_tree_idx_map_rectified(img_t *img, error_metric_t em, int *nhood_radius, int *idx_map);

/* Create a tree for part of an image based on the given error metric */
ANNkd_tree_t *img_create_ann_subtree(img_t *img, error_metric_t em, int *nhood_radius,
				     int xmin, int xmax, int ymin, int ymax, int rectified, int downsize,
				     ANNkd_tree_t *parent, int *idx_map);

/* Compute the (asymmetric) distance between image a and image b */
double img_distance(img_t *a, img_t *b);

/* Produces an image of the chamfer metric between two images */
img_t *img_error(img_t *a, img_t *b, error_metric_t em);

void set_error_x_weight(double w);
double get_error_x_weight();

void set_error_y_weight(double w);
double get_error_y_weight();

#endif /* __error_h__ */
