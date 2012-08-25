/* pyramid.h */
/* Routines for creating a gaussian pyramid from an image */

#ifndef __pyramid_h__
#define __pyramid_h__

#ifdef __cplusplus
extern "C" {
#endif

#include <sys/types.h>

#include "anniface.h"
#include "dmap.h"
#include "error.h"
#include "image.h"
#include "vector.h"

typedef struct {
    u_int16_t num_levels;  /* Number of levels of the pyramid */
    u_int16_t w, h;        /* Width and height of the bottom of the pyramid */
    img_t **imgs;          /* Images in the pyramid */
} img_pyr_t;

typedef struct {
    u_int16_t num_levels;
    ANNkd_tree_t **trees;
} ANNkd_pyramid_t;

typedef struct {

    /* Information about the metric used */
    u_int16_t em;
    u_int16_t nhood_radius;
    double zweight;
    u_int16_t num_levels;

    img_dmap_t *dmaps;    /* The pyramid of distance maps */
    u_int16_t w, h;       /* Width and height of the bottom of the pyramid */

} img_dist_pyr_t;

/* Create a gaussian pyramid from an image, with a number of extra levels */
img_pyr_t *img_create_gaussian_pyramid(img_t *img, int extra);

/* Create a pyramid of kd trees for the given image */
ANNkd_pyramid_t *img_create_ann_pyramid(img_pyr_t *pyr, error_metric_t em, int *nhood_radius);

/* Create a pyramid of distance maps for the given images */
img_dist_pyr_t *img_create_distance_pyramid(img_t *a, img_t *b, error_metric_t em, double eps);

/* Create a pyramid of distance maps for the given images */
img_dist_pyr_t *img_create_hierarchical_distance_pyramid(img_t *a, img_t *b, 
							 error_metric_t em, double eps, int min_size,
							 int rectified);

img_dist_pyr_t *img_create_hierarchical_distance_pyramid_2(img_t *a, img_t *b, 
							   error_metric_t em, double eps, int min_size,
							   int rectified);

/* Shrink-wrap an image distance pyramid, making each level as small
 * as possible */
img_dist_pyr_t *img_dist_pyr_shrink_wrap(img_dist_pyr_t *pyr);

/* Draw A using neighborhoods from B
 * a, b: input images
 * rgb: true if the result should be rgb, false if grayscale
 * t: interpolation parameter -- 0 = first image, 1 = second image */
img_t *img_pyramid_remap(img_dist_pyr_t *dist_pyr, int rgb, int rad);

/* Turn on verbose printing for pyramid creation */
void set_verbose_pyramid(int on);

#ifdef __cplusplus
}
#endif

#endif /* __pyramid_h__ */
