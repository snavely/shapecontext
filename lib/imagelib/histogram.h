/* histogram.h */
/* Compute image color histograms */

#ifndef __histogram_h__
#define __histogram_h__

#include "anniface.h"
#include "image.h"

typedef enum {
    HIST_RED,
    HIST_GREEN,
    HIST_BLUE,
    HIST_INTENSITY,
} histogram_value_t;

void precompute_erf_table(double sigma, int num_buckets);

/* Compute the color histogram of a square neighborhood (with the
 * given length) of an image centered at the given pixel.  The results
 * will be stored in `vector'.  Note that num_buckets must be >= 3,
 * because there is an extra bucket on each side of the interval
 * [0,255] */
void img_compute_histogram(img_t *img, int xmin, int ymin, int xmax, int ymax, int num_buckets, 
			   double *vector, double *weights, histogram_value_t htype, int normalize);

/* Create a tree using histograms of r-neighborhoods surrounding each
 * pixel */
ANNkd_tree_t *img_create_kd_tree_histogram(img_t *img, double r, int num_buckets, int rgb);

/* Create a tree using histograms of rings surrounding each
 * pixel */
ANNkd_tree_t *img_create_kd_tree_ring_histogram(img_t *img, double rmax, int num_histograms, int num_buckets, int rgb);

/* Count the number of buckets actually needed to store the (r,g,b) histogram
 * of the given image (before smoothing) */
int img_count_buckets(img_t *img, int num_buckets);

#endif /* __histogram_h__ */
