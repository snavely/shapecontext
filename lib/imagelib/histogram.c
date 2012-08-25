/* histogram.c */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "anniface.h"
#include "anntool.h"
#include "circle.h"
#include "color.h"
#include "histogram.h"
#include "image.h"
#include "util.h"

#define HIST_SIGMA 10.0

static double *erf_table = NULL;
static int global_bucket_size = 0;
static double integrate_gaussian(double mean, double sigma, double a, double b);

void precompute_erf_table(double sigma, int num_buckets) {
    int bucket_size, i;
    double a, b;

    /* The number of buckets must be a power of two plus 2 */
    if (!is_power_of_two(num_buckets - 2)) {
	printf("The number of buckets must be two more than a power of two.\n");
	return;
    }

    /* Prepare the table */
    bucket_size = 256 / (num_buckets - 2);
    global_bucket_size = bucket_size;
    
    if (erf_table)
	free(erf_table);

    erf_table = malloc(sizeof(double) * (3 * (bucket_size + 256) + 1));

    /* Calculate the integral for relative mean over all possible
     * intensities (i.e. any n / 3 in the range -bucket_size to 255 */
    
    a = -bucket_size;
    b = 0.0;
    for (i = -3 * bucket_size; i <= 3 * 256; i++) {
	erf_table[i + 3 * bucket_size] = 
	    integrate_gaussian(((double) i) / 3.0, sigma, a, b);
    }
}

static double erf_table_lookup(int r, int g, int b, int bucket) {
    int sum = r + g + b;
    int idx = sum - 3 * global_bucket_size * bucket;

    if (idx < 0)
	idx = -idx + 3 * global_bucket_size;

    return erf_table[idx];

#if 0
    check = integrate_gaussian((r + g + b) / 3.0, 10.0, 
			       bucket * global_bucket_size, (bucket + 1) * global_bucket_size);
    
    if (fabs(ret - check) > 0.0001) 
	printf("error!\n");
#endif
}

/* Return the area under the gaussian with the given mean and standard
 * deviation between points a and b */
static double integrate_gaussian(double mean, double sigma, double a, double b) {
    double aint = erf(((a - mean) / sigma) / M_SQRT2);
    double bint = erf(((b - mean) / sigma) / M_SQRT2);

    return 0.5 * (bint - aint);
}

#if 1
/* Compute the color histogram of a square neighborhood (with the
 * given length) of an image centered at the given pixel.  The results
 * will be stored in `vector'.  Note that num_buckets must be >= 3,
 * because there is an extra bucket on each side of the interval
 * [0,255]. */
void img_compute_histogram(img_t *img, int xmin, int ymin, int xmax, int ymax, int num_buckets, 
			   double *vector, double *weights, histogram_value_t htype, int normalize) {
    int i, x, y, idx = 0;
    // double bucket_size = global_bucket_size; // 256.0 / (num_buckets - 2);
    /* double sigma = 10.0; */
    
    double *cvec = NULL;
    double sum = 0.0;
    double avg = 0.0;
    double var = 0.0;
    double stddev = 0.0;
    double wt_sum = 0.0;
    int num_pixels = 0;

    /* Initialize the buckets */
    for (i = 0; i < num_buckets; i++)
	vector[i] = 0.0;

    if (normalize) {
	/* Compute the mean and variance of the samples in the
	 * window */

	cvec = (double *)malloc(sizeof(double) * (ymax - ymin + 1) * (xmax - xmin + 1));
	for (y = ymin; y <= ymax; y++) {
	    for (x = xmin; x <= xmax; x++) {
		color_t c = img_get_pixel(img, x, y);


		switch (htype) {
		case HIST_RED:
		    cvec[idx] = (double) c.r;
		    break;
		case HIST_GREEN:
		    cvec[idx] = (double) c.g;
		    break;
		case HIST_BLUE:
		    cvec[idx] = (double) c.b;
		    break;
		case HIST_INTENSITY:
		    cvec[idx] = color_intensity(c);
		    break;
		}

		sum += cvec[idx];
	    
		num_pixels++;
		idx++;
	    }
	}

	avg = sum / num_pixels;

	idx = 0;
	for (y = ymin; y <= ymax; y++) {
	    for (x = xmin; x <= xmax; x++) {
		var += (cvec[idx] - avg) * (cvec[idx] - avg);
		idx++;
	    }
	}

	var /= num_pixels;
	stddev = sqrt(var);
    }

    idx = 0;

    for (y = ymin; y <= ymax; y++) {
	for (x = xmin; x <= xmax; x++) {
	    double weight = weights[idx++];

	    /* Check if the pixel is defined */
	    if (!img_pixel_is_valid(img, x, y))
		continue;

	    wt_sum += weight;

#define EPSILON 1.0e-4
	    if (weight < EPSILON) {
		continue;
	    } else {
		/* Get the pixel color */
		color_t c = img_get_pixel(img, x, y);

#define DISTRIBUTION_WIDTH 2.0
		if (normalize) {
		    double val = 127.5 * ((cvec[idx] - avg) / (DISTRIBUTION_WIDTH * stddev)) + 127.5;
		    c.r = (unsigned char) clamp(val, 0.0, 255.0);
		    c.g = (unsigned char) clamp(val, 0.0, 255.0);
		    c.b = (unsigned char) clamp(val, 0.0, 255.0);
		}

		/* Add the area under the gaussian with mean intensity and
		 * standard deviation sigma squared */
		for (i = 0; i < num_buckets; i++) {
		    switch (htype) {
			case HIST_RED:
			    vector[i] += weight * erf_table_lookup(c.r, c.r, c.r, i - 1);
			    break;
			case HIST_GREEN:
			    vector[i] += weight * erf_table_lookup(c.g, c.g, c.g, i - 1);
			    break;
			case HIST_BLUE:
			    vector[i] += weight * erf_table_lookup(c.b, c.b, c.b, i - 1);
			    break;
			case HIST_INTENSITY:
			    vector[i] += weight * erf_table_lookup(c.r, c.g, c.b, i - 1);
			    break;
		    }

		    if (vector[i] < 0.0)
			printf("error!\n");
		}
	    }
	}
    }

    if (normalize) {
	free(cvec);
    }

    /* Divide by the total weight */
    for (i = 0; i < num_buckets; i++) {
	vector[i] /= wt_sum;
    }
}

/* Compute a sequence of n histograms inside the given region using an
 * array of weight arrays and storing the results in the array of
 * vector arrays */
void img_compute_histogram_seq(img_t *img, int n, int xmin, int ymin, int xmax, int ymax, int *num_buckets,
			       double **vector, double **weights, histogram_value_t htype) 
{
    int i;
    
    for (i = 0; i < n; i++)
	img_compute_histogram(img, xmin, ymin, xmax, ymax, num_buckets[i], vector[i], weights[i], htype, 0);
}

/* Count the number of buckets actually needed to store the (r,g,b) histogram
 * of the given image (before smoothing) */
int img_count_buckets(img_t *img, int num_buckets) {
    int *bucket_counts;
    int x, y, i;
    int bucket_size;
    int buckets_used = 0;

    if (!is_power_of_two(num_buckets)) {
	printf("Number of buckets must be a power of two.\n");
	return 0;
    }

    bucket_counts = (int *)calloc(num_buckets * num_buckets * num_buckets, sizeof(int));
    bucket_size = 256 / num_buckets;

    for (y = 0; y < img->h; y++) {
	for (x = 0; x < img->w; x++) {
	    color_t c = img_get_pixel(img, x, y);
	    int br = c.r / bucket_size;
	    int bg = c.g / bucket_size;
	    int bb = c.b / bucket_size;
	    int bidx = num_buckets * (num_buckets * br + bg) + bb;

	    bucket_counts[bidx]++;
	}
    }

    for (i = 0; i < num_buckets * num_buckets * num_buckets; i++) {
	if (bucket_counts[i] > 0) 
	    buckets_used++;
    }

    free(bucket_counts);

    return buckets_used;
}

/* Create a tree using histograms of r-neighborhoods surrounding each
 * pixel */
ANNkd_tree_t *img_create_kd_tree_histogram(img_t *img, double r, int num_buckets, int rgb) {
    ANNkd_tree_t *tree;
    int grid_length = circle_grid_side_length(r);
    int grid_size = circle_grid_size(r);
    int w = img->w, h = img->h, n = w * h; // n = (w - grid_length + 1) * (h - grid_length + 1);
    int x, y, i, idx = 0;
    int d = rgb ? (2 + 3 * num_buckets) : (2 + num_buckets);
    vec_t *pts = malloc(sizeof(vec_t) * n), axis_weights;

    double *buckets = malloc(sizeof(double) * (rgb ? (3 * num_buckets) : (num_buckets)));
    double *weights = malloc(sizeof(double) * grid_size);

    /* The grid should have odd length */
    assert(grid_length % 2 == 1);

    /* Fill in the weights */
    compute_circle_weights(r, r, weights);

    for (y = 0; y < h; y++) {
	for (x = 0; x < w; x++) {
    // for (y = grid_length / 2; y < h - grid_length / 2; y++) {
    //	for (x = grid_length / 2; x < w - grid_length / 2; x++) {
	    if (!rgb) {
		/* Compute the histogram of intensities */
		img_compute_histogram(img, 
				      x - grid_length / 2, y - grid_length / 2,
				      x + grid_length / 2, y + grid_length / 2, 
				      num_buckets, buckets, weights, HIST_INTENSITY, 0);

		pts[idx] = vec_new(d);

		Vn(pts[idx], 0) = x;
		Vn(pts[idx], 1) = y;

		/* Fill in the buckets */
		for (i = 0; i < num_buckets; i++)
		    Vn(pts[idx], i+2) = buckets[i];
		
		idx++;
	    } else {
		/* Compute the histogram for each channel separately */
		img_compute_histogram(img, 
				      x - grid_length / 2, y - grid_length / 2,
				      x + grid_length / 2, y + grid_length / 2, 
				      num_buckets, buckets, weights, HIST_RED, 0);

		img_compute_histogram(img, 
				      x - grid_length / 2, y - grid_length / 2,
				      x + grid_length / 2, y + grid_length / 2, 
				      num_buckets, buckets + num_buckets, weights, HIST_GREEN, 0);

		img_compute_histogram(img, 
				      x - grid_length / 2, y - grid_length / 2,
				      x + grid_length / 2, y + grid_length / 2, 
				      num_buckets, buckets + 2 * num_buckets, weights, HIST_BLUE, 0);

		pts[idx] = vec_new(d);

		Vn(pts[idx], 0) = x;
		Vn(pts[idx], 1) = y;

		/* Fill in the buckets */
		for (i = 0; i < 3 * num_buckets; i++)
		    Vn(pts[idx], i+2) = buckets[i];
		
		idx++;
	    }
	}
    }

    axis_weights = vec_new(d);

    Vn(axis_weights, 0) = 1.0;
    Vn(axis_weights, 1) = 1.0;

    if (!rgb) {
	for (i = 0; i < num_buckets; i++)
	    Vn(axis_weights, i+2) = get_ann_z_weight();
    } else {
	for (i = 0; i < 3 * num_buckets; i++)
	    Vn(axis_weights, i+2) = get_ann_z_weight();
    }

    tree = create_ann_tree(n, d, pts, axis_weights);

    /* Cleanup */
    for (i = 0; i < n; i++)
        vec_free(pts[i]);

    vec_free(axis_weights);
    free(pts);

    free(buckets);
    free(weights);

    return tree;
}

/* Create a tree using histograms of rings surrounding each
 * pixel */
ANNkd_tree_t *img_create_kd_tree_ring_histogram(img_t *img, double rmax, int num_histograms, int num_buckets, int rgb) {
    ANNkd_tree_t *tree;
    int grid_length = circle_grid_side_length(rmax);
    int grid_size = circle_grid_size(rmax);
    int w = img->w, h = img->h, n = w * h;
    int x, y, i, idx = 0;
    int d = rgb ? (2.0 + 3 + 3 * num_histograms * num_buckets) : (2.0 + 1 + num_histograms * num_buckets);
    vec_t *pts = malloc(sizeof(vec_t) * n), axis_weights;

    /* The first histogram is a circle of radius 2 -- the rest are rings */
    double ring_width = (rmax - 2.0) / (num_histograms - 1);

    /* We need buckets for each of num_histogram histograms */
    double *buckets = 
	malloc(sizeof(double) * (rgb ? (3 * num_histograms * num_buckets) : (num_histograms * num_buckets)));
    double *weights = malloc(sizeof(double) * num_histograms * grid_size);

    /* The grid should have odd length (why?) */
    assert(grid_length % 2 == 1);

    /* Fill in the weights.  The first histogram is circular and has
     * radius 2.0 */
    compute_circle_weights(rmax, 2.0, weights);

    /* The other histograms are rings */
    for (i = 1; i < num_histograms; i++)
	compute_ring_weights(rmax, 2.0 + i * ring_width, 2.0 + (i-1) * ring_width, weights + i * grid_size);

    for (y = 0; y < h; y++) {
	for (x = 0; x < w; x++) {
    // for (y = grid_length / 2; y < h - grid_length / 2; y++) {
    //	for (x = grid_length / 2; x < w - grid_length / 2; x++) {
	    if (!rgb) {
		/* The first bucket is just the pixel intensity */
		color_t c = img_get_pixel(img, x, y);
		double intensity = (double) (c.r + c.g + c.b) / 3.0;

		/* Compute the histograms of intensities */
		for (i = 0; i < num_histograms; i++) {
		    img_compute_histogram(img, 
					  x - grid_length / 2, y - grid_length / 2,
					  x + grid_length / 2, y + grid_length / 2, 
					  num_buckets, buckets + i * num_buckets, 
					  weights + i * grid_size, HIST_INTENSITY, 0);
		}
		
		pts[idx] = vec_new(d);

		Vn(pts[idx], 0) = x;
		Vn(pts[idx], 1) = y;
		Vn(pts[idx], 2) = intensity;

		/* Fill in the buckets */
		for (i = 0; i < num_histograms * num_buckets; i++)
		    Vn(pts[idx], i+3) = buckets[i];

		idx++;
	    } else {
		/* The first bucket is just the pixel color */
		color_t c = img_get_pixel(img, x, y);

		for (i = 0; i < num_histograms; i++) {
		    /* Compute the histogram for each channel separately */
		    img_compute_histogram(img, 
					  x - grid_length / 2, y - grid_length / 2,
					  x + grid_length / 2, y + grid_length / 2, 
					  num_buckets, buckets + 3 * i * num_buckets, 
					  weights + i * grid_size, HIST_RED, 0);

		    img_compute_histogram(img, 
					  x - grid_length / 2, y - grid_length / 2,
					  x + grid_length / 2, y + grid_length / 2, 
					  num_buckets, buckets + 3 * i * num_buckets + num_buckets, 
					  weights + i * grid_size, HIST_GREEN, 0);

		    img_compute_histogram(img, 
					  x - grid_length / 2, y - grid_length / 2,
					  x + grid_length / 2, y + grid_length / 2, 
					  num_buckets, buckets + 3 * i * num_buckets + 2 * num_buckets, 
					  weights + i * grid_size, HIST_BLUE, 0);
		}
		
		pts[idx] = vec_new(d);

		Vn(pts[idx], 0) = x;
		Vn(pts[idx], 1) = y;
		Vn(pts[idx], 2) = (double) c.r;
		Vn(pts[idx], 3) = (double) c.g;
		Vn(pts[idx], 4) = (double) c.b;

		/* Fill in the buckets */
		for (i = 0; i < 3 * num_histograms * num_buckets; i++)
		    Vn(pts[idx], i+5) = buckets[i];
		
		idx++;
	    }
	}
    }

    axis_weights = vec_new(d);

    Vn(axis_weights, 0) = 1.0;
    Vn(axis_weights, 1) = 1.0;

    if (rgb) {
	Vn(axis_weights, 2) = Vn(axis_weights, 3) = Vn(axis_weights, 4) = 0.001 * get_ann_z_weight();

	for (i = 5; i < d; i++)
	    Vn(axis_weights, i) = get_ann_z_weight();
    } else {
	Vn(axis_weights, 2) = 5.0;

	for (i = 3; i < d; i++)
	    Vn(axis_weights, i) = get_ann_z_weight();
    }

    tree = create_ann_tree(n, d, pts, axis_weights);

    /* Cleanup */
    for (i = 0; i < n; i++)
        vec_free(pts[i]);

    vec_free(axis_weights);
    free(pts);

    free(buckets);
    free(weights);

    return tree;
}

#endif

#if 0
int main() {
    double area1, area2, a = 64.0, b = 80.0;

    precompute_erf_table(10.0, 18);

    area1 = integrate_gaussian(67.0, 10.0, a, b);
    area2 = erf_table_lookup(67, 67, 67, 4);

    printf("Area from %0.3f to %0.3f is %0.3f\n", a, b, area1);
    printf("Area from %0.3f to %0.3f is %0.3f\n", a, b, area2);

    return 0;
}
#endif
