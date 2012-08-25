/* error.c */
/* Routines for computing various error metrics */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>

#include "anntool.h"
#include "circle.h"
#include "color.h"
#include "defines.h"
#include "error_metric.h"
#include "fileio.h"
#include "histogram.h"
#include "matrix.h"
#include "resample.h"
#include "util.h"

// #define EM_HGRAM_SIGMA 25.0
#define EM_HGRAM_SIGMA 15.0

/* Defines for histograms */
// #define BIG_CENTRAL_RING
// #define ADD_CENTER_PIXEL

/* Size of neighborhoods used in the EM_NHOOD_GS and EM_NHOOD_RGB
 * metrics */
static int em_nhood_size = 3;
static int em_nhood_w = 3;
static int em_nhood_h = 3;

/* Size of the radius of the neighborhood used in histogram
 * calculations */
static double em_hgram_nhood_radius = 4.0;

/* Number of buckets to store the histogram in */
static int em_num_hgram_buckets = 18;

/* The amount by which an image is oversampled, in order to increase
 * accuracy during the final rounds of fitting */ 
static int em_grid_fill_ratio = 1;

/* Array storing the current per-pixel histogram weights */
static double *em_histogram_weights = NULL;

/* Number of histograms to compute in the ring histogram method */
static int em_num_ring_hgrams = 4;

void em_set_nhood_size(int size) {
    if (size <= 0)
	printf("[em_set_nhood_size] The size of a neighborhood must be a positive number.\n");
    else
	em_nhood_size = size;
}

void em_set_nhood_dims(int w, int h) {
    if (w <= 0 || h <= 0) {
	printf("[em_set_nhood_dims] The dimensions of a neighborhood must be positive.\n");
    } else {
	em_nhood_w = w;
	em_nhood_h = h;
    }
}

void em_set_hgram_nhood_radius(double radius) { 
   if (radius <= 0.5 * M_SQRT2)
	printf("[em_set_hgram_nhood_radius] The radius of a histogram neighborhood must be at least sqrt(2).\n");
    else
	em_hgram_nhood_radius = radius;
}

void em_set_num_hgram_buckets(int num_buckets) {
    if (!is_power_of_two(num_buckets - 2))
	printf("[em_set_num_hgram_buckets] The number of buckets in a histogram must be 2 + a power of 2.\n");
    else
	em_num_hgram_buckets = num_buckets;
}

void em_set_grid_fill_ratio(int ratio) {
    if (ratio < 1)
	printf("[em_set_grid_fill_ratio] The grid fill ratio must be at least 1.\n");
    else
	em_grid_fill_ratio = ratio;
}

void em_set_num_ring_hgrams(int n) {
    if (n <= 0)
	printf("[em_set_num_ring_histograms] Number of histograms must be at least 1\n");
    else
	em_num_ring_hgrams = n;
}

static color_space_t em_color_space = CS_RGB;

void em_use_RGB() {
    em_color_space = CS_RGB;
}

void em_use_LUV() {
    em_color_space = CS_LUV;
}

color_space_t em_get_color_space() {
    return em_color_space;
}

/* Create a query point for the nearest-neighbor search for the
 * current error metric */
vec_t em_create_query_pt(img_t *img, error_metric_t em, double x, double y) {
    return em_create_query_pt_downsize(img, em, x, y, 1); /* No downsize */
}

/* Create a query point for the nearest-neighbor search for the
 * current error metric */
vec_t em_create_query_pt_downsize(img_t *img, error_metric_t em, double x, double y, int downsize) {
    vec_t pt;

    switch (em) {
	case EM_GRAYSCALE:
	case EM_GRAYSCALE_FILL: {
	    fcolor_t col = pixel_lerp(img, x, y);

	    pt = vec_new(3);
	    Vn(pt, 0) = (double) x + Vx(img->origin);
	    Vn(pt, 1) = (double) y + Vy(img->origin);
	    Vn(pt, 2) = (double) (col.r + col.g + col.b) / 3.0;

	    break;
	}
	case EM_RGB: {
	    fcolor_t col;

	    if (downsize == 1) {
		col = pixel_lerp(img, x, y);
	    } else {
		/* Average colors from a larger region */
		if (downsize != 2) {
		    printf("[em_create_query_pt_downsize] Can currently only downsize by two\n");
		    col = fcolor_new(0.0, 0.0, 0.0);
		} else {
		    col = img_average_2x2_square(img, (int) 2 * x, (int) 2 * y);
		}
	    }

	    pt = vec_new(5);
	    Vn(pt, 0) = (double) x + Vx(img->origin);
	    Vn(pt, 1) = (double) y + Vy(img->origin);

	    if (em_color_space == CS_RGB) {
		Vn(pt, 2) = (double) col.r;
		Vn(pt, 3) = (double) col.g;
		Vn(pt, 4) = (double) col.b;
	    } else {
		double L, U, V;

		color_RGBtoLUV(col.r, col.g, col.b, &L, &U, &V);
		Vn(pt, 2) = L;
		Vn(pt, 3) = U;
		Vn(pt, 4) = V;
	    }

	    break;
	}
	case EM_GRADIENT: {
	    fcolor_t col = pixel_lerp(img, x, y);

	    pt = vec_new(4);
	    Vn(pt, 0) = (double) x + Vx(img->origin);
	    Vn(pt, 1) = (double) y + Vy(img->origin);
	    Vn(pt, 2) = (double) col.r;
	    Vn(pt, 3) = img_gradient(img, x, y);

	    break;
	} 
	case EM_NHOOD_GS: {
	    int i, j;

	    pt = vec_new(2 + em_nhood_size * em_nhood_size);

	    Vn(pt, 0) = (double) x + Vx(img->origin);
	    Vn(pt, 1) = (double) y + Vy(img->origin);

	    for (i = 0; i < em_nhood_size; i++) {
		for (j = 0; j < em_nhood_size; j++) {
		    fcolor_t c = pixel_lerp(img, x + j - em_nhood_size / 2, y + i - em_nhood_size / 2);
		    Vn(pt, 2 + i * em_nhood_size + j) = fcolor_intensity(c); // (double) (c.r + c.g + c.b) / 3.0;
		}
	    }
	    break;
	}
	case EM_NHOOD_RGB: {
	    int i, j;

	    pt = vec_new(2 + 3 * em_nhood_size * em_nhood_size);

	    Vn(pt, 0) = (double) x + Vx(img->origin) / downsize;
	    Vn(pt, 1) = (double) y + Vy(img->origin) / downsize;

	    for (i = 0; i < em_nhood_size; i++) {
		for (j = 0; j < em_nhood_size; j++) {
		    fcolor_t c;

		    if (downsize == 1) {
			c = pixel_lerp(img, x + j - em_nhood_size / 2, y + i - em_nhood_size / 2);
		    } else if (downsize == 2) {
			c = img_average_2x2_square(img, 
						   (int) 2 * (x + j - em_nhood_size / 2), 
						   (int) 2 * (y + i - em_nhood_size / 2));
		    } else {
			printf("[em_create_query_pt] Error: downsize must be equal to 1 or 2\n");
		    }
		    
		    if (em_color_space == CS_RGB) {
			Vn(pt, 2 + (i * em_nhood_size + j) * 3 + 0) = (double) c.r;
			Vn(pt, 2 + (i * em_nhood_size + j) * 3 + 1) = (double) c.g;
			Vn(pt, 2 + (i * em_nhood_size + j) * 3 + 2) = (double) c.b;
		    } else if (em_color_space == CS_LUV) {
			double L, U, V;
			color_RGBtoLUV(c.r, c.g, c.b, &L, &U, &V);
			    
			Vn(pt, 2 + (i * em_nhood_size + j) * 3 + 0) = L;
			Vn(pt, 2 + (i * em_nhood_size + j) * 3 + 1) = U;
			Vn(pt, 2 + (i * em_nhood_size + j) * 3 + 2) = V;
		    } else {
			printf("[em_create_query_pt] Error: bad color space\n");
		    }
		}
	    }

	    break;
	}
	case EM_NHOOD_RGB_NORM: {
	    int i, j;
	    double r_sum = 0.0, g_sum = 0.0, b_sum = 0.0;
	    double r_avg = 0.0, g_avg = 0.0, b_avg = 0.0;
	    double r_var = 0.0, g_var = 0.0, b_var = 0.0;
	    double r_stddev = 0.0, g_stddev = 0.0, b_stddev = 0.0;
	    int num_pixels = 0;

	    pt = vec_new(2 + 3 * em_nhood_size * em_nhood_size);

	    Vn(pt, 0) = (double) x + Vx(img->origin);
	    Vn(pt, 1) = (double) y + Vy(img->origin);

	    for (i = 0; i < em_nhood_size; i++) {
		for (j = 0; j < em_nhood_size; j++) {
		    fcolor_t c = pixel_lerp(img, x + j - em_nhood_size / 2, y + i - em_nhood_size / 2);
	    
		    Vn(pt, 2 + (i * em_nhood_size + j) * 3 + 0) = (double) c.r;
		    Vn(pt, 2 + (i * em_nhood_size + j) * 3 + 1) = (double) c.g;
		    Vn(pt, 2 + (i * em_nhood_size + j) * 3 + 2) = (double) c.b;

		    r_sum += (double) c.r;
		    g_sum += (double) c.g;
		    b_sum += (double) c.b;

		    num_pixels++;
		}
	    }

	    r_avg = r_sum / num_pixels;
	    g_avg = g_sum / num_pixels;
	    b_avg = b_sum / num_pixels;

	    for (i = 0; i < em_nhood_size; i++) {
		for (j = 0; j < em_nhood_size; j++) {
		    double r = Vn(pt, 2 + (i * em_nhood_size + j) * 3 + 0);
		    double g = Vn(pt, 2 + (i * em_nhood_size + j) * 3 + 1);
		    double b = Vn(pt, 2 + (i * em_nhood_size + j) * 3 + 2);
	    
		    r_var += (r - r_avg) * (r - r_avg);
		    g_var += (g - g_avg) * (g - g_avg);
		    b_var += (b - b_avg) * (b - b_avg);
		}
	    }

	    r_stddev = sqrt(r_var);
	    g_stddev = sqrt(g_var);
	    b_stddev = sqrt(b_var);

	    if (r_stddev == 0.0)
		r_stddev = 1.0;

	    if (g_stddev == 0.0)
		g_stddev = 1.0;

	    if (b_stddev == 0.0)
		b_stddev = 1.0;

#if 0
	    for (i = 0; i < em_nhood_size; i++) {
		for (j = 0; j < em_nhood_size; j++) {
		    Vn(pt, 2 + (i * em_nhood_size + j) * 3 + 0) /= r_sum;
		    Vn(pt, 2 + (i * em_nhood_size + j) * 3 + 1) /= g_sum;
		    Vn(pt, 2 + (i * em_nhood_size + j) * 3 + 2) /= b_sum;
		}
	    }
#else
	    for (i = 0; i < em_nhood_size; i++) {
		for (j = 0; j < em_nhood_size; j++) {
		    double r = Vn(pt, 2 + (i * em_nhood_size + j) * 3 + 0);
		    double g = Vn(pt, 2 + (i * em_nhood_size + j) * 3 + 1);
		    double b = Vn(pt, 2 + (i * em_nhood_size + j) * 3 + 2);
	    
		    Vn(pt, 2 + (i * em_nhood_size + j) * 3 + 0) = (r - r_avg) / r_stddev;
		    Vn(pt, 2 + (i * em_nhood_size + j) * 3 + 1) = (g - g_avg) / g_stddev;
		    Vn(pt, 2 + (i * em_nhood_size + j) * 3 + 2) = (b - b_avg) / b_stddev;
		}
	    }
#endif

	    break;
	}
	case EM_HISTOGRAM_GS: {
	    double *buckets = malloc(sizeof(double) * em_num_hgram_buckets);
	    int grid_length = circle_grid_side_length(em_hgram_nhood_radius);
	    int i;

	    /* Initialize histogram weights if necessary */
	    if (em_histogram_weights == NULL) {
		int grid_size = circle_grid_size(em_hgram_nhood_radius);
		em_histogram_weights = malloc(sizeof(double) * grid_size);
		compute_circle_weights(em_hgram_nhood_radius, em_hgram_nhood_radius, em_histogram_weights);
	    }

	    img_compute_histogram(img, 
				  (int)x - grid_length / 2, (int)y - grid_length / 2,
				  (int)x + grid_length / 2, (int)y + grid_length / 2, 
				  em_num_hgram_buckets, buckets, em_histogram_weights,
				  HIST_INTENSITY, 0);

	    pt = vec_new(2 + em_num_hgram_buckets);
	    Vn(pt, 0) = (double) x + Vx(img->origin);
	    Vn(pt, 1) = (double) y + Vy(img->origin);

	    for (i = 0; i < em_num_hgram_buckets; i++)
		Vn(pt, i+2) = buckets[i];

	    free(buckets);

	    break;
	}
	case EM_HISTOGRAM_RGB: {
	    double *buckets = malloc(sizeof(double) * 3 * em_num_hgram_buckets);
	    int grid_length = circle_grid_side_length(em_hgram_nhood_radius);
	    int i;

	    /* Initialize histogram weights if necessary */
	    if (em_histogram_weights == NULL) {
		int grid_size = circle_grid_size(em_hgram_nhood_radius);
		em_histogram_weights = malloc(sizeof(double) * grid_size);
		compute_circle_weights(em_hgram_nhood_radius, em_hgram_nhood_radius, em_histogram_weights);
	    }

	    img_compute_histogram(img, 
				  (int) x - grid_length / 2, (int) y - grid_length / 2,
				  (int) x + grid_length / 2, (int) y + grid_length / 2, 
				  em_num_hgram_buckets, buckets, em_histogram_weights,
				  HIST_RED, 0);

	    img_compute_histogram(img, 
				  (int) x - grid_length / 2, (int) y - grid_length / 2,
				  (int) x + grid_length / 2, (int) y + grid_length / 2, 
				  em_num_hgram_buckets, 
				  buckets + em_num_hgram_buckets, 
				  em_histogram_weights,
				  HIST_GREEN, 0);

	    img_compute_histogram(img, 
				  (int) x - grid_length / 2, (int) y - grid_length / 2,
				  (int) x + grid_length / 2, (int) y + grid_length / 2, 
				  em_num_hgram_buckets, 
				  buckets + 2 * em_num_hgram_buckets, 
				  em_histogram_weights,
				  HIST_BLUE, 0);

	    pt = vec_new(2 + 3 * em_num_hgram_buckets);
	    Vn(pt, 0) = (double) x + Vx(img->origin);
	    Vn(pt, 1) = (double) y + Vy(img->origin);

	    for (i = 0; i < 3 * em_num_hgram_buckets; i++)
		Vn(pt, i+2) = buckets[i];

	    free(buckets);

	    break;
	}

	case EM_RING_HISTOGRAM_GS: {
	    double *buckets = malloc(sizeof(double) * em_num_ring_hgrams * em_num_hgram_buckets);
	    int grid_length = circle_grid_side_length(em_hgram_nhood_radius);
	    int grid_size = circle_grid_size(em_hgram_nhood_radius);
	    color_t c = img_get_pixel(img, x, y);
	    int i;

	    /* Initialize histogram weights if necessary */
	    if (em_histogram_weights == NULL) {
#ifdef BIG_CENTRAL_RING
		double ring_width = (em_hgram_nhood_radius - 2.0) / (em_num_ring_hgrams - 1);
#else
		double ring_width = (em_hgram_nhood_radius) / (em_num_ring_hgrams - 1);
#endif

		em_histogram_weights = malloc(sizeof(double) * em_num_ring_hgrams * grid_size);

		// #define USE_RINGS
#ifdef USE_RINGS
#ifdef BIG_CENTRAL_RING
		compute_circle_weights(em_hgram_nhood_radius, 2.0, em_histogram_weights);
		for (i = 1; i < em_num_ring_hgrams; i++)
		    compute_ring_weights(em_hgram_nhood_radius, 2.0 + i * ring_width, 2.0 + (i-1) * ring_width, 
					 em_histogram_weights + i * grid_size);
#else
		for (i = 0; i < em_num_ring_hgrams; i++)
		    compute_ring_weights(em_hgram_nhood_radius, (i+1) * ring_width, i * ring_width, 
					 em_histogram_weights + i * grid_size);
#endif /* BIG_CENTRAL_RING */
#else /* !USE_RINGS */
		/* Use gaussians */
		printf("ring_width = %0.3f\n", ring_width);
		for (i = 0; i < em_num_ring_hgrams; i++) {
		    compute_gaussian_ring_weights(em_hgram_nhood_radius, ((double) i) * ring_width,
						  0.8 * ring_width, em_histogram_weights + i * grid_size, 0);
		    
		    matrix_print(sqrt(grid_size), sqrt(grid_size), em_histogram_weights + i * grid_size);
		    printf("\n");
		}
#endif /* USE_RINGS */
	    }

	    for (i = 0; i < em_num_ring_hgrams; i++) {
		img_compute_histogram(img, 
				      (int) x - grid_length / 2, (int) y - grid_length / 2,
				      (int) x + grid_length / 2, (int) y + grid_length / 2, 
				      em_num_hgram_buckets, 
				      buckets + i * em_num_hgram_buckets, 
				      em_histogram_weights + i * grid_size,
				      HIST_INTENSITY, 1);
	    }
	    
#ifdef ADD_CENTER_PIXEL
	    pt = vec_new(2 + 1 + em_num_ring_hgrams * em_num_hgram_buckets);
#else
	    pt = vec_new(2 + em_num_ring_hgrams * em_num_hgram_buckets);
#endif

	    Vn(pt, 0) = (double) x + Vx(img->origin);
	    Vn(pt, 1) = (double) y + Vy(img->origin);

#ifdef ADD_CENTER_PIXEL
	    Vn(pt, 2) = (double) color_intensity(c);
	    for (i = 0; i < em_num_ring_hgrams * em_num_hgram_buckets; i++)
		Vn(pt, i+3) = buckets[i];
#else
	    for (i = 0; i < em_num_ring_hgrams * em_num_hgram_buckets; i++)
		Vn(pt, i+2) = buckets[i];
#endif	    


	    free(buckets);

	    break;
	}

	case EM_RING_HISTOGRAM_RGB: {
	    double *buckets = malloc(sizeof(double) * 3 * em_num_ring_hgrams * em_num_hgram_buckets);
	    int grid_length = circle_grid_side_length(em_hgram_nhood_radius);
	    int grid_size = circle_grid_size(em_hgram_nhood_radius);
	    color_t c = img_get_pixel(img, x, y);
	    int i;

	    /* Initialize histogram weights if necessary */
	    if (em_histogram_weights == NULL) {
#ifdef BIG_CENTRAL_RING
		double ring_width = (em_hgram_nhood_radius - 2.0) / (em_num_ring_hgrams - 1);
#else
		double ring_width = (em_hgram_nhood_radius) / (em_num_ring_hgrams);
#endif

		em_histogram_weights = malloc(sizeof(double) * em_num_ring_hgrams * grid_size);

#ifdef USE_RINGS
#ifdef BIG_CENTRAL_RING
		compute_circle_weights(em_hgram_nhood_radius, 2.0, em_histogram_weights);

		for (i = 1; i < em_num_ring_hgrams; i++)
		    compute_ring_weights(em_hgram_nhood_radius, 2.0 + i * ring_width, 2.0 + (i-1) * ring_width, 
					 em_histogram_weights + i * grid_size);
#else
		for (i = 0; i < em_num_ring_hgrams; i++)
		    compute_ring_weights(em_hgram_nhood_radius, (i+1) * ring_width, i * ring_width, 
					 em_histogram_weights + i * grid_size);
#endif /* BIG_CENTRAL_RING */
#else /* !USE_RINGS */

		for (i = 0; i < em_num_ring_hgrams; i++) {
		    compute_gaussian_ring_weights(em_hgram_nhood_radius, ((double) i) * ring_width,
						  0.8 * ring_width, em_histogram_weights + i * grid_size, 0);
		    
		    matrix_print(sqrt(grid_size), sqrt(grid_size), em_histogram_weights + i * grid_size);
		    printf("\n");
		}
#endif /* USE_RINGS */
	    }

	    for (i = 0; i < em_num_ring_hgrams; i++) {
		img_compute_histogram(img, 
				      (int) x - grid_length / 2, (int) y - grid_length / 2,
				      (int) x + grid_length / 2, (int) y + grid_length / 2, 
				      em_num_hgram_buckets, 
				      buckets + 3 * i * em_num_hgram_buckets, 
				      em_histogram_weights + i * grid_size,
				      HIST_RED, 1);

		img_compute_histogram(img, 
				      (int) x - grid_length / 2, (int) y - grid_length / 2,
				      (int) x + grid_length / 2, (int) y + grid_length / 2, 
				      em_num_hgram_buckets, 
				      buckets + 3 * i * em_num_hgram_buckets + em_num_hgram_buckets, 
				      em_histogram_weights + i * grid_size,
				      HIST_GREEN, 1);

		img_compute_histogram(img, 
				      (int) x - grid_length / 2, (int) y - grid_length / 2,
				      (int) x + grid_length / 2, (int) y + grid_length / 2, 
				      em_num_hgram_buckets, 
				      buckets + 3 * i * em_num_hgram_buckets + 2 * em_num_hgram_buckets, 
				      em_histogram_weights + i * grid_size,
				      HIST_BLUE, 1);
	    }
	    
#ifdef ADD_CENTER_PIXEL
	    pt = vec_new(2 + 3 + 3 * em_num_ring_hgrams * em_num_hgram_buckets);
#else
	    pt = vec_new(2 + 3 * em_num_ring_hgrams * em_num_hgram_buckets);
#endif

	    Vn(pt, 0) = (double) x + Vx(img->origin);
	    Vn(pt, 1) = (double) y + Vy(img->origin);

#ifdef ADD_CENTER_PIXEL
	    Vn(pt, 2) = (double) c.r;
	    Vn(pt, 3) = (double) c.g;
	    Vn(pt, 4) = (double) c.b;

	    for (i = 0; i < 3 * em_num_ring_hgrams * em_num_hgram_buckets; i++)
		Vn(pt, i+5) = buckets[i];

#else
	    for (i = 0; i < 3 * em_num_ring_hgrams * em_num_hgram_buckets; i++)
		Vn(pt, i+2) = buckets[i];
#endif

	    free(buckets);

	    break;
	}

    default:
	printf("[em_create_query_pt] Error: unknown error metric\n");
	break;
    }
    
    return pt;
}

int em_get_dimension(error_metric_t em) 
{    
    switch (em) {
    case EM_GRAYSCALE:
    case EM_GRAYSCALE_FILL:
	return 3;
	    
    case EM_RGB:
	return 5;

    case EM_GRADIENT:
	return 4;

    case EM_NHOOD_GS:
	return 2 + em_nhood_size * em_nhood_size;

    case EM_NHOOD_RGB:
    case EM_NHOOD_RGB_NORM:
	return 2 + 3 * em_nhood_size * em_nhood_size;
	
    case EM_HISTOGRAM_GS:
	return 2 + em_num_hgram_buckets;

    case EM_HISTOGRAM_RGB:
	return 2 + 3 * em_num_hgram_buckets;

    case EM_RING_HISTOGRAM_GS:
#ifdef ADD_CENTER_PIXEL
	return 2 + 1 + em_num_ring_hgrams * em_num_hgram_buckets;
#else
	return 2 + em_num_ring_hgrams * em_num_hgram_buckets;
#endif

    case EM_RING_HISTOGRAM_RGB:
#ifdef ADD_CENTER_PIXEL
	return 2 + 3 + 3 * em_num_ring_hgrams * em_num_hgram_buckets;
#else
	return 2 + 3 * em_num_ring_hgrams * em_num_hgram_buckets;
#endif

    default:
	printf("[em_get_dimension] Error: invalid metric\n");
	return 0;
    }    
}

double error_x_weight = 1.0;
double error_y_weight = 1.0;

void set_error_x_weight(double w) 
{
    if (w >= 0.0)
	error_x_weight = w;
}

double get_error_x_weight() 
{
    return error_x_weight;
}

void set_error_y_weight(double w) 
{
    if (w >= 0.0)
	error_y_weight = w;
}

double get_error_y_weight() 
{
    return error_y_weight;
}

vec_t em_create_axis_weights(error_metric_t em, int rectified) 
{
    int d = em_get_dimension(em);
    int i;

    vec_t axis_weights = vec_new(d);
    
    Vn(axis_weights, 0) = error_x_weight;

    if (rectified) {
	/* Weight the y-axis high */
	double other = 32.0 * (d - 2) * get_ann_z_weight() * get_ann_z_weight();
	Vn(axis_weights, 1) = sqrt(other);
    } else {
	Vn(axis_weights, 1) = error_y_weight;
    }
    
    /* Special case for ring histograms */
    if (em == EM_RING_HISTOGRAM_GS) {
#ifdef ADD_CENTER_PIXEL
	Vn(axis_weights, 2) = 0.001 * get_ann_z_weight();
	for (i = 3; i < d; i++)
	    Vn(axis_weights, i) = get_ann_z_weight();
#else
	for (i = 2; i < d; i++)
	    Vn(axis_weights, i) = get_ann_z_weight();
#endif
    } else if (em == EM_RING_HISTOGRAM_RGB) {
#ifdef ADD_CENTER_PIXEL
	Vn(axis_weights, 2) = Vn(axis_weights, 3) = Vn(axis_weights, 4) = 0.001 * get_ann_z_weight();
	for (i = 5; i < d; i++)
	    Vn(axis_weights, i) = get_ann_z_weight();
#else
	// Vn(axis_weights, 2) = Vn(axis_weights, 3) = Vn(axis_weights, 4) = 0.001 * get_ann_z_weight();
	for (i = 2; i < d; i++)
	    Vn(axis_weights, i) = get_ann_z_weight();
#endif
    } else if (em == EM_NHOOD_GS) {
	double *weights = (double *)malloc(sizeof(double) * em_nhood_size * em_nhood_size);
	compute_circle_weights(0.5 * ((double) em_nhood_size), 0.5 * ((double) em_nhood_size), weights);

	for (i = 2; i < d; i++)
	    Vn(axis_weights, i) = MAX(1.0e-4, weights[i] * get_ann_z_weight());

	free(weights);
    } else if (em == EM_NHOOD_RGB || em == EM_NHOOD_RGB_NORM) {
	double *weights = (double *)malloc(sizeof(double) * em_nhood_size * em_nhood_size);
	compute_circle_weights(0.5 * ((double) em_nhood_size), 0.5 * ((double) em_nhood_size), weights);

	for (i = 2; i < d; i++)
	    Vn(axis_weights, i) = MAX(1.0e-4, weights[(i-2) / 3] * get_ann_z_weight());

	free(weights);
    } else {
	for (i = 2; i < d; i++)
	    Vn(axis_weights, i) = get_ann_z_weight();
    }

    return axis_weights;
}

static int verbose_create_tree = 0;

/* Create a tree for part of an image based on the given error metric */
ANNkd_tree_t *img_create_ann_subtree(img_t *img, error_metric_t em, int *nhood_radius,
				     int xmin, int xmax, int ymin, int ymax, int rectified, int downsize,
				     ANNkd_tree_t *parent, int *idx_map) 
{
    ANNkd_tree_t *tree;
    clock_t start, end;
    int diam = 0, rad;
    int w, h, d, num_pts;
    double dx, dy;
    vec_t *pts, axis_weights;
    int i;

    if (idx_map && em_grid_fill_ratio != 1.0) {
	printf("[img_create_ann_subtree] Cannot handle em_grid_fill_ratio > 1.0 and idx_map != NULL\n");
	return NULL;
    }

    /* Compute the neighborhood radius */
    switch (em) {
	case EM_NHOOD_GS:
	    *nhood_radius = em_nhood_size / 2;
	    break;
	case EM_NHOOD_RGB:
	case EM_NHOOD_RGB_NORM:
	    *nhood_radius = em_nhood_size / 2;
	    break;
	case EM_HISTOGRAM_GS:
	    *nhood_radius = circle_grid_side_length(em_hgram_nhood_radius) / 2;
	    break;
	case EM_HISTOGRAM_RGB:
	    // *nhood_radius = circle_grid_side_length(em_hgram_nhood_radius) / 2;
	    *nhood_radius = 0;
	    break;
	case EM_RING_HISTOGRAM_GS:
	case EM_RING_HISTOGRAM_RGB:
	    // *nhood_radius = circle_grid_side_length(em_hgram_nhood_radius) / 2;
	    *nhood_radius = 0;
	    break;
	default:
	    *nhood_radius = 0;
	    break;
    }

    rad = *nhood_radius;
    diam = rad * 2 + 1;

    /* Metric-dependent setup */
    switch (em) {
	case EM_HISTOGRAM_GS:
	case EM_HISTOGRAM_RGB:
	case EM_RING_HISTOGRAM_GS:
	case EM_RING_HISTOGRAM_RGB:
	    precompute_erf_table(EM_HGRAM_SIGMA, em_num_hgram_buckets);
	    // *nhood_radius = circle_grid_side_length(em_hgram_nhood_radius) / 2;
	    *nhood_radius = 0;
	    break;

	default:
	    break;
    }

    start = clock();

    /* Build the kd-tree */
    w = xmax - xmin + 1;
    h = ymax - ymin + 1;
    num_pts = (em_grid_fill_ratio * (w - diam) + 1) * (em_grid_fill_ratio * (h - diam) + 1); //  / (downsize * downsize);

    d = em_get_dimension(em);

    dx = 1.0 / (double) em_grid_fill_ratio;
    dy = dx; 

    if (verbose_create_tree) {
	printf("[img_create_ann_subtree] Creating ANN tree:\n"
	       "[img_create_ann_subtree]     Points:      %d\n"
	       "[img_create_ann_subtree]     Dimensions:  %d\n"
	       "[img_create_ann_subtree]     Radius:      %d\n",
	       num_pts, d, *nhood_radius);
    }
    
    /* The computation above is not accurate when the image is too
     * small */
    if (w < diam || h < diam)
	num_pts = 0;

    if (num_pts == 0) {
	/* The image is too small -- create a tree with a single point */
	num_pts = 1;

	pts = malloc(sizeof(vec_t));
	pts[0] = vec_new(d);

	for (i = 0; i < d; i++)
	    Vn(pts[0], i) = 0.0;

	printf("[img_create_ann_subtree] Error: no points\n");
    } else {
	int idx = 0;
	double x, y;
	int rad2 = downsize * rad;

	pts = malloc(sizeof(vec_t) * num_pts);

	for (y = ymin + rad2; y <= ymin + h - 2 - rad2; y += dy) {
	    for (x = xmin + rad2; x <= xmin + w - 2 - rad2; x += dx) {
		/* Check if all the points in the neighborhood are
		 * valid */

		if (img_region_is_valid(img, x - rad2, x + rad2, y - rad2, y + rad2)) {
		    if (parent) {
			/* Compute the correct index into the parent
			 * tree */
			int parent_idx = idx_map[(int) y * img->w + (int) x];
			pts[idx++] = extract_ann_tree_pt(parent, parent_idx);
		    } else {
			pts[idx++] = em_create_query_pt_downsize(img, em, x / downsize, y / downsize, downsize);
			if (idx_map) {
			    idx_map[(int) y * img->w + (int) x] = idx - 1;
			}
		    }
		}
	    }
	}

	/* We may have created fewer points than we expected */
	if (idx > num_pts)
	    printf("[img_create_ann_subtree] Error: Created more points than expected.\n");

	num_pts = idx;

	if (num_pts == 0) {
	    printf("[img_create_ann_subtree] No points were created!\n");
	}
    }

    /* Create axis weights */
    axis_weights = em_create_axis_weights(em, rectified);

    /* Create the tree! */
    tree = create_ann_tree(num_pts, d, pts, axis_weights);

    /* Cleanup */
    for (i = 0; i < num_pts; i++)
	vec_free(pts[i]);

    vec_free(axis_weights);
    free(pts);

    end = clock();

    // printf("[Tick-Tock] kd-tree creation took %0.3fs\n", ((double) (end - start)) / CLOCKS_PER_SEC);

    return tree;
}

/* Create a tree for an image based on the given error metric */
ANNkd_tree_t *img_create_ann_tree(img_t *img, error_metric_t em, int *nhood_radius) {
    return img_create_ann_subtree(img, em, nhood_radius, 0, img->w - 1, 0, img->h - 1, 0, 1, NULL, NULL);
}

ANNkd_tree_t *img_create_ann_tree_params(img_t *img, error_metric_t em, int *nhood_radius, int rectified, int downsize) 
{
    return img_create_ann_subtree(img, em, nhood_radius, 0, img->w - 1, 0, img->h - 1, rectified, downsize, NULL, NULL);
}

ANNkd_tree_t *img_create_ann_tree_rectified(img_t *img, error_metric_t em, int *nhood_radius) {
    printf("rect\n");    
    return img_create_ann_subtree(img, em, nhood_radius, 0, img->w - 1, 0, img->h - 1, 1, 1, NULL, NULL);
}

/* Create a tree for an image based on the given error metric, and
 * make an index map along with it */
ANNkd_tree_t *img_create_ann_tree_idx_map(img_t *img, error_metric_t em, int *nhood_radius, int *idx_map) {
    return img_create_ann_subtree(img, em, nhood_radius, 0, img->w - 1, 0, img->h - 1, 0, 1, NULL, idx_map);
}

ANNkd_tree_t *img_create_ann_tree_idx_map_rectified(img_t *img, error_metric_t em, int *nhood_radius, int *idx_map) {
    printf("rect\n");    
    return img_create_ann_subtree(img, em, nhood_radius, 0, img->w - 1, 0, img->h - 1, 1, 1, NULL, idx_map);
}

/* Produces an image of the chamfer metric between two images */
img_t *img_error(img_t *a, img_t *b, error_metric_t em) {
    ANNkd_tree_t *tree;
    int rad;
    double *dist = NULL, max_dist = 0.0;
    int x, y, idx;
    vec_t sample;

    img_t *img = img_new(b->w, b->h);

    dist = malloc(sizeof(double) * b->w * b->h);
    tree = img_create_ann_tree(a, em, &rad);

    idx = 0;
    for (y = 0; y < b->h; y++) {
	for (x = 0; x < b->w; x++) {
	    sample = em_create_query_pt(b, em, x, y);
	    dist[idx] = query_ann_tree_color_dist(tree, sample, 0.0);

	    if (dist[idx] > max_dist)
		max_dist = dist[idx];

	    vec_free(sample);

	    idx++;
	}
    }

    /* Fill in the pixels */
    for (idx = 0; idx < b->h * b->w; idx++) {
	int c = (int) (256 * dist[idx] / max_dist);
	img->pixels[idx].r = img->pixels[idx].g = img->pixels[idx].b = c;
    }

    free(dist);
    free_ann_tree(tree);

    return img;
}

