/* alignicp.c */
/* Routines for aligning pairs of images using the Iterative Closest Point
 * algorithm */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>

#include "minpack.h"

#include "align.h"
#include "anniface.h"
#include "anntool.h"
#include "binimage.h"
#include "binvolume.h"
#include "circle.h"
#include "defines.h"
#include "dist.h"
#include "distutil.h"
#include "error.h"
#include "filter.h"
#include "histogram.h"
#include "horn.h"
#include "image.h"
#include "lerp.h"
#include "matrix.h"
#include "metric.h"
#include "qsort.h"
#include "ransac.h"
#include "resample.h"
#include "transform.h"
#include "transform3D.h"
#include "util.h"

/* If this is defined a bitmap will be written for each iteration
 * of ICP, showing how the alignment is doing */
#define ICP_OUTPUT_INTERMEDIATE

/* If this is defined the ICP algorithm will not consider matches 
 * between a point on the target image and one on the boundary of
 * the object image */
#define ICP_REJECT_BOUNDARY_MATCHES

/* Exactly one of the following should be defined.  If the first is
 * defined then all points in the image will make a contribution to
 * the ICP computation.  If the second is defined then
 * ICP_PAIRS_SAMPLE_PERCENT percent of points will be selected at
 * random */

// #define ICP_SELECT_ALL_PAIRS
#define ICP_SELECT_RANDOM_PAIRS

/* If this is defined the ICP algorithm will throw out all the worst
 * N% of matches when finding the optimal transform */
#define ICP_REJECT_WORST_N_PERCENT

/* These definitions specify the method by which distances between
 * points should be computed.  Only one of these should be defined.
 * Leave the first one defined if the distance transform should be
 * precomputed.  Leave the second defined if the approximate nearest
 * neighbor library should be used */
// #define ICP_USE_DISTANCE_TRANSFORM
#define ICP_USE_ANN

#define REPORT_TIMES

static double icp_reject_percent = 0.0;
static double icp_ann_epsilon = 0.0;
static double icp_pairs_sample_percent = 0.05;
static int icp_max_rounds = 64;

static error_metric_t icp_error_metric = EM_RGB;

static int icp_sample_intersection = 0;

/* Set the error metric functions used by icp */
int icp_set_error_metric(error_metric_t em) {
    icp_error_metric = em;
    return 0;
}

/* Set whether we just want to sample the intersection of two images */
void set_icp_sample_intersection(int val) 
{
    if (val != 0 && val != 1) {
	printf("[set_icp_sample_intersection] Invalid value\n");
	return;
    }
    
    icp_sample_intersection = val;
}

/* Set the percent of worst matches that will be rejected during ICP */
void set_icp_reject_percent(double p) {
    if (p < 0.0 || p > 1.0) {
        printf("[set_icp_reject_percent] Value must be in range [0.0, 1.0]\n");
        return;
    }

    icp_reject_percent = p;
}

/* Set the epsilon used in the nearest neighbor calculations */
void set_icp_ann_epsilon(double eps) {
    if (eps < 0.0) {
        printf("[set_icp_ann_epsilon] Epsilon must be non-negative\n");
        return;
    }
    
    icp_ann_epsilon = eps;
}

/* Set the percent of pairs which are randomly sampled during ICP */
void set_icp_pairs_sample_percent(double p) {
    if (p < 0.0 || p > 1.0) {
        printf("[set_icp_pairs_sample_percent] Value must be in range [0.0, 1.0]\n");
        return;
    }

    icp_pairs_sample_percent = p;
}

/* Set the maximum number of iterations for ICP to execute */
void set_icp_max_rounds(int max_rounds) {
    if (max_rounds <= 0) {
	printf("ICP must execute at least one round.\n");
	return;
    }
    
    icp_max_rounds = max_rounds;
}

#define MAX_TRIES 128
static int check_nhood = 1;

#define MAX_ROUNDS_SINCE_MIN 64
#define ERROR_FILTER_LENGTH 8

#define ICP_NHOOD_SIZE 5
#define ICP_GRID_FILL_RATIO 1
#define ICP_NHOOD_RADIUS 4.0
#define ICP_NUM_HISTOGRAM_BUCKETS 18

/* Find a linear transformation that aligns image B to image A.
 * Uses ICP and resamples the image after each iteration. */
trans2D_t *align_image_ICP(img_t *a, img_t *b, transform_class_t tclass) {
    trans2D_t *Tin = new_identity_transform();
    trans2D_t *Tout = align_Timage_ICP(a, b, Tin, tclass);
    transform_free(Tin);
    return Tout;
}

static int icp_recurse = 1;

void icp_recurse_off() {
    icp_recurse = 0;
}

void icp_recurse_on() {
    icp_recurse = 1;
}

/* Find a linear transformation that aligns image B to image A.
 * Uses ICP and resamples the image after each iteration.  The
 * transform that describes the initial guess is given. */
trans2D_t *align_Timage_ICP(img_t *a, img_t *b, trans2D_t *Tin, transform_class_t tclass) {
    int num_pts = a->w * a->h;

    int x, y, bit, max_bits, num_bits, i, j, k;
    trans2D_t *T = transform_copy(Tin); // new_identity_transform();
    trans2D_t *Tinv = transform_invert(Tin);
    trans2D_t *Tnew = new_identity_transform(), *Tprod;
    trans2D_t *Tout = new_identity_transform();
    trans2D_t *Tbest = transform_copy(Tin);
    trans2D_t *Tident = new_identity_transform();
    double error = DBL_MAX, min_error = DBL_MAX,
        olderror, errorx, errory, error_new;
    double scale;

    v3_t *r_pts = malloc(sizeof(v3_t) * num_pts);
    v3_t *l_pts = malloc(sizeof(v3_t) * num_pts);

    // double eps = 1.0e-4;
    double error_eps = 1.0e-2;
    double *errors = malloc(sizeof(double) * num_pts); /* Stores the errors 
                                                        * for each matching */
    vec_t *disps = malloc(sizeof(vec_t) * num_pts);    /* Stores displacements */
    vec_t *pts = malloc(sizeof(vec_t) * num_pts);      /* Stores points */
    double *grads = malloc(sizeof(double) * num_pts);  /* Stores gradients */
    int *perm = malloc(sizeof(int) * num_pts); /* Permutation after sorting 
                                                * errors */

    int num_samples;
    double sample_percent;
    double grad_max, grad_max_inv;
    double xtotal = 0.0, ytotal = 0.0;

    double error_hist[ERROR_FILTER_LENGTH];
    int error_idx = 0;

#ifdef REPORT_TIMES
    clock_t start, end;
#endif

    ANNkd_tree_t *tree;
    img_t *itmp;

    int round = 0, rounds_since_min = 0;

    int nhood_radius = 0;

    /* Set the error metric parameters */
    // em_set_nhood_size(ICP_NHOOD_SIZE);
    em_set_grid_fill_ratio(ICP_GRID_FILL_RATIO);
    em_set_hgram_nhood_radius(ICP_NHOOD_RADIUS);
    em_set_num_hgram_buckets(ICP_NUM_HISTOGRAM_BUCKETS);

    printf("creating tree...\n");

    tree = img_create_ann_tree(a, icp_error_metric, &nhood_radius);
    itmp = img_resample_bbox(b, T); // img_copy(b);

    printf("done\n");

#ifdef REPORT_TIMES
    start = clock();
#endif
    do {

	double max_error = 0.0;
	int failed_samples = 0;

#ifdef ICP_OUTPUT_INTERMEDIATE
	img_t *iout;
        bmp_t *bmp;
        FILE *f;
        char fname[32];
	double x0, y0, x1, y1, x2, y2, x3, y3;
#endif

        img_t *inext;

	check_nhood = 1;
	olderror = error;

        error = errorx = errory = 0.0;
        bit = 0;

#ifdef ICP_SELECT_ALL_PAIRS
        for (y = 0; y < itmp->h; y++) {
           for (x = 0; x < itmp->w; x++) {
#endif

#ifdef ICP_SELECT_RANDOM_PAIRS

#define MAX_SAMPLE_PERCENT 0.30
        /* Ramp up to a higher sample percent as we get closer to a
	 * solution */
	/* sample_percent = icp_pairs_sample_percent + 
	    (MAX_SAMPLE_PERCENT - icp_pairs_sample_percent) *
	    ((double) rounds_since_min / MAX_ROUNDS_SINCE_MIN); */
		
	sample_percent = icp_pairs_sample_percent;
	num_samples = ((double) itmp->h * itmp->w) * sample_percent;

        for (k = 0; k < num_samples; k++) {
	    x = rand() % itmp->w;
            y = rand() % itmp->h;
	    {
#endif

                double dx, dy, dz;
                  vec_t pt, nn, disp;


#if 0
		/* Make sure the neighborhood lies within the image */
		if (x - nhood_radius < 0 || y - nhood_radius < 0 ||
		    x + nhood_radius > b->w - 1 ||
		    y + nhood_radius > b->h - 1) {
		    k--;
		    continue;
		}

#define MAX_FAILED_SAMPLES 20000
		if (failed_samples > MAX_FAILED_SAMPLES) {
		    /* Image has probably been translated too far */
		    check_nhood = 0;
		    failed_samples = 0;
		}

		if (check_nhood && !nhood_image_within_image(b->w, b->h, b->origin, Tinv, 
							     x + Vx(itmp->origin), y + Vy(itmp->origin), nhood_radius)) {
		    k--;
		    failed_samples++;
		    continue;
		}

		failed_samples = 0;
#else 
		failed_samples = 0;

		do {
		    img_sample_random_pt(itmp, Tinv, b, nhood_radius, &x, &y, 0.0);
		    failed_samples++;
		} while (icp_sample_intersection && 
			 !img_region_is_valid(a, 
					      x + Vx(itmp->origin) - nhood_radius, x + Vx(itmp->origin) + nhood_radius, 
					      y + Vy(itmp->origin) - nhood_radius, y + Vy(itmp->origin) + nhood_radius));
		// !nhood_image_within_image(a->w, a->h, a->origin, Tident, x + Vx(itmp->origin), y + Vy(itmp->origin), nhood_radius) && failed_samples < 256);
#endif
                pt = em_create_query_pt(itmp, icp_error_metric, x, y);
		
		nn = query_ann_tree_pt(tree, pt, icp_ann_epsilon, NULL, 0, 1);
                disp = vec_sub(nn, pt);

		dx = Vx(disp);
                dy = Vy(disp);
		dz = Vz(disp);

                error += dx * dx + dy * dy; // vec_norm(disp);

		errors[bit] = vec_norm(disp) - dx * dx - dy * dy; // dx * dx + dy * dy; // vec_norm(disp);
                disps[bit] = disp;
                pts[bit] = pt;
		// grads[bit] = img_gradient(itmp, x, y);
		// grads[bit] = img_gradient_sobel(itmp, x, y, NULL);

		errorx += dx;
                errory += dy;
    
		if (errors[bit] > max_error)
		    max_error = errors[bit];

		r_pts[bit] = v3_new(x + dx, y + dy, 1);
		l_pts[bit] = v3_new(x, y, 1);

                bit++;

                vec_free(nn);
            }
        }

        max_bits = num_bits = bit;

        error /= num_bits;
        errorx /= num_bits;
        errory /= num_bits;

#ifdef ICP_REJECT_WORST_N_PERCENT
        if (icp_reject_percent > 0.0) {
            int pt_max;

            /* Sort the errors and create a permutation that
             * when applied to the original array results in the
             * sorted array */

	    // qsort_descending();
            // qsort_perm(num_bits, errors, perm);
	    
	    qsort_ascending();
	    qsort_perm(num_bits, grads, perm);
            permute(num_bits, sizeof(vec_t), disps, perm);
            permute(num_bits, sizeof(vec_t), pts, perm);
	    permute(num_bits, sizeof(double), grads, perm);

            pt_max = (int) num_bits * (1.0 - icp_reject_percent) + 1;
            if (pt_max > num_bits)
                printf("Error: too many points\n");

            num_bits = pt_max;
        }
#endif

	error_new = 0.0;
	errorx = errory = 0.0;

#if 0
	/* Compute the maximum gradient */
	grad_max = 0.0;

	for (bit = 0; bit < num_bits; bit++) {
	    if (grads[bit] > grad_max)
		grad_max = grads[bit];
	}

	grad_max_inv = 1.0 / grad_max;
	// grad_max_inv = 1.0;
#endif

        for (bit = 0; bit < num_bits; bit++) {
            double dx, dy, dz;
            double x, y, z;

	    error_new += Vx(disps[bit]) * Vx(disps[bit]) +
		Vy(disps[bit]) * Vy(disps[bit]);

	    errorx += Vx(disps[bit]);
	    errory += Vy(disps[bit]);

	    /* Scale the displacement by the magnitude of the gradient */
	    dx = Vx(disps[bit]);
            dy = Vy(disps[bit]);
            dz = Vz(disps[bit]);

            x = Vx(pts[bit]);
            y = Vy(pts[bit]);
            z = Vz(pts[bit]);

#if 0
	    grads[bit] *= grad_max_inv;
	    grads[bit] = pow(grads[bit], 0.1);
#else
	    grads[bit] = 1.0;
#endif

	    // printf("pt: (%d, %d), grad: %0.3f\n", (int) x, (int) y, grads[bit]);
        }

	error_new /= num_bits;
	errorx /= num_bits;
	errory /= num_bits;

	find_optimal_transform(tclass, num_bits, r_pts, l_pts, (double *) Tnew->T);

        Tprod = transform_product(Tnew, T);
        transform_free(T);
        T = Tprod;
        transform_free(Tinv);
        Tinv = transform_invert(T);

        print_transform(Tprod);

        /* Resample the image for the next iteration */
	printf("resampling...\n");

        inext = img_resample_bbox(b, T);
        img_free(itmp);
        itmp = inext;

	printf("done.\n");

#ifdef ICP_OUTPUT_INTERMEDIATE
        sprintf(fname, "icpout-%03d.bmp", round);

        f = fopen(fname, "w");
        if (f == NULL) {
            printf("Cannot open %s for writing\n", fname);
            return NULL;
        }

	iout = img_fuse(2, 1, a, itmp); // img_copy(itmp);

	transform_point(T, 0.0, 0.0, &x0, &y0);
	transform_point(T, 0.0, b->h - 1, &x1, &y1);
	transform_point(T, b->w - 1, b->h - 1, &x2, &y2);
	transform_point(T, b->w - 1, 0.0, &x3, &y3);

	img_draw_line(iout, x0, y0, x1, y1, 0x00, 0x00, 0x00);
	img_draw_line(iout, x1, y1, x2, y2, 0x00, 0x00, 0x00);
	img_draw_line(iout, x2, y2, x3, y3, 0x00, 0x00, 0x00);
	img_draw_line(iout, x3, y3, x0, y0, 0x00, 0x00, 0x00);

#if 1
	for (i = 0; i < num_bits; i++) {
	    // img_draw_pt(iout, Vx(pts[i]), Vy(pts[i]), 1, 255, 0, 0);
	    // img_draw_line(iout, Vx(pts[i]), Vy(pts[i]), Vx(pts[i]) + Vx(disps[i]), Vy(pts[i]) + Vy(disps[i]), 255, 0, 0);
	    // img_draw_pt(iout, Vx(pts[i]), Vy(pts[i]), 1, 0, 0, 255);

	    // printf("(%0.3f, %0.3f) ==> (%0.3f, %0.3f)\n", Vx(pts[i]), Vy(pts[i]), Vx(pts[i]) + Vx(disps[i]), Vy(pts[i]) + Vy(disps[i]));
	}
#endif

        bmp = img2bmp(iout);
        write_bmp(f, bmp);
        fclose(f);

        free_bmp(bmp);
	img_free(iout);
#endif     

        /* Free the points and displacement arrays */
        for (i = 0; i < max_bits; i++) {
            vec_free(pts[i]);
            vec_free(disps[i]);
        }
  
        /* Filter the error to reduce noise due to random sampling */
        error_hist[error_idx] = error;
	error = 0.0;
	
	for (j = 0; j < MIN(round + 1, ERROR_FILTER_LENGTH); j++)
	    error += error_hist[j];
	
	error /= MIN(round + 1, ERROR_FILTER_LENGTH);
        
	xtotal += errorx;
	ytotal += errory;

        printf("Round %3d: error == %0.5f, (unweighted = %0.5f, new = %0.5f), num_bits == %d\n", 
               round, error, error_hist[error_idx], error_new, num_bits);
        printf("Round %3d: average disp == (%0.3e, %0.3e), (%0.3f, %0.3f)\n", 
               round, errorx, errory, xtotal, ytotal);
	fflush(stdout);

        error_idx = (error_idx + 1) % ERROR_FILTER_LENGTH;

        if (round > ERROR_FILTER_LENGTH && error < min_error) {
            min_error = error;
            rounds_since_min = 0;
	    transform_free(Tbest);
	    Tbest = transform_copy(T);
        } else if (round <= ERROR_FILTER_LENGTH) {
	    transform_free(Tbest);
	    Tbest = transform_copy(T);
	}

        round++;

        /* Wait until the filter is full to start incrementing
         * the termination counter */
        if (round > ERROR_FILTER_LENGTH)
            rounds_since_min++;

    } while (round < icp_max_rounds && error > error_eps &&
             rounds_since_min < MAX_ROUNDS_SINCE_MIN);


#ifdef REPORT_TIMES
    end = clock();
    // printf("Alignment time: %0.3fs\n", (float)(end - start) / CLOCKS_PER_SEC);
    printf("Number of rounds: %d\n", round);
#endif

    transform_free(Tnew);
    transform_free(Tinv);

    free_ann_tree(tree);

    transform_free(T);
    T = Tbest;

    free(r_pts);
    free(l_pts);

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            Tout->T[i][j] = T->T[i][j];
        }
    }

    transform_free(T);
    transform_free(Tident);
    free(errors);
    free(disps);
    free(pts);
    free(grads);
    free(perm);

    if (0 && icp_recurse) {
	int old_max_rounds = icp_max_rounds;
	double old_z_weight = get_ann_z_weight();
	error_metric_t old_error_metric = icp_error_metric;

	set_ann_z_weight(0.5);
	icp_set_error_metric(EM_RGB);
	set_icp_max_rounds(35);

	icp_recurse = 0;

	/* Run a few rounds of ICP to correct small errors */
	Tnew = align_Timage_ICP(a, b, Tout, tclass);
	transform_free(Tout);
	Tout = Tnew;

	icp_recurse = 1;

	set_ann_z_weight(old_z_weight);
	icp_error_metric = old_error_metric;
	set_icp_max_rounds(old_max_rounds);
    }
    
    return Tout;
}

