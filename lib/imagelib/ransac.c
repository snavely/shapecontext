/* ransac.c */
/* Perform the RANSAC (RANdom SAmple Consensus) algorithm to align a
 * pair of images */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>

#include "affine.h"
#include "alignicp.h"
#include "anniface.h"
#include "anntool.h"
#include "bmp.h"
#include "error.h"
#include "histogram.h"
#include "homography.h"
#include "horn.h"
#include "image.h"
#include "ransac.h"
#include "resample.h"
#include "transform-opt.h"
#include "vector.h"

// #define RANSAC_OUTPUT_INTERMEDIATE

// #define RIGID_TRANSFORM
// #define AFFINE_TRANSFORM
#define HOMOGRAPHY_TRANSFORM

error_metric_t ransac_error_metric = EM_NHOOD_RGB;

#if defined (RIGID_TRANSFORM) || defined (AFFINE_TRANSFORM)
#define NUM_SAMPLE_POINTS 3
#else /* HOMOGRAPHY_TRANSFORM */
#define NUM_SAMPLE_POINTS 4
#endif

#define RANSAC_ANN_EPSILON 2.0

#define MAX_TRIES 256
static int sample_random_pt(int w, int h,             /* width, height of new image */
			    v2_t origin,              /* origin of new image */
			    trans2D_t *Tinv,          /* transform back to old image */
			    int w_old, int h_old,     /* w, h of old image */
			    v2_t origin_old,          /* origin of old image */
			    int nhood_radius, int *xout, int *yout) 
{
    int x, y;
    int tries = 0;
    int fox, foy;

    fox = (int) floor(Vx(origin));
    foy = (int) floor(Vy(origin));

    while (1) {
	x = fox + rand() % (w - 2 * nhood_radius) + nhood_radius;
	y = foy + rand() % (h - 2 * nhood_radius) + nhood_radius;

	/* Make sure the neighborhood lies within the original image */
	if (!nhood_image_within_image(w_old, h_old, origin_old, Tinv, x, y, nhood_radius)) {
	    tries++;

	    if (tries > MAX_TRIES)
		return 1;

	    continue;
	}

	break;
    }

    *xout = x;
    *yout = y;

    return 0;
}

#define MAX_ROUNDS 1
#define RANSAC_SAMPLE_RATE 0.10

static int ransac_max_trials = 1024;
static int verbose_ransac = 1;

/* Set the maximum number of trials that RANSAC will use */
void set_ransac_max_trials(int n) 
{
    if (n <= 0) {
	printf("[set_ransac_max_trials] Maximum number of trials must be positive\n");
	return;
    }
    
    ransac_max_trials = n;
}

trans2D_t *align_image_ransac(img_t *a, img_t *b, int use_icp, double *error_out) {
    trans2D_t *Tin = new_identity_transform();
    trans2D_t *Tout = align_Timage_ransac(a, b, Tin, use_icp, error_out);

    transform_free(Tin);
    return Tout;
}

trans2D_t *align_image_ransac_map(img_t *a, img_t *b, int diameter, 
				  img_dmap_t *amap, img_dmap_t *bmap, int use_icp, double *error_out) 
{
    trans2D_t *Tin = new_identity_transform();
    trans2D_t *Tout = align_Timage_ransac_map(a, b, Tin, diameter, amap, bmap, use_icp, error_out);

    transform_free(Tin);
    return Tout;
}

/* **** NOTE:  l_pts come from b, r_pts come from a **** */

trans2D_t *align_Timage_ransac(img_t *a, img_t *b, trans2D_t *Tin, int use_icp, double *error_out) {
    static int num_ransac_calls = 0;

    int round, trial;
    ANNkd_tree_t *tree_small = NULL, *tree_large = NULL;
    vec_t sample, nn, *samples;
    v3_t *r_pts = NULL, *l_pts = NULL, *r_pts_best = NULL, *l_pts_best = NULL;
    int nhood_radius_large, nhood_radius_small;
    int w = b->w, h = b->h, i, j, x, y;
    // double trans[3][3], rot[3][3];
    trans2D_t *Tnew = NULL, *T = NULL, *Tinv = NULL, *Tprod = NULL;
    img_t *Timg = NULL, *itmp = NULL, *a_cpy, *b_cpy;
    double *weights = NULL;
    int num_test_points = w * h * RANSAC_SAMPLE_RATE;
    double errors[NUM_SAMPLE_POINTS];
    double min_dist = DBL_MAX;

    bmp_t *bmp;
    FILE *f;
    char fname[128];
    img_t *Timg_test;

    clock_t start, end;
    unsigned long sample_time = 0, error_time = 0, correct_time = 0;

    start = clock();

    T = transform_copy(Tin); // new_identity_transform();
    Tinv = transform_invert(T); // new_identity_transform();
    Tnew = new_identity_transform();

    samples = malloc(sizeof(double) * NUM_SAMPLE_POINTS);
    weights = malloc(sizeof(double) * NUM_SAMPLE_POINTS);
    r_pts = malloc(sizeof(v3_t) * NUM_SAMPLE_POINTS);
    l_pts = malloc(sizeof(v3_t) * NUM_SAMPLE_POINTS);
    r_pts_best = malloc(sizeof(v3_t) * NUM_SAMPLE_POINTS);
    l_pts_best = malloc(sizeof(v3_t) * NUM_SAMPLE_POINTS);

    for (i = 0; i < NUM_SAMPLE_POINTS; i++)
	weights[i] = 1.0;

    em_set_num_hgram_buckets(18);
    em_set_hgram_nhood_radius(4.0);
     // em_set_grid_fill_ratio(4);

    /* Create the nearest-neighbor search structure for the first
     * image */
    set_ann_z_weight(1000.0);
    em_set_nhood_size(5);
    tree_large = img_create_ann_tree(a, ransac_error_metric, &nhood_radius_large);
    // set_ann_z_weight(0.5);
    set_ann_z_weight(0.1);

    em_set_nhood_size(3);
    tree_small = img_create_ann_tree(a, EM_RGB, &nhood_radius_small);

    printf("[ransac] Done creating tree\n");

    Timg = img_resample_bbox(b, T); // img_copy(b);

    for (round = 0; round < MAX_ROUNDS; round++) {
	double min_error = DBL_MAX;
	trans2D_t *Tbest = NULL;

	for (trial = 0; trial < ransac_max_trials; trial++) {
	    // img_t *Timg_test;     /* The temporary transformed image */
	    trans2D_t *Tinv_test; /* The temporary inverse transform */
	    double error = 0.0;
	    double min_sample_dist = DBL_MAX;
	    // v2_t crs[4], min, max;
	    // int w_new, h_new;

	    /* Sample four (three?) points at random from the second
	     * image, and compute their nearest neighbors in the first
	     * image */

	    em_set_nhood_size(5);

	    start = clock();
	    for (i = 0; i < NUM_SAMPLE_POINTS; i++) {
		int reselect;
		do {
		    reselect = 0;
		    img_sample_random_pt(Timg, Tinv, b, nhood_radius_large, &x, &y, 0.0);

		    /* Make sure we didn't pick this point already */
		    for (j = 0; j < i; j++) {
			if (x == Vx(l_pts[j]) && y == Vy(l_pts[j])) {
			    reselect = 1;
			    break;
			}
		    }

		    if (reselect)
			continue;
		    
		    sample = em_create_query_pt(Timg, ransac_error_metric, x, y);

		    /* Check that the query point is not too uniform
		     * in color */

		    min_sample_dist = 0.0;

		    if (reselect) {
			vec_free(sample);
			continue;
		    }

		    nn = query_ann_tree_pt(tree_large, sample, RANSAC_ANN_EPSILON, errors + i, 0, 1);
		    errors[i] -= (Vx(sample) - Vx(nn)) * (Vx(sample) - Vx(nn));
		    errors[i] -= (Vy(sample) - Vy(nn)) * (Vy(sample) - Vy(nn));

		    /* Check that two points don't map to the same
		     * nearest neighbor */
		    for (j = 0; j < i; j++) {
			if (Vx(nn) == Vx(r_pts[j]) && Vy(nn) == Vy(r_pts[j])) {
			    reselect = 1;
			    vec_free(sample);
			    vec_free(nn);
			    break;
			}
		    }
		} while (reselect);

		if (errors[i] < min_dist) {
		    min_dist = errors[i];
		}
		
		r_pts[i] = v3_new(Vx(nn), Vy(nn), 0.0);
		l_pts[i] = v3_new(Vx(sample), Vy(sample), 0.0);

		samples[i] = sample;
		vec_free(nn);
	    }
	    end = clock();
	    sample_time += end - start;

	    /* Free all the samples */
	    for (i = 0; i < NUM_SAMPLE_POINTS; i++) 
		vec_free(samples[i]);

#ifdef RIGID_TRANSFORM
	    /* Compute the optimal translation and rotation using the Horn
	     * solver */ 
	    align_horn(NUM_SAMPLE_POINTS, r_pts, l_pts, (double *)rot, (double *)trans, (double *)Tnew->T, weights);
#endif

#ifdef AFFINE_TRANSFORM
	    align_affine(NUM_SAMPLE_POINTS, r_pts, l_pts, (double *)Tnew->T);
#endif

#ifdef HOMOGRAPHY_TRANSFORM
	    align_homography(NUM_SAMPLE_POINTS, r_pts, l_pts, (double *)Tnew->T, 0);
#endif


#if 0
	    for (j = 0; j < NUM_SAMPLE_POINTS; j++)
		printf("(%0.3f, %0.3f) ==> (%0.3f, %0.3f) [%0.3f, %0.3f] %0.3f\n",
		       Vx(l_pts[j]), Vy(l_pts[j]), Vx(r_pts[j]), Vy(r_pts[j]),
		       Vx(r_pts[j]) - Vx(l_pts[j]), Vy(r_pts[j]) - Vy(l_pts[j]),
		       errors[j]);

	    printf("%0.3f %0.3f\n", error, min_sample_dist);
#endif

	    if (isnan(Tnew->T[0][0])) {
		printf("result is nan\n");
		trial--;
		continue;
	    }

	    Tprod = transform_product(Tnew, T);
	    Tinv_test = transform_invert(Tprod);

	    /* Select points inside the transformed image */

	    error = 0.0;
	    sample = vec_new(5);

	    start = clock();
	    for (j = 0; j < num_test_points; j++) {
		double dist;
		fcolor_t col;

		if (sample_random_pt(w, h, b->origin, 
				     Tinv_test, 
				     w, h, b->origin,
				     nhood_radius_small, &x, &y) == 1) 
		    {
			error = DBL_MAX;
			break;
		    }

		/* Compute the pixel at (x, y) in the transformed
		 * image */

		col = pixel_transform(b, Tinv_test, x, y);

		Vx(sample) = x;
		Vy(sample) = y;
		Vn(sample, 2) = (float) col.r;
		Vn(sample, 3) = (float) col.g;
		Vn(sample, 4) = (float) col.b;

		// sample = em_transform_create_query_pt(Timg, Tprod, EM_RGB, x, y);
		dist = query_ann_tree_dist(tree_small, sample, 2.0, NULL);

		error += dist;
	    }
	    end = clock();
	    error_time += end - start;


	    vec_free(sample);

	    /* Check if this is the new best transform */
	    if (error < min_error) {
		if (Tbest != NULL) 
		    transform_free(Tbest);
		Tbest = transform_copy(Tprod);
		min_error = error;

		memcpy(r_pts_best, r_pts, sizeof(v3_t) * NUM_SAMPLE_POINTS);
		memcpy(l_pts_best, l_pts, sizeof(v3_t) * NUM_SAMPLE_POINTS);
	    }

	    transform_free(Tprod);
	    transform_free(Tinv_test);

	    if (min_error < 100.0) {
		printf("exiting after %d trials\n", trial);
		break;
	    }

#ifdef RANSAC_OUTPUT_INTERMEDIATE
	    Timg_test = img_resample_bbox(b, Tprod);

	    sprintf(fname, "ransac_%03d-%03d.bmp", num_ransac_calls, trial);

	    f = fopen(fname, "w");
	    if (f == NULL) {
		printf("Cannot open %s for writing\n", fname);
		return NULL;
	    }

	    bmp = img2bmp(Timg_test);
	    write_bmp(f, bmp);
	    fclose(f);

	    free_bmp(bmp);
	    img_free(Timg_test);
#endif
	}

	printf("Minimum error: %0.3f\n", min_error);
	printf("Minimum dist: %0.3f\n", min_dist);

	/* Transform the image by Tbest for the next round */
	itmp = img_resample(b, Tbest);
	img_free(Timg);
	Timg = itmp;

	transform_free(T);
	T = Tbest;

	transform_free(Tinv);
	Tinv = transform_invert(T);

	if (error_out != NULL)
	    *error_out = min_error;
    }
    
    free(r_pts);
    free(l_pts);
    free(weights);
    free(samples);

    transform_free(Tinv);
    transform_free(Tnew);

    free_ann_tree(tree_small);
    free_ann_tree(tree_large);

    print_transform(T);

    Timg_test = img_resample_bbox(b, T);
    a_cpy = img_copy(a);
    b_cpy = img_resample_bbox(b, Tin);

    printf("Matches:\n");

    Tinv = transform_invert(Tin);

    for (i = 0; i < NUM_SAMPLE_POINTS; i++) {
	int r = ((i % 2) == 0 || (i % 3) == 0) ? 255 : 0;
	int g = ((i % 2) == 1 || (i % 3) == 1) ? 255 : 0;
	int b = (i % 3) == 2 ? 255 : 0;
	double xout, yout;

	printf("        (%0.3f, %0.3f) ==> (%0.3f, %0.3f) [%0.3f, %0.3f]\n", 
	       Vx(l_pts_best[i]), Vy(l_pts_best[i]),
	       Vx(r_pts_best[i]), Vy(r_pts_best[i]),
	       Vx(r_pts_best[i]) - Vx(l_pts_best[i]), Vy(r_pts_best[i]) - Vy(l_pts_best[i]));

	img_draw_pt(a_cpy, Vx(r_pts_best[i]), Vy(r_pts_best[i]), 5, r, g, b);
	img_draw_pt(b_cpy, Vx(l_pts_best[i]), Vy(l_pts_best[i]), 5, r, g, b);

	transform_point(Tinv, Vx(l_pts_best[i]), Vy(l_pts_best[i]), &xout, &yout);
	transform_point(T, xout, yout, &xout, &yout);
	img_draw_pt(Timg_test, (int) rint(xout), (int) rint(yout), 5, r, g, b);
    }

    transform_free(Tinv);

    sprintf(fname, "ransac_post-%03d.bmp", num_ransac_calls);
    f = fopen(fname, "w");
    if (f == NULL) {
	printf("Cannot open %s for writing\n", fname);
	return NULL;
    }

    bmp = img2bmp(Timg_test);
    write_bmp(f, bmp);
    fclose(f);

    free_bmp(bmp);
    img_free(Timg_test);

    sprintf(fname, "ransac_pre0-%03d.bmp", num_ransac_calls);
    f = fopen(fname, "w");
    if (f == NULL) {
	printf("Cannot open %s for writing\n", fname);
	return NULL;
    }

    bmp = img2bmp(a_cpy);
    write_bmp(f, bmp);
    fclose(f);

    free_bmp(bmp);
    img_free(a_cpy);

    sprintf(fname, "ransac_pre1-%03d.bmp", num_ransac_calls);
    f = fopen(fname, "w");
    if (f == NULL) {
	printf("Cannot open %s for writing\n", fname);
	return NULL;
    }

    bmp = img2bmp(b_cpy);
    write_bmp(f, bmp);
    fclose(f);

    free_bmp(bmp);
    img_free(b_cpy);

    free(r_pts_best);
    free(l_pts_best);

    img_free(Timg);

    num_ransac_calls++;

    start = clock();

    if (use_icp) {
	icp_recurse_off();
	// set_ann_z_weight(0.1); /* For nhood */
	set_ann_z_weight(0.01); /* For nhood */
	// set_ann_z_weight(20.0); /* For ring */
	set_icp_sample_intersection(1);
	// set_icp_reject_percent(0.80);
	set_icp_pairs_sample_percent(0.05);
	set_icp_ann_epsilon(0.0);
	em_set_nhood_size(1);
	icp_set_error_metric(EM_NHOOD_RGB);
	set_icp_max_rounds(64);

	/* Run a few rounds of ICP to correct small errors */
	Tnew = align_Timage_ICP(a, b, T, TRANSFORM_HOMOGRAPHY);
	transform_free(T);

	end = clock();
	correct_time = end - start;
    } else {
	correct_time = 0.0;
	Tnew = transform_copy(T);
	transform_free(T);
    }

    printf("[Tick-Tock] Sample time:      %0.3fs\n"
	   "[Tick-Tock] Error time:       %0.3fs\n"
	   "[Tick-Tock] Correction time:  %0.3fs\n",
	   (double) sample_time / CLOCKS_PER_SEC,
	   (double) error_time / CLOCKS_PER_SEC,
	   (double) correct_time / CLOCKS_PER_SEC);

    print_transform(Tnew);

    return Tnew;
}




/* **** NOTE:  l_pts come from b, r_pts come from a **** */
/* This is the same function as above but takes a distance map for a
 * and Tin(b) as input */
trans2D_t *align_Timage_ransac_map(img_t *a, img_t *b, trans2D_t *Tin, int diameter, 
				   img_dmap_t *amap, img_dmap_t *bmap, int use_icp, double *error_out)
{
    static int num_ransac_calls = 0;

    int round, trial;
    ANNkd_tree_t *tree_small = NULL;
    vec_t sample, nn, *samples;
    v3_t *r_pts = NULL, *l_pts = NULL, *r_pts_best = NULL, *l_pts_best = NULL;
    int nhood_radius_large = diameter / 2, nhood_radius_small;
    int w = b->w, h = b->h, i, j, x, y;
    // double trans[3][3], rot[3][3];
    trans2D_t *Tnew = NULL, *T = NULL, *Tinv = NULL, *Tprod = NULL;
    img_t *Timg = NULL, *itmp = NULL, *a_cpy, *b_cpy, *img_test;
    double *weights = NULL;
    int num_test_points = w * h * RANSAC_SAMPLE_RATE;
    double errors[NUM_SAMPLE_POINTS];
    double min_dist = DBL_MAX;

    bmp_t *bmp;
    FILE *f;
    char fname[128];
    img_t *Timg_test, *img_fused;
    int disconnected = 0;

    clock_t start, end;
    unsigned long sample_time = 0, error_time = 0, correct_time = 0;

    start = clock();

    T = transform_copy(Tin); // new_identity_transform();
    Tinv = transform_invert(T); // new_identity_transform();
    Tnew = new_identity_transform();

    samples = malloc(sizeof(double) * NUM_SAMPLE_POINTS);
    weights = malloc(sizeof(double) * NUM_SAMPLE_POINTS);
    r_pts = malloc(sizeof(v3_t) * NUM_SAMPLE_POINTS);
    l_pts = malloc(sizeof(v3_t) * NUM_SAMPLE_POINTS);
    r_pts_best = malloc(sizeof(v3_t) * NUM_SAMPLE_POINTS);
    l_pts_best = malloc(sizeof(v3_t) * NUM_SAMPLE_POINTS);

    for (i = 0; i < NUM_SAMPLE_POINTS; i++)
	weights[i] = 1.0;

#if 0 /* Don't create the big tree anymore */
    em_set_num_hgram_buckets(18);
    em_set_hgram_nhood_radius(4.0);
     // em_set_grid_fill_ratio(4);

    /* Create the nearest-neighbor search structure for the first
     * image */
    set_ann_z_weight(1000.0);
    em_set_nhood_size(5);
    tree_large = img_create_ann_tree(a, ransac_error_metric, &nhood_radius_large);
#endif

    set_ann_z_weight(0.5);
    em_set_nhood_size(3);
    tree_small = img_create_ann_tree(a, EM_RGB, &nhood_radius_small);

    if (verbose_ransac)
	printf("[ransac] Done creating tree\n");

    Timg = img_resample_bbox(b, T); // img_copy(b);

    for (round = 0; round < MAX_ROUNDS; round++) {
	double min_error = DBL_MAX;
	trans2D_t *Tbest = NULL;

	for (trial = 0; trial < ransac_max_trials; trial++) {
	    // img_t *Timg_test;     /* The temporary transformed image */
	    trans2D_t *Tinv_test; /* The temporary inverse transform */
	    double error = 0.0;
	    double min_sample_dist = DBL_MAX;

	    int w_new, h_new;
	    v2_t origin_new;

	    int bad = 0;

	    /* Sample four (three?) points at random from the second
	     * image, and compute their nearest neighbors in the first
	     * image */

	    em_set_nhood_size(5);

	    start = clock();
	    for (i = 0; i < NUM_SAMPLE_POINTS; i++) {
		int reselect;
		do {
		    double x_nn, y_nn, x_nn_back, y_nn_back, back_disp;

		    reselect = 0;
		    img_sample_random_pt(Timg, Tinv, b, nhood_radius_large, &x, &y, 0.0);

		    /* Make sure we didn't pick this point already */
		    for (j = 0; j < i; j++) {
			if (x == Vx(l_pts[j]) && y == Vy(l_pts[j])) {
			    reselect = 1;
			    break;
			}
		    }

		    if (reselect)
			continue;
		    
		    sample = em_create_query_pt(Timg, ransac_error_metric, x, y);

		    /* Check that the query point is not too uniform
		     * in color */

		    min_sample_dist = 0.0;

		    x_nn = Vx(bmap->nns[y * bmap->w + x]);
		    y_nn = Vy(bmap->nns[y * bmap->w + x]);

		    /* Make sure that the relation is (more or less)
		     * symmetric */

		    x_nn_back = Vx(amap->nns[(int)(rint(y_nn) * amap->w + rint(x_nn))]);
		    y_nn_back = Vy(amap->nns[(int)(rint(y_nn) * amap->w + rint(x_nn))]);

		    back_disp = (x - x_nn_back) * (x - x_nn_back) + (y - y_nn_back) * (y - y_nn_back);

#define SIMILARITY_THRESHOLD (4.0 * 4.0)
		    if (back_disp > SIMILARITY_THRESHOLD) {
			reselect = 1;
			vec_free(sample);
			continue;
		    }

		    nn = em_create_query_pt(a, ransac_error_metric, x_nn, y_nn); // query_ann_tree_pt(tree_large, sample, RANSAC_ANN_EPSILON, errors + i);
		    errors[i] = bmap->dists[y * bmap->w + x];
		    errors[i] -= (Vx(sample) - Vx(nn)) * (Vx(sample) - Vx(nn));
		    errors[i] -= (Vy(sample) - Vy(nn)) * (Vy(sample) - Vy(nn));

		    /* Check that two points don't map to the same
		     * nearest neighbor */
		    for (j = 0; j < i; j++) {
			if (Vx(nn) == Vx(r_pts[j]) && Vy(nn) == Vy(r_pts[j])) {
			    reselect = 1;
			    vec_free(sample);
			    vec_free(nn);
			    break;
			}
		    }
		} while (reselect);

		if (errors[i] < min_dist) {
		    min_dist = errors[i];
		}
		
		r_pts[i] = v3_new(Vx(nn), Vy(nn), 0.0);
		l_pts[i] = v3_new(Vx(sample), Vy(sample), 0.0);

		samples[i] = sample;
		vec_free(nn);
	    }
	    end = clock();
	    sample_time += end - start;

	    /* Free all the samples */
	    for (i = 0; i < NUM_SAMPLE_POINTS; i++) 
		vec_free(samples[i]);

#ifdef RIGID_TRANSFORM
	    /* Compute the optimal translation and rotation using the Horn
	     * solver */ 
	    align_horn(NUM_SAMPLE_POINTS, r_pts, l_pts, (double *)rot, (double *)trans, (double *)Tnew->T, weights);
#endif

#ifdef AFFINE_TRANSFORM
	    align_affine(NUM_SAMPLE_POINTS, r_pts, l_pts, (double *)Tnew->T);
#endif

#ifdef HOMOGRAPHY_TRANSFORM
	    align_homography(NUM_SAMPLE_POINTS, r_pts, l_pts, (double *)Tnew->T, 0);
#endif


	    if (isnan(Tnew->T[0][0])) {
		printf("result is nan\n");
		trial--;
		continue;
	    }

	    Tprod = transform_product(Tnew, T);
	    Tinv_test = transform_invert(Tprod);

	    img_compute_Tstats(b, Tprod, &w_new, &h_new, &origin_new);


	    /* Select points inside the transformed image */
	    
	    if (img_disconnected_under_homography(b, (double *)Tnew->T)) {
		disconnected++;
	    }

	    /* Check if the image is too small or too big.  If so
	     * something strange probably happened */
	    if (!img_disconnected_under_homography(b, (double *)Tnew->T) && 
		w_new > 0.2 * b->w && w_new < 8 * b->w && h_new > 0.2 * b->h && h_new < 8 * b->h) {

		bad = 0;
		
		// #define RESAMPLE_TEST_IMAGE
#ifdef RESAMPLE_TEST_IMAGE
		img_test = img_resample_bbox(b, Tprod);
#else
		sample = vec_new(5);
#endif

		error = 0.0;

		start = clock();

		// num_test_points = img_test->w * img_test->h * RANSAC_SAMPLE_RATE;

		for (j = 0; j < num_test_points; j++) {
		    double dist;
		    fcolor_t col;

#ifndef RESAMPLE_TEST_IMAGE
		    if (sample_random_pt(w_new, h_new, origin_new, 
					 Tinv_test, 
					 w, h, b->origin,
					 nhood_radius_small, &x, &y) == 1) 
			{
			    error = DBL_MAX;
			    break;
			}

#if 1
		    /* If (x, y) lies outside the other image, try
		     * again (we only compare pixels in the
		     * intersection of both images */
		    if (!img_pixel_is_valid(a, x, y)) {
			j--;
			continue;
		    }
#endif

		    /* Compute the pixel at (x, y) in the transformed
		     * image */

		    col = pixel_transform(b, Tinv_test, x, y);
		    
		    Vx(sample) = x;
		    Vy(sample) = y;
		    Vn(sample, 2) = col.r;
		    Vn(sample, 3) = col.g;
		    Vn(sample, 4) = col.b;
#else	    
		    img_sample_random_pt(img_test, Tinv_test, b, nhood_radius_small, &x, &y, 0.0);

		    if (!img_pixel_is_valid(a, Vx(img_test->origin) + x, Vy(img_test->origin) + y)) {
			j--;
			continue;
		    }
		    
		    sample = em_create_query_pt(img_test, EM_RGB, x, y);
#endif

		    dist = query_ann_tree_dist(tree_small, sample, 0.0 /* 2.0 */, NULL);

		    error += dist;

#ifdef RESAMPLE_TEST_IMAGE
		    vec_free(sample);
#endif
		}
		end = clock();
		error_time += end - start;

		error /= num_test_points;
		
#ifdef RESAMPLE_TEST_IMAGE
		img_free(img_test);
#else
		vec_free(sample);
#endif
	    } else {
		bad = 1;
		error = DBL_MAX;
		trial--;
	    }

#if 1
	    if (verbose_ransac && !bad) {
		for (j = 0; j < NUM_SAMPLE_POINTS; j++)
		    printf("(%0.3f, %0.3f) ==> (%0.3f, %0.3f) [%0.3f, %0.3f] [%0.3f, %0.3f] %0.3f\n",
			   Vx(l_pts[j]), Vy(l_pts[j]), Vx(r_pts[j]), Vy(r_pts[j]),
			   Vx(r_pts[j]) - Vx(l_pts[j]), Vy(r_pts[j]) - Vy(l_pts[j]),
			   Vx(r_pts[j]) / Vx(l_pts[j]), Vy(r_pts[j]) / Vy(l_pts[j]),
			   errors[j]);

		printf("[Trial %04d] %0.3f %0.3f\n", trial, error, min_sample_dist);
	    }
#endif

	    // vec_free(sample);

	    /* Check if this is the new best transform */
	    if (error < min_error) {
		if (Tbest != NULL) 
		    transform_free(Tbest);
		Tbest = transform_copy(Tprod);
		min_error = error;

		memcpy(r_pts_best, r_pts, sizeof(v3_t) * NUM_SAMPLE_POINTS);
		memcpy(l_pts_best, l_pts, sizeof(v3_t) * NUM_SAMPLE_POINTS);
	    }

#ifdef RANSAC_OUTPUT_INTERMEDIATE
	    if (!img_disconnected_under_homography(b, (double *)Tnew->T) &&
		w_new > 0.2 * b->w && w_new < 8 * b->w && h_new > 0.2 * b->h && h_new < 8 * b->h) {
		Timg_test = img_resample_bbox(b, Tprod);

		// img_light_invalid_pixels(Timg_test);
		img_fused = img_fuse(2, 1, a, Timg_test);

		for (i = 0; i < NUM_SAMPLE_POINTS; i++)
		    img_draw_pt(img_fused, Vx(r_pts[i]), Vy(r_pts[i]), 4, 0x0, 0xff, 0);

		sprintf(fname, "ransac_%03d-%03d.bmp", num_ransac_calls, trial);
		
		f = fopen(fname, "w");
		if (f == NULL) {
		    printf("Cannot open %s for writing\n", fname);
		    return NULL;
		}
		
		bmp = img2bmp(img_fused);
		write_bmp(f, bmp);
		fclose(f);

		free_bmp(bmp);
		img_free(Timg_test);
		img_free(img_fused);
	    }
#endif
	    transform_free(Tprod);
	    transform_free(Tinv_test);

	    if (min_error < 1.0) {
		printf("exiting after %d trials\n", trial);
		break;
	    }
	}

	if (verbose_ransac) {
	    printf("Minimum error: %0.3f\n", min_error);
	    printf("Minimum dist: %0.3f\n", min_dist);
	    printf("Number of disconnected images: %d\n", disconnected);
	}
	
	/* Transform the image by Tbest for the next round */
	itmp = img_resample(b, Tbest);
	img_free(Timg);
	Timg = itmp;

	transform_free(T);
	T = Tbest;

	transform_free(Tinv);
	Tinv = transform_invert(T);

	if (error_out != NULL)
	    *error_out = min_error;
    }
    
    free(r_pts);
    free(l_pts);
    free(weights);
    free(samples);

    transform_free(Tinv);
    transform_free(Tnew);

    free_ann_tree(tree_small);
    // free_ann_tree(tree_large);

    print_transform(T);

    Timg_test = img_resample_bbox(b, T);
    a_cpy = img_copy(a);
    b_cpy = img_resample_bbox(b, Tin);

    if (verbose_ransac)
	printf("Matches:\n");

    Tinv = transform_invert(Tin);

    for (i = 0; i < NUM_SAMPLE_POINTS; i++) {
	int r = ((i % 2) == 0 || (i % 3) == 0) ? 255 : 0;
	int g = ((i % 2) == 1 || (i % 3) == 1) ? 255 : 0;
	int b = (i % 3) == 2 ? 255 : 0;
	double xout, yout;

	if (verbose_ransac) 
	    printf("        (%0.3f, %0.3f) ==> (%0.3f, %0.3f) [%0.3f, %0.3f]\n", 
		   Vx(l_pts_best[i]), Vy(l_pts_best[i]),
		   Vx(r_pts_best[i]), Vy(r_pts_best[i]),
		   Vx(r_pts_best[i]) - Vx(l_pts_best[i]), Vy(r_pts_best[i]) - Vy(l_pts_best[i]));

	img_draw_pt(a_cpy, Vx(r_pts_best[i]), Vy(r_pts_best[i]), 5, r, g, b);
	img_draw_pt(b_cpy, Vx(l_pts_best[i]), Vy(l_pts_best[i]), 5, r, g, b);

	transform_point(Tinv, Vx(l_pts_best[i]), Vy(l_pts_best[i]), &xout, &yout);
	transform_point(T, xout, yout, &xout, &yout);
	img_draw_pt(Timg_test, (int) rint(xout), (int) rint(yout), 5, r, g, b);
    }

    fflush(stdout);

    transform_free(Tinv);

    sprintf(fname, "ransac_post-%03d.bmp", num_ransac_calls);
    f = fopen(fname, "w");
    if (f == NULL) {
	printf("Cannot open %s for writing\n", fname);
	return NULL;
    }

    // img_light_invalid_pixels(Timg_test);
    bmp = img2bmp(Timg_test);
    write_bmp(f, bmp);
    fclose(f);

    free_bmp(bmp);


    sprintf(fname, "ransac_fuse-%03d.bmp", num_ransac_calls);
    f = fopen(fname, "w");
    if (f == NULL) {
	printf("Cannot open %s for writing\n", fname);
	return NULL;
    }

    img_fused = img_fuse(2, 0, a, Timg_test);
    bmp = img2bmp(img_fused);
    write_bmp(f, bmp);
    fclose(f);
    
    free_bmp(bmp);
    img_free(img_fused);
    img_free(Timg_test);


    sprintf(fname, "ransac_pre0-%03d.bmp", num_ransac_calls);
    f = fopen(fname, "w");
    if (f == NULL) {
	printf("Cannot open %s for writing\n", fname);
	return NULL;
    }

    bmp = img2bmp(a_cpy);
    write_bmp(f, bmp);
    fclose(f);

    free_bmp(bmp);
    img_free(a_cpy);

    sprintf(fname, "ransac_pre1-%03d.bmp", num_ransac_calls);
    f = fopen(fname, "w");
    if (f == NULL) {
	printf("Cannot open %s for writing\n", fname);
	return NULL;
    }

    bmp = img2bmp(b_cpy);
    write_bmp(f, bmp);
    fclose(f);

    free_bmp(bmp);
    img_free(b_cpy);

    free(r_pts_best);
    free(l_pts_best);

    img_free(Timg);

    num_ransac_calls++;

    start = clock();    

    use_icp = 0;
    if (use_icp) {
	icp_recurse_off();
	set_ann_z_weight(0.1); /* For nhood */
	// set_ann_z_weight(0.01); /* For nhood */
	// set_ann_z_weight(20.0); /* For ring */
	set_icp_sample_intersection(1);
	// set_icp_reject_percent(0.9);
	set_icp_pairs_sample_percent(0.05);
	set_icp_ann_epsilon(0.0);
	em_set_nhood_size(1);
	icp_set_error_metric(EM_NHOOD_RGB);
	set_icp_max_rounds(64);

	/* Run a few rounds of ICP to correct small errors */
	printf("[ransac] Running ICP...\n");
	Tnew = transform_copy(T); 
	Tnew = align_Timage_ICP(a, b, T, TRANSFORM_HOMOGRAPHY);
	transform_free(T);

	end = clock();
	correct_time = end - start;
    } else {
	correct_time = 0.0;
	Tnew = transform_copy(T);
	transform_free(T);
    }
    
    printf("[Tick-Tock] Sample time:      %0.3fs\n"
	   "[Tick-Tock] Error time:       %0.3fs\n"
	   "[Tick-Tock] Correction time:  %0.3fs\n",
	   (double) sample_time / CLOCKS_PER_SEC,
	   (double) error_time / CLOCKS_PER_SEC,
	   (double) correct_time / CLOCKS_PER_SEC);

    print_transform(Tnew);

    return Tnew;
}
