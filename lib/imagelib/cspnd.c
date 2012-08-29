/* cspnd.c */
/* Find correspondences between two images */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>

#include "cspnd.h"
#include "defines.h"
#include "dmap.h"
#include "dmap-io.h"
#include "error.h"
#include "fmatrix.h"
#include "image.h"
#include "matrix.h"
#include "pyramid.h"
#include "pyramid-io.h"
#include "qsort.h"
#include "ransac.h"
#include "resample.h"
#include "transform.h"
#include "util.h"

static double similarity_threshold = 4.0;

/* Set the similarity threshold for finding symmetric matches */
void cspnd_set_similarity_threshold(double t) {
    if (t < 0) {
	printf("[cspnd_set_similarity_threshold] Error: similarity threshold must be >= 0\n");
    } else {
	similarity_threshold = t;
    }
}

/* Find correspondences from image B to image A */
img_dmap_t *img_find_correspondence(img_t *a, img_t *b, 
                                    int diameter, img_dmap_t *amap, img_dmap_t *bmap) {
    trans2D_t *T;
    img_t *b_warp;
    img_dist_pyr_t *apyr_new, *bpyr_new;
    img_dmap_t *bmap_new;

    int x, y;

    printf("[img_find_correspondence] Aligning images\n");

    /* Align B to A using ransac */
    T = align_image_ransac_map(a, b, diameter, amap, bmap, 1, NULL);

    printf("[img_find_correspondence] Warping\n");

    /* Warp B */
    b_warp = img_resample_bbox(b, T);

    printf("[img_find_correspondence] Computing new pyramids\n");

    set_ann_z_weight(10000.0);
    em_set_nhood_size(5);

    /* Compute new distance pyramids */
    apyr_new = img_create_hierarchical_distance_pyramid(a, b_warp, EM_NHOOD_RGB, 0.0, MAX(least_larger_power_of_two(b_warp->w), least_larger_power_of_two(b_warp->h)) / 8, 0);
    // bpyr_new = img_create_hierarchical_distance_pyramid(b_warp, a, EM_NHOOD_RGB, 0.0, 32);

    bmap_new = img_dmap_new(b->w, b->h);

    printf("[img_find_correspondence] Computing correspondence\n");

    /* Compute correspondence of original B */
    for (y = 0; y < b->h; y++) {
	for (x = 0; x < b->w; x++) {
	    double x_new, y_new, dist;
	    int x_img, y_img, x_nn, y_nn;

	    transform_point(T, x + Vx(b->origin), y + Vy(b->origin), &x_new, &y_new);

	    x_img = rint(x_new - Vx(b_warp->origin));
	    y_img = rint(y_new - Vy(b_warp->origin));

	    if (x_img < 0 || x_img >= b_warp->w || y_img < 0 || y_img >= b_warp->h)
		continue;

	    if (!img_pixel_is_valid(b_warp, x_img, y_img))
		continue;

	    x_nn = Vx(apyr_new->dmaps[0].nns[y_img * apyr_new->w + x_img]);
	    y_nn = Vy(apyr_new->dmaps[0].nns[y_img * apyr_new->w + x_img]);
	    
	    dist = apyr_new->dmaps[0].dists[y_img * apyr_new->w + x_img];

	    bmap_new->dists[y * b->w + x] = dist;
	    bmap_new->nns[y * b->w + x] = v2_new(x_nn, y_nn);
	}
    }

    printf("[img_find_correspondence] Cleaning up\n");

    img_free_distance_pyramid(apyr_new);
    // img_free_distance_pyramid(bpyr_new);
    img_free(b_warp);
    transform_free(T);

    return bmap_new;
}


/* Estimate correspondences from image B to image A by detecting
 * symmetric matches.  The correspondences are returned in the form of
 * a distance map */
img_dmap_t *img_estimate_correspondence(img_t *a, img_t *b, 
                                        img_dmap_t *amap, img_dmap_t *bmap) 
{
    img_dmap_t *bmap_new = img_dmap_new(bmap->w, bmap->h);
    int x, y;

    /* Compute correspondence of original B */
    for (y = 0; y < bmap->h; y++) {
	for (x = 0; x < bmap->w; x++) {
	    /* Extract the nearest neighbor */
	    v2_t nn = bmap->nns[y * bmap->w + x];
	    double dist = bmap->dists[y * bmap->w + x];
	    double x_nn = Vx(nn);
	    double y_nn = Vy(nn);
	    double x_nn_back, y_nn_back, back_disp;
	    color_t col;
	    int black;

	    if (dist == DBL_MAX) {
		bmap_new->nns[y * bmap->w + x] = v2_new(0.0, 0.0);
		bmap_new->dists[y * bmap->w + x] = DBL_MAX;
	    } else {
		/* Hack -- check if the pixel is black */
		col = img_get_pixel(b, x, y);
		black = (col.r == 0 && col.g == 0 && col.b == 0);

		x_nn_back = Vx(amap->nns[(int)(rint(y_nn) * amap->w + rint(x_nn))]);
		y_nn_back = Vy(amap->nns[(int)(rint(y_nn) * amap->w + rint(x_nn))]);

		back_disp = sqrt((x - x_nn_back) * (x - x_nn_back) + (y - y_nn_back) * (y - y_nn_back));

		/* Check the reverse nearest neighbor */
		if (back_disp < similarity_threshold && !black) {
		    bmap_new->nns[y * bmap->w + x] = nn;
		    bmap_new->dists[y * bmap->w + x] = dist;
		} else {
		    bmap_new->nns[y * bmap->w + x] = v2_new(0.0, 0.0);
		    bmap_new->dists[y * bmap->w + x] = DBL_MAX;
		}
	    }
	}
    }

    return bmap_new;
}

#if 0
static int double_compare(const void *a, const void *b) 
{
    double ad = *((double *)a);
    double bd = *((double *)b);
    
    if (ad < bd)
	return -1;
    else if (ad == bd)
	return 0;
    else
	return 1;
}
#endif

static int num_fmatrix_trials = 2048;

/* Set the number of F-matrix estimation RANSAC trials */
void set_num_fmatrix_trials(int n) {
    if (n >= 0)
	num_fmatrix_trials = n;
}

#if 0
/* Estimate the fundamental matrix for the images `a' and `b' given an
 * initial set of correspondences `bmap', then filter the set to only
 * contain matches that satisfy the epipolar constraint */
void img_estimate_correspondence_epipolar(img_t *a, img_t *b, img_dmap_t *amap, img_dmap_t *bmap, 
					  double *F, double *e1, double *e2, img_dmap_t **b2a, img_dmap_t **a2b) {
    int i, j, k, x, y, count, idx, num_cspnd;

    v3_t *a_cspnd, *b_cspnd;
    v3_t l_pts_best[8], r_pts_best[8];

    double Fbest[9];
    double *resid;
    double error_min;

    double threshold = 1.0e-10;

    img_dmap_t *bmap_new = img_dmap_new(bmap->w, bmap->h);

    img_t *a_cpy, *b_cpy;
    bmp_t *a_bmp, *b_bmp;

    FILE *f;

    srand(time(0));

    /* Make an array of all good correspondences */
    num_cspnd = 0;
    for (i = 0; i < bmap->w * bmap->h; i++) {
	if (bmap->dists[i] != DBL_MAX) 
	    num_cspnd++;
    }

    if (num_cspnd < 8) {
	printf("[img_estimate_correspondence_epipolar] Could not find 8 good correspondences, F-matrix estimation failed\n");

	img_dmap_free(bmap_new);

	return;
    }

    a_cspnd = (v3_t *)malloc(sizeof(v3_t) * num_cspnd);
    b_cspnd = (v3_t *)malloc(sizeof(v3_t) * num_cspnd);

    count = idx = 0;
    for (y = 0; y < bmap->h; y++) {
	for (x = 0; x < bmap->w; x++) {
	    if (bmap->dists[idx] != DBL_MAX) {
		v2_t nn = bmap->nns[idx];

		a_cspnd[count] = v3_new(Vx(nn), Vy(nn), 1.0);
		b_cspnd[count] = v3_new(x, y, 1.0);

		count++;
	    }

	    idx++;
	}
    }

    error_min = DBL_MAX;
    resid = (double *)malloc(sizeof(double) * num_cspnd);

    /* Estimate the F-matrix using RANSAC */
    for (i = 0; i < num_fmatrix_trials; i++) {
	int idxs[8];
	v3_t l_pts[8], r_pts[8];
	double Ftmp[9], e1_tmp[3], e2_tmp[3];
	double error;

	/* Sample 8 random correspondences */
	for (j = 0; j < 8; j++) {
	    int reselect = 0;

	    idx = rand() % num_cspnd;
	    
	    /* Make sure we didn't sample this index yet */
	    for (k = 0; k < j - 1; k++) {
		if (idx == idxs[k]) {
		    reselect = 1;
		    break;
		}
	    }

	    if (reselect) {
		i--;
		continue;
	    }

	    idxs[j] = idx;
	}

	/* Fill in the left and right points */
	for (j = 0; j < 8; j++) {
	    l_pts[j] = b_cspnd[idxs[j]];
	    r_pts[j] = a_cspnd[idxs[j]];
	}

	/* Estimate the F-matrix */
	estimate_fmatrix_linear(8, r_pts, l_pts, Ftmp, e1_tmp, e2_tmp);

	/* Compute residuals */
	for (j = 0; j < num_cspnd; j++)
	    resid[j] = fmatrix_compute_residual(Ftmp, a_cspnd[j], b_cspnd[j]);

	/* Find the median */
	error = kth_element_copy(num_cspnd, num_cspnd / 2, resid); // median(num_cspnd, resid);

	if (error < error_min) {
	    error_min = error;
	    memcpy(Fbest, Ftmp, sizeof(double) * 9);
	    memcpy(e1, e1_tmp, sizeof(double) * 3);
	    memcpy(e2, e2_tmp, sizeof(double) * 3);
	    memcpy(l_pts_best, l_pts, sizeof(v3_t) * 8);
	    memcpy(r_pts_best, r_pts, sizeof(v3_t) * 8);
	}

	if (error < threshold)
	    break;
    }

    printf("Number of trials: %d\n", i);
    printf("Minimum error: %0.5e\n", error_min);

    matrix_print(3, 3, Fbest);

    /* Now use non-linear least-squares to refine the F-matrix */
    printf("Refining...\n");
#if 1
    memcpy(F, Fbest, 9 * sizeof(double));
#else
    printf("num_cspnd = %d\n", num_cspnd);
    refine_fmatrix_nonlinear(bmap, Fbest, F);
#endif

    printf("Epipoles:\n");

    e1[0] /= e1[2];
    e1[1] /= e1[2];
    e1[2] /= e1[2];

    e2[0] /= e2[2];
    e2[1] /= e2[2];
    e2[2] /= e2[2];

    matrix_print(1, 3, e1);
    matrix_print(1, 3, e2);

    free(a_cspnd);
    free(b_cspnd);
    free(resid);

    /* Output temporary images */
    a_cpy = img_copy(a);
    b_cpy = img_copy(b);

    for (i = 0; i < 8; i++) {
	img_draw_pt(a_cpy, Vx(r_pts_best[i]), Vy(r_pts_best[i]), 2, 255, 0, 0);
	img_draw_pt(b_cpy, Vx(l_pts_best[i]), Vy(l_pts_best[i]), 2, 255, 0, 0);
    }

    a_bmp = img2bmp(a_cpy);
    b_bmp = img2bmp(b_cpy);

    f = fopen("a_cspnd.bmp", "w");
    write_bmp(f, a_bmp);
    fclose(f);

    f = fopen("b_cspnd.bmp", "w");
    write_bmp(f, b_bmp);
    fclose(f);

    img_free(a_cpy);
    img_free(b_cpy);
    free_bmp(a_bmp);
    free_bmp(b_bmp);

    /* Copy out the F-matrix */
    // memcpy(F, Fbest, sizeof(double) * 9);

    /* Filter the correspondences using the epipolar constraint */
    *b2a = img_dmap_new(bmap->w, bmap->h);
    *a2b = img_dmap_new(amap->w, amap->h);

    /* Initialize all of a2b's distances to DBL_MAX */
    for (idx = 0; idx < amap->w * amap->h; idx++)
	(*a2b)->dists[idx] = DBL_MAX;

    idx = 0;
    for (y = 0; y < bmap->h; y++) {
	for (x = 0; x < bmap->w; x++) {
	    if (bmap->dists[idx] != DBL_MAX) {
		v2_t a_pt2d = bmap->nns[idx];

		v3_t a_pt = v3_new(Vx(a_pt2d), Vy(a_pt2d), 1.0);
		v3_t b_pt = v3_new(x, y, 1.0);

		double error = fmatrix_compute_residual(F, a_pt, b_pt);

		if (error > error_min) {
		    (*b2a)->dists[idx] = DBL_MAX;
		} else {
		    int ax, ay, aidx;

		    ax = Vx(a_pt2d);
		    ay = Vy(a_pt2d);

		    aidx = ay * amap->w + ax;

		    (*b2a)->dists[idx] = bmap->dists[idx];
		    (*b2a)->nns[idx] = bmap->nns[idx];
		    (*a2b)->dists[aidx] = bmap->dists[idx];
		    (*a2b)->nns[aidx] = v2_new(x, y);
		}
	    } else {
		(*b2a)->dists[idx] = DBL_MAX;
	    }

	    idx++;
	}
    }

    // return bmap_new;
}
#endif
