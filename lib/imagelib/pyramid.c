/* pyramid.c */
/* Routines for creating a gaussian pyramid from an image */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "anniface.h"
#include "dmap-io.h"
#include "error.h"
#include "fileio.h"
#include "image.h"
#include "pyramid.h"
#include "pyramid-io.h"
#include "resample.h"
#include "util.h"

/* Initialize a new image pyramid with the given width and height */
img_pyr_t *img_pyramid_new(int w, int h) {
    int n, num_levels;
    img_pyr_t *pyr;

    if (!is_power_of_two(w) || w != h) {
	printf("[img_pyramid_new] Error: invalid width and height\n");
	return NULL;
    }
    
    n = w;
    num_levels = ilog2(n);

    /* Initialize the pyramid */
    pyr = (img_pyr_t *)malloc(sizeof(img_pyr_t));
    pyr->w = n;  pyr->h = n;
    pyr->num_levels = num_levels;
    pyr->imgs = (img_t **)malloc(sizeof(img_t *) * num_levels);

    return pyr;
}

/* Creates a gaussian pyramid from an image */
img_pyr_t *img_create_gaussian_pyramid(img_t *img, int extra) {
    /* Assumes that the image size is NxN where N=2^k */
    int num_levels;
    img_pyr_t *pyr;
    int i, n, x, y;
    color_t *pxls[4];

    img_t *i_new = img_upsize_square_power_of_two(img);

    printf("Creating pyramid for (%d,%d)\n", img->w, img->h);

    pyr = img_pyramid_new(i_new->w, i_new->h);
    num_levels = pyr->num_levels;
    printf("num_levels = %d\n", num_levels);
    n = i_new->w;

    /* Copy over the first image */
    pyr->imgs[extra] = img_copy(i_new);

    /* Resample the bigger images */
    for (i = 0; i < extra; i++) {
	trans2D_t *T = new_scaling_transform((1 << (extra - i)), (1 << (extra - i)));
	pyr->imgs[i] = img_resample_bbox(img, T);
	pyr->imgs[i]->origin = v2_scale((1 << (extra - i)), img->origin);
	transform_free(T);
    }

    /* Resample the smaller images */
    for (i = extra + 1; i < num_levels; i++) {
	pyr->imgs[i] = img_new((n >> (i - extra)), (n >> (i - extra)));
	pyr->imgs[i]->origin = v2_scale(1.0 / (1 << (i - extra)), img->origin);

	for (y = 0; y < (n >> (i - extra)); y++) {
	    for (x = 0; x < (n >> (i - extra)); x++) {
		color_t p;

		int xmin = (x << 1);
		int ymin = (y << 1);

		if (img_region_is_valid(pyr->imgs[i-1], xmin, xmin + 1, ymin, ymin + 1)) {
		    img_get_pixel_square(pyr->imgs[i-1], x << 1, y << 1, pxls);
		    
		    p.r = (pxls[0]->r + pxls[1]->r + pxls[2]->r + pxls[3]->r) >> 2;
		    p.g = (pxls[0]->g + pxls[1]->g + pxls[2]->g + pxls[3]->g) >> 2;
		    p.b = (pxls[0]->b + pxls[1]->b + pxls[2]->b + pxls[3]->b) >> 2;
		    
		    img_set_pixel(pyr->imgs[i], x, y, p.r, p.g, p.b);
		}
	    }
	}
    }

#if 1
    for (i = 0; i < num_levels; i++) {
	img_t *img_tmp = pyr->imgs[i];
	pyr->imgs[i] = img_shrink_wrap(pyr->imgs[i]);
	img_free(img_tmp);
    }
#endif

    return pyr;
}

/* Create a pyramid of kd trees for the given image */
ANNkd_pyramid_t *img_create_ann_pyramid(img_pyr_t *pyr, error_metric_t em, int *nhood_radius) {
    int i, num_levels;
    ANNkd_pyramid_t *kd_pyr;

    num_levels = pyr->num_levels;

    kd_pyr = (ANNkd_pyramid_t *)malloc(sizeof(ANNkd_pyramid_t));
    kd_pyr->num_levels = num_levels;
    kd_pyr->trees = (ANNkd_tree_t **)malloc(sizeof(ANNkd_tree_t *) * num_levels);

    for (i = 0; i < num_levels; i++)
	kd_pyr->trees[i] = img_create_ann_tree(pyr->imgs[i], em, nhood_radius);
    
    return kd_pyr;
}

/* Free an ANN kd pyramid */
void ann_pyramid_free(ANNkd_pyramid_t *pyr) {
    int i;
    for (i = 0; i < pyr->num_levels; i++) {
	free_ann_tree(pyr->trees[i]);
    }

    free(pyr->trees);
    free(pyr);
}

/* Create a pyramid of distance maps for the given images */
img_dist_pyr_t *img_create_distance_pyramid(img_t *a, img_t *b, error_metric_t em, double eps) {
    int rad;
    img_pyr_t *apyr = img_create_gaussian_pyramid(a, 0);
    img_pyr_t *bpyr = img_create_gaussian_pyramid(b, 0);
    ANNkd_pyramid_t *akdpyr = img_create_ann_pyramid(apyr, em, &rad);
    img_dist_pyr_t *dist_pyr;
    int x, y, i, num_levels;

    num_levels = akdpyr->num_levels;

    /* Initialize the distance pyramid */
    dist_pyr = (img_dist_pyr_t *)malloc(sizeof(img_dist_pyr_t));

    // dist_pyr->dists = (double **)malloc(sizeof(double *) * num_levels);
    // dist_pyr->nns = (v2_t **)malloc(sizeof(v2_t *) * num_levels);

    dist_pyr->dmaps = (img_dmap_t *)malloc(sizeof(img_dmap_t) * num_levels);

    dist_pyr->num_levels = num_levels;
    dist_pyr->w = b->w;
    dist_pyr->h = b->h;

    dist_pyr->em = em;
    dist_pyr->zweight = get_ann_z_weight();
    dist_pyr->nhood_radius = rad;

    /* Fill in each level of the distance pyramid */
    for (i = 0; i < num_levels; i++) {
	dist_pyr->dmaps[i].w = bpyr->imgs[i]->w;
	dist_pyr->dmaps[i].h = bpyr->imgs[i]->h;

	dist_pyr->dmaps[i].dists = malloc(sizeof(double) * bpyr->imgs[i]->w * bpyr->imgs[i]->h);
	dist_pyr->dmaps[i].nns = malloc(sizeof(v2_t) * bpyr->imgs[i]->w * bpyr->imgs[i]->h);
	dist_pyr->dmaps[i].uppers = NULL;

	for (y = 0; y < bpyr->imgs[i]->h; y++) {
	    for (x = 0; x < bpyr->imgs[i]->w; x++) {
		double dist;
		vec_t sample = em_create_query_pt(bpyr->imgs[i], em, x, y);
		vec_t nn = query_ann_tree_pt(akdpyr->trees[i], sample, eps, &dist, 0, 1);

		dist_pyr->dmaps[i].nns[y * bpyr->imgs[i]->w + x] = v2_new(Vx(nn), Vy(nn));
		dist_pyr->dmaps[i].dists[y * bpyr->imgs[i]->w + x] = dist;
		
		vec_free(sample);
		vec_free(nn);
	    }
	}
    }

    img_pyr_free(apyr);
    img_pyr_free(bpyr);
    ann_pyramid_free(akdpyr);

    return dist_pyr;
}

#define UPDATE_ALL_NHOOD_PIXELS
// #define CHAMFER

static int verbose_pyramid = 0;

/* Turn on verbose printing for pyramid creation */
void set_verbose_pyramid(int on) {
    if (on)
	verbose_pyramid = 1;
    else
	verbose_pyramid = 0;
}

/* Create a pyramid of distance maps for the given images */
img_dist_pyr_t *img_create_hierarchical_distance_pyramid(img_t *a, img_t *b, 
							 error_metric_t em, double eps, int min_size,
							 int rectified) 
{
    img_pyr_t *apyr;
    img_pyr_t *bpyr;
    ANNkd_tree_t *atree, *atree_full;
    img_dist_pyr_t *dist_pyr;
    int x, y, i, j, num_levels, max_level = -1, n, rad = 0, num_samples_freed;

    apyr = img_create_gaussian_pyramid(a, 0);
    bpyr = img_create_gaussian_pyramid(b, 0);

    num_levels = bpyr->num_levels;
    n = bpyr->w;

    printf("num_levels is %d\n", num_levels);

    /* Compute the maximum level of the pyramid */
    for (i = 0; i < num_levels; i++) {
	if ((n >> i) < min_size) {
	    max_level = i - 1;
	    break;
	}
    }

    if (max_level < 0) {
	printf("[img_pyramid_hierarchical_remap] Error: min_size is too large\n");
	img_pyr_free(apyr);
	img_pyr_free(bpyr);

	return NULL;
    }

    num_levels = max_level + 1;

    /* Initialize the distance pyramid */
    dist_pyr = (img_dist_pyr_t *)malloc(sizeof(img_dist_pyr_t));

    dist_pyr->dmaps = (img_dmap_t *)malloc(sizeof(img_dmap_t) * num_levels);

    dist_pyr->num_levels = num_levels;
    dist_pyr->w = bpyr->imgs[0]->w;
    dist_pyr->h = bpyr->imgs[0]->h;

    dist_pyr->em = (u_int16_t) em;
    dist_pyr->zweight = get_ann_z_weight();
    dist_pyr->nhood_radius = rad;

    /* Fill in each level of the distance pyramid from the highest
     * level */
    for (i = max_level; i >= 0; i--) {
	int w = bpyr->imgs[i]->w;
	int h = bpyr->imgs[i]->h;

	int w_other = apyr->imgs[i]->w;
	int h_other = apyr->imgs[i]->h;

	dist_pyr->dmaps[i].w = bpyr->imgs[i]->w;
	dist_pyr->dmaps[i].h = bpyr->imgs[i]->h;

	dist_pyr->dmaps[i].dists = malloc(sizeof(double) * bpyr->imgs[i]->w * bpyr->imgs[i]->h);
	dist_pyr->dmaps[i].nns = calloc(bpyr->imgs[i]->w * bpyr->imgs[i]->h, sizeof(v2_t));
	dist_pyr->dmaps[i].uppers = NULL;

	if (verbose_pyramid)
	    printf("[pyramid] Creating level %d, size = [%d x %d]\n", i, w, h);

	if (i == max_level) {

	    /* Initialize pixel distances */
	    for (y = 0; y < h; y++)
		for (x = 0; x < w; x++)
		    dist_pyr->dmaps[i].dists[y * w + x] = DBL_MAX;

	    /* Create the entire tree for the smallest level */
	    if (rectified)
		atree = img_create_ann_tree_rectified(apyr->imgs[i], em, &rad);
	    else
		atree = img_create_ann_tree(apyr->imgs[i], em, &rad);
	    
	    for (y = rad; y < bpyr->imgs[i]->h - rad; y++) {
		if (verbose_pyramid) {
		    printf(".");
		    fflush(stdout);
		}

		for (x = rad; x < bpyr->imgs[i]->w - rad; x++) {
		    if (img_region_is_valid(bpyr->imgs[i], x - rad, x + rad, y - rad, y + rad)) {
			double dist;
			vec_t sample = em_create_query_pt(bpyr->imgs[i], em, x, y);
#ifdef CHAMFER
			vec_t nn = query_ann_tree_pt_chamfer(atree, sample, &dist);
#else
			vec_t nn = query_ann_tree_pt(atree, sample, eps, &dist, 0, 1);
#endif

			dist_pyr->dmaps[i].nns[y * bpyr->imgs[i]->w + x] = v2_new(Vx(nn), Vy(nn));
			dist_pyr->dmaps[i].dists[y * bpyr->imgs[i]->w + x] = dist;
			
			vec_free(sample);
			vec_free(nn);
		    }
		}

	    }
	    
	    free_ann_tree(atree);

	    if (verbose_pyramid) {
		printf("\n");
		fflush(stdout);
	    }
	} else {
	    int w_prev = bpyr->imgs[i+1]->w;
	    int h_prev = bpyr->imgs[i+1]->h;

	    int *idx_map = (int *)malloc(sizeof(int) * w * h);
	    vec_t *samples = (vec_t *)calloc(w * h, sizeof(vec_t)); /* Cache for samples */

	    /* Initialize pixel distances */
	    for (y = 0; y < h; y++)
		for (x = 0; x < w; x++)
		    dist_pyr->dmaps[i].dists[y * w + x] = DBL_MAX;

	    /* Now build the full tree for this level */
	    if (rectified)
		atree_full = img_create_ann_tree_idx_map_rectified(apyr->imgs[i], em, &rad, idx_map);
	    else
		atree_full = img_create_ann_tree_idx_map(apyr->imgs[i], em, &rad, idx_map);

	    for (y = rad; y < h_prev - rad; y++) {
		if (verbose_pyramid) {
		    printf(".");
		    fflush(stdout);
		}

		for (x = rad; x < w_prev - rad; x++) {

		    /* Check if the region of the bigger image is valid */
		    if (img_region_is_valid(bpyr->imgs[i+1], x - rad, x + rad, y - rad, y + rad)) {
		    
			/* Find the matching neighborhood one level higher */
			int mx = (int) rint(Vx(dist_pyr->dmaps[i+1].nns[y * w_prev + x]) - Vx(apyr->imgs[i+1]->origin));
			int my = (int) rint(Vy(dist_pyr->dmaps[i+1].nns[y * w_prev + x]) - Vy(apyr->imgs[i+1]->origin));
			int xp, yp;
			
			int xmin = 2 * (mx - rad), xmax = 2 * (mx + rad) + 1;
			int ymin = 2 * (my - rad), ymax = 2 * (my + rad) + 1;

			/* Create a search tree */
#ifdef REBUILD_ANN_TREE
			atree = img_create_ann_subtree(apyr->imgs[i], em, &rad, xmin, xmax, ymin, ymax);
#else
			atree = img_create_ann_subtree(apyr->imgs[i], em, &rad, xmin, xmax, ymin, ymax, rectified, 1, atree_full, idx_map);
#endif

#ifdef UPDATE_ALL_NHOOD_PIXELS
			/* Update the information for a window of
			 * pixels that are contained in the larger
			 * neighborhood */
			for (yp = -rad; yp <= rad /* + 1 */; yp++) {
			    int ynew = 2 * y + yp;

			    for (xp = -rad; xp <= rad /* + 1 */; xp++) {

				int xnew = 2 * x + xp;

				double dist;
				vec_t sample = (samples[ynew * w + xnew].d == 0) ? 
				    em_create_query_pt(bpyr->imgs[i], em, xnew, ynew) : samples[ynew * w + xnew];

#ifdef CHAMFER
				vec_t nn = query_ann_tree_pt_chamfer(atree, sample, &dist);
#else
				vec_t nn = query_ann_tree_pt(atree, sample, eps, &dist, 0, 1);
				// vec_t nn = query_ann_tree_pt(atree, sample, &dist);
#endif

				double old_dist = dist_pyr->dmaps[i].dists[ynew * w + xnew];

				if (dist < old_dist) {
				    dist_pyr->dmaps[i].nns[ynew * w + xnew] = v2_new(Vx(nn), Vy(nn));
				    dist_pyr->dmaps[i].dists[ynew * w + xnew] = dist;
				}
			    
				samples[ynew * w + xnew] = sample;

				vec_free(nn);
			    }
			}
#else
			/* Update only the four pixels that make up
			 * the center of the higher level */
			for (yp = 0; yp <= 1; yp++) {
			    int ynew = 2 * y + yp;

			    for (xp = 0; xp <= 1; xp++) {

				int xnew = 2 * x + xp;

				double dist;
				vec_t sample = (samples[ynew * w + xnew].d == 0) ? 
				    em_create_query_pt(bpyr->imgs[i], em, xnew, ynew) : samples[ynew * w + xnew];

				// vec_t sample = em_create_query_pt(bpyr->imgs[i], em, xnew, ynew);
#ifdef CHAMFER
				vec_t nn = query_ann_tree_pt_chamfer(atree, sample, &dist);
#else
				vec_t nn = query_ann_tree_pt(atree, sample, eps, &dist);
#endif	
			
				double old_dist = dist_pyr->dmaps[i].dists[ynew * w + xnew];

				if (dist < old_dist) {
				    dist_pyr->dmaps[i].nns[ynew * w + xnew] = v2_new(Vx(nn), Vy(nn));
				    dist_pyr->dmaps[i].dists[ynew * w + xnew] = dist;
				}
			    
				// vec_free(sample);
				samples[ynew * w + xnew] = sample;

				vec_free(nn);
			    }
			}

#endif
			free_ann_tree(atree);
		    }
		}
	    }

	    if (verbose_pyramid) {
		printf("\n");
		fflush(stdout);
	    }

	    /* Fill in any holes in the map */
	    for (y = rad; y < h - rad; y++) {
		for (x = rad; x < w - rad; x++) {
		    if (img_region_is_valid(bpyr->imgs[i], x - rad, x + rad, y - rad, y + rad)) {
			if (dist_pyr->dmaps[i].dists[y * w + x] == DBL_MAX) {

#if 1
			    double dist;
			    vec_t sample = em_create_query_pt(bpyr->imgs[i], em, x, y);
			    vec_t nn = query_ann_tree_pt(atree_full, sample, eps, &dist, 0, 1);
				
			    double old_dist = dist_pyr->dmaps[i].dists[y * w + x];

			    // printf("Filling in (%d, %d)\n", x, y);

			    if (dist < old_dist) {
				dist_pyr->dmaps[i].nns[y * w + x] = v2_new(Vx(nn), Vy(nn));
				dist_pyr->dmaps[i].dists[y * w + x] = dist;
			    }
			    
			    vec_free(sample);
			    vec_free(nn);			    
#else
			    dist_pyr->dmaps[i].dists[y * w + x] = 0.0;
			    dist_pyr->dmaps[i].nns[y * w + x] = v2_new(x, y);
#endif

			}
		    }
		}
	    }

	    /* Free samples */
	    num_samples_freed = 0;
	    for (j = 0; j < w * h; j++) {
		if (samples[j].d != 0) {
		    vec_free(samples[j]);
		    num_samples_freed++;
		}
	    }
	    printf("Freed %d samples\n", num_samples_freed);
	    free(samples);

	    free_ann_tree(atree_full);
	    free(idx_map);
	}

#if 1
	/* Check that we updated all valid entries */
	for (y = rad; y < h - rad; y++) {
	    for (x = rad; x < w - rad; x++) {
		v2_t nn = dist_pyr->dmaps[i].nns[y * w + x];

		nn = v2_sub(nn, apyr->imgs[i]->origin);
		    
		if (img_region_is_valid(bpyr->imgs[i], x - rad, x + rad, y - rad, y + rad)) {
		    if (Vx(nn) < rad || Vy(nn) < rad ||
			Vx(nn) > w_other - rad - 1 || Vy(nn) > h_other - rad - 1)
			printf("[img_create_hdp] Error: nn [%0.3f, %0.3f] out of range at (%d, %d)\n", 
			       Vx(nn), Vy(nn), x, y);
		    
		    if (dist_pyr->dmaps[i].dists[y * w + x] == DBL_MAX)
			printf("[img_create_hdp] Error: (%d, %d)!\n", x, y);
		}
	    }
	}
#endif

    }

    img_pyr_free(apyr);
    img_pyr_free(bpyr);

    return dist_pyr;
}

/* Create a pyramid of distance maps for the given images (take two) */
img_dist_pyr_t *img_create_hierarchical_distance_pyramid_2(img_t *a, img_t *b, 
							   error_metric_t em, double eps, int min_size,
							   int rectified) 
{
    img_pyr_t *apyr;
    img_pyr_t *bpyr;
    ANNkd_tree_t *atree, *atree_full;
    img_dist_pyr_t *dist_pyr;
    int x, y, i, num_levels, max_level = -1, n, rad = 0;

    apyr = img_create_gaussian_pyramid(a, 0);
    bpyr = img_create_gaussian_pyramid(b, 0);

    num_levels = bpyr->num_levels;
    n = bpyr->w;

    printf("num_levels is %d\n", num_levels);

    /* Compute the maximum level of the pyramid */
    for (i = 0; i < num_levels; i++) {
	printf("(n >> i) = %d, min_size = %d\n", (n >> i), min_size);
	if ((n >> i) < min_size) {
	    max_level = i - 1;
	    break;
	}
    }

    if (max_level < 0) {
	printf("[img_pyramid_hierarchical_remap] Error: min_size is too large\n");
	img_pyr_free(apyr);
	img_pyr_free(bpyr);

	return NULL;
    }

    num_levels = max_level + 1;

    /* Initialize the distance pyramid */
    dist_pyr = (img_dist_pyr_t *)malloc(sizeof(img_dist_pyr_t));

    dist_pyr->dmaps = (img_dmap_t *)malloc(sizeof(img_dmap_t) * num_levels);

    dist_pyr->num_levels = num_levels;
    dist_pyr->w = bpyr->imgs[0]->w;
    dist_pyr->h = bpyr->imgs[0]->h;

    dist_pyr->em = (u_int16_t) em;
    dist_pyr->zweight = get_ann_z_weight();
    dist_pyr->nhood_radius = rad;

    /* Fill in each level of the distance pyramid from the highest
     * level */
    for (i = max_level; i >= 0; i--) {
	int w = bpyr->imgs[i]->w;
	int h = bpyr->imgs[i]->h;

	int w_other = apyr->imgs[i]->w;
	int h_other = apyr->imgs[i]->h;

	dist_pyr->dmaps[i].w = bpyr->imgs[i]->w;
	dist_pyr->dmaps[i].h = bpyr->imgs[i]->h;

	dist_pyr->dmaps[i].dists = malloc(sizeof(double) * bpyr->imgs[i]->w * bpyr->imgs[i]->h);
	dist_pyr->dmaps[i].nns = calloc(bpyr->imgs[i]->w * bpyr->imgs[i]->h, sizeof(v2_t));
	dist_pyr->dmaps[i].uppers = NULL;

	if (verbose_pyramid)
	    printf("[pyramid] Creating level %d, size = [%d x %d]\n", i, w, h);

	if (i == max_level) {

	    /* Initialize pixel distances */
	    for (y = 0; y < h; y++)
		for (x = 0; x < w; x++)
		    dist_pyr->dmaps[i].dists[y * w + x] = DBL_MAX;

	    /* Create the entire tree for the smallest level */
	    if (rectified) {
		if (i == 0) {
		    atree = img_create_ann_tree_rectified(apyr->imgs[i], em, &rad);
		} else {
		    atree = img_create_ann_tree_params(apyr->imgs[i-1], em, &rad, 1, 2);
		}
	    } else {
		atree = img_create_ann_tree(apyr->imgs[i], em, &rad);
	    }

	    for (y = rad; y < bpyr->imgs[i]->h - rad; y++) {
		if (verbose_pyramid) {
		    printf(".");
		    fflush(stdout);
		}

		for (x = rad; x < bpyr->imgs[i]->w - rad; x++) {
		    if (img_region_is_valid(bpyr->imgs[i], x - rad, x + rad, y - rad, y + rad)) {
			double dist;
			vec_t sample = em_create_query_pt(bpyr->imgs[i], em, x, y);
#ifdef CHAMFER
			vec_t nn = query_ann_tree_pt_chamfer(atree, sample, &dist);
#else
			// vec_t nn = query_ann_tree_pt_brute_force(atree, sample, &dist);
			vec_t nn = query_ann_tree_pt(atree, sample, eps, &dist, 0, 1);
#endif

#if 0
			if (y == 4 && x == 25) {
			    int iter;
			    
			    printf("sample:\n");
			    for (iter = 0; iter < sample.d; iter++) {
				printf("%0.3f\n", sample.p[iter]);
			    }
			    
			    printf("nn:\n");
			    for (iter = 0; iter < nn.d; iter++) {
				printf("%0.3f\n", nn.p[iter]);
			    }			    
			}
#endif

			dist_pyr->dmaps[i].nns[y * bpyr->imgs[i]->w + x] = v2_new(Vx(nn), Vy(nn));
			dist_pyr->dmaps[i].dists[y * bpyr->imgs[i]->w + x] = dist;
			
			vec_free(sample);
			vec_free(nn);
		    }
		}

	    }
	    
	    free_ann_tree(atree);

	    if (verbose_pyramid) {
		printf("\n");
		fflush(stdout);
	    }
	} else {
	    int w_prev = bpyr->imgs[i+1]->w;
	    int h_prev = bpyr->imgs[i+1]->h;

	    /* Allocate upper level pointers */
	    img_dmap_allocate_uppers(&(dist_pyr->dmaps[i]));

	    // int *idx_map = (int *)malloc(sizeof(int) * w * h);
	    // vec_t *samples = (vec_t *)calloc(w * h, sizeof(vec_t)); /* Cache for samples */

	    /* Initialize pixel distances */
	    for (y = 0; y < h; y++) {
		for (x = 0; x < w; x++) {
		    dist_pyr->dmaps[i].dists[y * w + x] = DBL_MAX;
		    dist_pyr->dmaps[i].uppers[y * w + x] = iv2_new(0, 0);
		}
	    }

	    /* Now build the full tree for this level */
	    if (em == EM_NHOOD_RGB)
		set_ann_z_weight(0.1);
	    else if (em == EM_RING_HISTOGRAM_RGB) {
		set_ann_z_weight(20.0);
	    } else if (em == EM_RING_HISTOGRAM_GS) {
		set_ann_z_weight(100.0);
	    } else {
		set_ann_z_weight(0.1);
	    }

	    if (rectified) {
		if (i == 0) {
		    atree_full = img_create_ann_tree_idx_map_rectified(apyr->imgs[i], em, &rad, NULL);
		} else {
		    /* Downsize by 2 */
		    atree_full = img_create_ann_tree_params(apyr->imgs[i-1], em, &rad, 1, 2);
		}
	    } else {
		atree_full = img_create_ann_tree_idx_map(apyr->imgs[i], em, &rad, NULL);
	    }
	    
	    for (y = rad; y < h - rad; y++) {
		if (verbose_pyramid) {
		    printf(".");
		    fflush(stdout);
		}

		for (x = rad; x < w - rad; x++) {
		    int x_upper = x / 2, y_upper = y / 2;
		    int x2, y2;
		    vec_t sample, nn;
		    double min_match = DBL_MAX;
		    int x2_min = 0, y2_min = 0, x2_min2 = 0, y2_min2 = 0;
		    double dist;
		    v2_t nbr;
		    double nnx_min = 0.0, nny_min = 0.0;

		    /* Check if the region of the bigger image is valid */
		    if (!img_region_is_valid(bpyr->imgs[i], x - rad, x + rad, y - rad, y + rad)) {
			continue;
		    }
		    

		    /* Find the best nearby match on the last level */
		    for (y2 = y_upper - 1; y2 <= y_upper + 1; y2++) {
			for (x2 = x_upper - 1; x2 <= x_upper + 1; x2++) {
			    double dist;

			    if (x2 < 0 || x2 >= w_prev || y2 < 0 || y2 >= h_prev)
				continue;

			    dist = dist_pyr->dmaps[i+1].dists[y2 * w_prev + x2];

			    if (dist == DBL_MAX)
				continue;

			    if (dist < min_match) {
				min_match = dist;
				x2_min = x2;
				y2_min = y2;
			    }
			}
		    }
		    
		    if (min_match == DBL_MAX) {
			/* No good match? */
			continue;
		    }

		    // if (i == 3 && x == 22 && y == 37) {
		    //    printf("hi\n");
		    // }

		    /* Create a sample */
		    sample = em_create_query_pt(bpyr->imgs[i], em, x, y);
		    nbr = dist_pyr->dmaps[i+1].nns[y2_min * w_prev + x2_min];

		    dist_pyr->dmaps[i].uppers[y * w + x] = iv2_new(x2_min, y2_min);

		    /* Center the sample around the best previous
		     * match */

		    /* Idea... maybe vary x over a small 2x2 window
		     * and pick the lowest cost result? */
		    min_match = DBL_MAX;

#if 0		    
		    for (y2 = rectified - 1; y2 <= 1 - rectified; y2++) {
			for (x2 = -1; x2 <= 1; x2++) {
			    Vx(sample) = (2 * Vx(nbr) + 0.5) + (x - (2 * x2_min + 0.5)) + x2;
		    
			    if (!rectified)
				Vy(sample) = 2 * Vy(nbr) - (2 * y2_min - y) + y2;

			    nn = query_ann_tree_pt(atree_full, sample, eps, &dist, 0);

			    if (dist < min_match) {
				min_match = dist;
				x2_min2 = x2;
				y2_min2 = y2;
				nnx_min = Vx(nn);
				nny_min = Vy(nn);
			    }

			    vec_free(nn);
			}
		    }

		    if (min_match == DBL_MAX) {
			/* No good match found? */
			printf("[img_create_hierarchical_distance_pyramid_2] Error: no good match\n");
		    } else {
			dist_pyr->dmaps[i].dists[y * w + x] = min_match;
			dist_pyr->dmaps[i].nns[y * w + x] = v2_new(nnx_min, nny_min);
		    }

#else
		    Vx(sample) = (2 * Vx(nbr)) + (x - (2 * x2_min));
		    
		    if (!rectified)
			Vy(sample) = 2 * Vy(nbr) - (2 * y2_min - y);

		    nn = query_ann_tree_pt(atree_full, sample, eps, &dist, 0, 1);

		    dist_pyr->dmaps[i].dists[y * w + x] = dist;
		    dist_pyr->dmaps[i].nns[y * w + x] = v2_new(Vx(nn), Vy(nn));
#endif

		    vec_free(sample);
		}
	    }

	    if (verbose_pyramid) {
		printf("\n");
		fflush(stdout);
	    }

	    /* Fill in any holes in the map */
	    for (y = rad; y < h - rad; y++) {
		for (x = rad; x < w - rad; x++) {
		    if (img_region_is_valid(bpyr->imgs[i], x - rad, x + rad, y - rad, y + rad)) {
			if (dist_pyr->dmaps[i].dists[y * w + x] == DBL_MAX) {

#if 1
			    double dist;
			    vec_t sample = em_create_query_pt(bpyr->imgs[i], em, x, y);
			    vec_t nn = query_ann_tree_pt(atree_full, sample, eps, &dist, 0, 1);
				
			    double old_dist = dist_pyr->dmaps[i].dists[y * w + x];

			    // printf("Filling in (%d, %d)\n", x, y);

			    if (dist < old_dist) {
				dist_pyr->dmaps[i].nns[y * w + x] = v2_new(Vx(nn), Vy(nn));
				dist_pyr->dmaps[i].dists[y * w + x] = dist;
			    }
			    
			    vec_free(sample);
			    vec_free(nn);			    
#else
			    dist_pyr->dmaps[i].dists[y * w + x] = 0.0;
			    dist_pyr->dmaps[i].nns[y * w + x] = v2_new(x, y);
#endif

			}
		    }
		}
	    }

	    free_ann_tree(atree_full);
	}

#if 1
	/* Check that we updated all valid entries */
	for (y = rad; y < h - rad; y++) {
	    for (x = rad; x < w - rad; x++) {
		v2_t nn = dist_pyr->dmaps[i].nns[y * w + x];

		nn = v2_sub(nn, apyr->imgs[i]->origin);
		    
		if (img_region_is_valid(bpyr->imgs[i], x - rad, x + rad, y - rad, y + rad)) {
		    if (Vx(nn) < rad || Vy(nn) < rad ||
			Vx(nn) > w_other - rad - 1 || Vy(nn) > h_other - rad - 1)
			printf("[img_create_hdp] Error: nn [%0.3f, %0.3f] out of range at (%d, %d)\n", 
			       Vx(nn), Vy(nn), x, y);
		    
		    if (dist_pyr->dmaps[i].dists[y * w + x] == DBL_MAX)
			printf("[img_create_hdp] Error: (%d, %d)!\n", x, y);
		}
	    }
	}
#endif

    }

    img_pyr_free(apyr);
    img_pyr_free(bpyr);

    return dist_pyr;
}


/* Draw A using neighborhoods from B
 * a, b: input images
 * rgb: true if the result should be rgb, false if grayscale
 * t: interpolation parameter -- 0 = first image, 1 = second image
 * rad: the radius of the neighborhoods used */
img_t *img_pyramid_remap(img_dist_pyr_t *dist_pyr, int rgb, int rad) 
{
    int x, y, w = dist_pyr->w, h = dist_pyr->h, i, j;
    img_pyr_t *pyr_out = img_pyramid_new(w, h);
    img_t *img_out;

    pyr_out->num_levels = dist_pyr->num_levels;

    for (i = 0; i < dist_pyr->num_levels; i++) {
	int num_entries = (w >> i) * (h >> i);
	int *pixel_counts = malloc(sizeof(int) * num_entries);

	pyr_out->imgs[i] = img_new(w >> i, h >> i);

	for (j = 0; j < num_entries; j++)
	    pixel_counts[j] = 0;

	/* Infer the image from the points */
	for (y = rad + 0; y < (h >> i) - rad; y += 1) { //2 * rad + 1) {
	    for (x = rad + 0; x < (w >> i) - rad; x += 1) { // * rad + 1) {
		int idx = y * (w >> i) + x;

		/* Smooth the neighborhood into the image */
		int nx, ny;
		for (ny = -rad; ny <= rad; ny++) {
		    for (nx = -rad; nx <= rad; nx++) {
			int vp = (rgb ? 3 : 1) * (((ny + rad) / 2) * (2 * rad + 1) + ((nx + rad) / 2)) + 2;
			int pidx = (y + ny) * (w >> i) + (x + nx);
			int xnn = Vx(dist_pyr->dmaps[i].nns[idx]) - (x + nx); // (w >> (i + 1));
			int ynn = Vy(dist_pyr->dmaps[i].nns[idx]) - (y + ny); // (h >> (i + 1));
			double dist = sqrt(xnn * xnn + ynn * ynn);
			double angle = 255.0 * (atan2(ynn, xnn) + M_PI) / (2 * M_PI);

#if 0
			double r = Vn(dist_pyr->nns[i][idx], vp + 0);
			double g = Vn(dist_pyr->nns[i][idx], vp + 1);
			double b = Vn(dist_pyr->nns[i][idx], vp + 2);
#else
			double r = 255 * dist / ((w >> (i+1)) + (h >> (i+1)));
			double g = 255 * dist / ((w >> (i+1)) + (h >> (i+1)));
			double b = 255 * dist / ((w >> (i+1)) + (h >> (i+1)));
#endif

			pixel_combine(pyr_out->imgs[i], x + nx, y + ny, r, g, b, pixel_counts[pidx]+1);
			pixel_counts[pidx]++;
		    }
		}
	    }
	}
	free(pixel_counts);
    }

    img_out = img_pyr_render(pyr_out);
    img_pyr_free(pyr_out);

    return img_out;
}

/* Shrink-wrap an image distance pyramid, making each level as small
 * as possible */
img_dist_pyr_t *img_dist_pyr_shrink_wrap(img_dist_pyr_t *pyr) {
    img_dist_pyr_t *pyr_new;
    int i;
    
    pyr_new = (img_dist_pyr_t *)malloc(sizeof(img_dist_pyr_t));
    
    pyr_new->em = pyr->em;
    pyr_new->nhood_radius = pyr->nhood_radius;
    pyr_new->zweight = pyr->zweight;
    pyr_new->num_levels = pyr->num_levels;

    pyr_new->dmaps = (img_dmap_t *)malloc(sizeof(img_dmap_t) * pyr_new->num_levels);

    /* Shrink-wrap all the dmaps */
    for (i = 0; i < pyr_new->num_levels; i++) {
	img_dmap_shrink_wrap(&(pyr->dmaps[i]), &(pyr_new->dmaps[i]));
    }

    pyr_new->w = pyr_new->dmaps[0].w;
    pyr_new->h = pyr_new->dmaps[0].h;

    return pyr_new;
}

