/* morph.c */
/* Routines for morphing two images to create a third */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>

#include "anniface.h"
#include "anntool.h"
#include "bmp.h"
#include "dmap.h"
#include "dmap-io.h"
#include "error.h"
#include "image.h"

#define ANN_EPS 0.0
#define ANN_Z_WEIGHT 6.5

#define MORPH_GS
// #define MORPH_RGB

#if 0
img_t *img_morph(img_t *img_a, img_t *img_b) {
    img_t *m = NULL;
    int w, h;
    ANNkd_tree_t *atree, *btree;
    int x, y;

#ifdef MORPH_GS
    int z, d = 3;
#endif

#ifdef MORPH_RGB
    int r, g, b, d = 5;
#endif

    vec_t pt;

    if (img_a->w != img_b->w || img_a->h != img_b->h) {
        printf("Image dimensions must match to perform morph\n");
        return NULL;
    }
    
    w = img_a->w, h = img_a->h;
    m = img_new(w, h);
    
    set_ann_z_weight(ANN_Z_WEIGHT);

    /* Create the trees for each image */
#ifdef MORPH_GS
    atree = img_create_kd_tree_grayscale(img_a);
    btree = img_create_kd_tree_grayscale(img_b);
#endif

#ifdef MORPH_RGB    
    atree = img_create_kd_tree_rgb(img_a);
    btree = img_create_kd_tree_rgb(img_b);
#endif

    pt = vec_new(d);
    
    for (y = 0; y < h; y++) {
        for (x = 0; x < w; x++) {
            double min = DBL_MAX;

#ifdef MORPH_GS
            int zmin = 0;
#endif

#ifdef MORPH_RGB
            int rmin = 0, gmin = 0, bmin = 0;
#endif

            // printf("[Morph] Generating pixel (%d, %d)...\n", x, y);

            Vx(pt) = x;
            Vy(pt) = y;

#ifdef MORPH_GS
            for (z = 0; z < 256; z++) {
                double adist, bdist;

                Vz(pt) = z;

                adist = query_ann_tree_dist(atree, pt, ANN_EPS, NULL);
                bdist = query_ann_tree_dist(btree, pt, ANN_EPS, NULL);

                if (adist + bdist <= min) {
                    min = adist + bdist;
                    zmin = z;
                }
            }

            img_set_pixel(m, x, y, zmin, zmin, zmin);
#endif

#ifdef MORPH_RGB
            /* Try all 16M possible colors to find the nearest point */
            for (r = 0; r < 256; r++) {
                Vn(pt, 2) = r;
                for (g = 0; g < 256; g++) {
                    Vn(pt, 3) = g;
                    for (b = 0; b < 256; b++) {
                        double adist, bdist;

                        Vn(pt, 4) = b;
                        adist = query_ann_tree_dist(atree, pt, ANN_EPS);
                        bdist = query_ann_tree_dist(btree, pt, ANN_EPS);

                        if (adist + bdist < min) {
                            min = adist + bdist;
                            rmin = r;
                            gmin = g;
                            bmin = b;
                        }
                    }
                }
            }

            img_set_pixel(m, x, y, rmin, gmin, bmin);
#endif

        }
    }

    vec_free(pt);

    return m;
}
#endif

// #define OUTPUT_INTERMEDIATE

/* Incrementally compute the morph between img_a and img_b */
img_t *img_morph_incremental(img_t *img_a, img_t *img_b) {
    int w, h, x, y, i;
    img_t *m = NULL, *m_curr, *m_next;
    ANNkd_tree_t *atree, *btree;
    vec_t pt, grad;
    int d = 5;
    int round = 0;

    if (img_a->w != img_b->w || img_a->h != img_b->h) {
	printf("Image dimensions must match to perform morph\n");
	return NULL;
    }
    
    w = img_a->w, h = img_b->h;
    m = img_new(w, h);

    /* Come up with a guess for the solution (in this case the
     * pixel-wise average) */

    for (y = 0; y < h; y++) {
	for (x = 0; x < w; x++) {
	    color_t c_a = img_get_pixel(img_a, x, y),
                    c_b = img_get_pixel(img_b, x, y);

	    color_t c = { rint(0.5 * ((double) (c_a.r + c_b.r))),
			  rint(0.5 * ((double) (c_a.g + c_b.g))),
			  rint(0.5 * ((double) (c_a.b + c_b.b))) };

	    img_set_pixel(m, x, y, c.r, c.g, c.b);
	}
    }

    atree = img_create_kd_tree_rgb(img_a, 1);
    btree = img_create_kd_tree_rgb(img_b, 1);
    pt = vec_new(d);
    grad = vec_new(d - 2);

    m_curr = m;

#define MAX_ROUNDS 255

    /* Incrementally fix the solution */
    do {

#ifdef OUTPUT_INTERMEDIATE
	bmp_t *bout;
	FILE *f;
	char fname[32];
#endif

	img_t *m_tmp;

	double total_error = 0.0;

	m_next = img_new(w, h);

	for (y = 0; y < h; y++) {
	    Vy(pt) = y;
	    for (x = 0; x < w; x++) {
		color_t c = img_get_pixel(m_curr, x, y);
		double adist[2], bdist[2];

		Vx(pt) = x;
		Vn(pt, 2) = c.r;
		Vn(pt, 3) = c.g;
		Vn(pt, 4) = c.b;

		/* Compute all the partial derivatives */
		for (i = 2; i < d; i++) {
		    Vn(pt, i) -= 1;
		    adist[0] = query_ann_tree_dist(atree, pt, ANN_EPS, NULL);
		    bdist[0] = query_ann_tree_dist(btree, pt, ANN_EPS, NULL);
		    Vn(pt, i) += 2;
		    adist[1] = query_ann_tree_dist(atree, pt, ANN_EPS, NULL);
		    bdist[1] = query_ann_tree_dist(btree, pt, ANN_EPS, NULL);
		    Vn(pt, i) -= 1;

		    Vn(grad, i-2) = (adist[1] + bdist[1] - (adist[0] + bdist[0]));
		}

		// if (vec_norm(grad) > 0)
		//     printf("|grad| = %0.3e\n", vec_norm(grad));

		total_error += vec_norm(grad);
		
#define ERROR_EPS 1.0
		/* Take a step in the direction of the gradient */
		if (Vn(grad, 0) > ERROR_EPS && c.r < 255)
		    c.r += 1;
		else if (Vn(grad, 0) < -ERROR_EPS && c.r > 0)
		    c.r -= 1;

		if (Vn(grad, 1) > ERROR_EPS && c.g < 255)
		    c.g += 1;
		else if (Vn(grad, 1) < -ERROR_EPS && c.g > 0)
		    c.g -= 1;

		if (Vn(grad, 2) > ERROR_EPS && c.b < 255)
		    c.b += 1;
		else if (Vn(grad, 2) < -ERROR_EPS && c.b > 0)
		    c.b -= 1;

		img_set_pixel(m_next, x, y, c.r, c.g, c.b);
	    }
	}

#ifdef OUTPUT_INTERMEDIATE
	bout = img2bmp(m_curr);
	sprintf(fname, "morph%03d.bmp", round);
	f = fopen(fname, "w");
	write_bmp(f, bout);
	fclose(f);
	free_bmp(bout);
#endif	
    
	m_tmp = m_curr;
	m_curr = m_next;
	img_free(m_tmp);

	printf("Error: %0.3e\n", total_error);

	round++;
    } while (round < MAX_ROUNDS);

    vec_free(pt);
    vec_free(grad);

    return m_curr;
}

img_t *img_morph_gradient_search(img_t *img_a, img_t *img_b) {
    int x, y, w, h, k, rad, round;
    img_t *m = NULL, *m_curr = NULL, *m_next = NULL;
    ANNkd_tree_t *atree, *btree;

    if (img_a->w != img_b->w || img_a->h != img_b->h) {
	printf("Image dimensions must match to perform morph\n");
	return NULL;
    }

    w = img_a->w, h = img_b->h;

    /* Create an initial guess using pixel-wise averages */
    for (y = 0; y < h; y++) {
	for (x = 0; x < w; x++) {
	    color_t c_a = img_get_pixel(img_a, x, y),
                    c_b = img_get_pixel(img_b, x, y);

	    color_t c = { rint(0.5 * ((double) (c_a.r + c_b.r))),
			  rint(0.5 * ((double) (c_a.g + c_b.g))),
			  rint(0.5 * ((double) (c_a.b + c_b.b))) };

	    img_set_pixel(m, x, y, c.r, c.g, c.b);
	}
    }

    atree = img_create_ann_tree(img_a, EM_NHOOD_GS, &rad);
    btree = img_create_ann_tree(img_b, EM_NHOOD_GS, &rad);

    m_curr = m;

    for (round = 0; round < MAX_ROUNDS; round++) {
	m_next = img_copy(m_curr);

	for (y = 0; y < h; y++) {
	    for (x = 0; x < w; x++) {
		color_t col = img_get_pixel(m_curr, x, y);
		int c = rint(((double) col.r + col.g + col.b) / 3.0);

		/* Create the query point */
		vec_t p = em_create_query_pt(m_curr, EM_NHOOD_GS, x, y);

		/* Find the distance */
		double dist = query_ann_tree_dist(atree, p, 0.0, NULL) +
		    query_ann_tree_dist(btree, p, 0.0, NULL);

		double min_dist = dist;
		int min_hop = 0;
		int dir;
		int end = 0;

		if (c == 0) /* Can only check higher values */
		    dir = 1;
		else if (c == 255) /* Can only check lower values */
		    dir = -1;
		else {
		    /* Compute the gradient */
		    double dist_up, dist_down;
		    vec_t p_up, p_down;

		    img_set_pixel(m_curr, x, y, c + 1, c + 1, c + 1);
		    p_up = em_create_query_pt(m_curr, EM_NHOOD_GS, x, y);
		    dist_up = query_ann_tree_dist(atree, p_up, 0.0, NULL) +
			query_ann_tree_dist(btree, p_up, 0.0, NULL);

		    img_set_pixel(m_curr, x, y, c - 1, c - 1, c - 1);
		    p_down = em_create_query_pt(m_curr, EM_NHOOD_GS, x, y);
		    dist_down = query_ann_tree_dist(atree, p_down, 0.0, NULL) +
			query_ann_tree_dist(btree, p_down, 0.0, NULL);

		    vec_free(p_up);
		    vec_free(p_down);

		    img_set_pixel(m_curr, x, y, c, c, c);
		    
		    if (dist_up - dist_down >= 0)
			dir = 1;
		    else
			dir = -1;
		}

		for (k = 0; k < 8; k++) {
		    int ctest = c + dir * (1 << k);
		    double dist_curr;
		    vec_t p_curr;

		    if (ctest <= 0) {
			ctest = 0;
			end = 1;
		    } else if (ctest >= 255) {
			ctest = 255;
			end = 1;
		    }

		    p_curr = em_create_query_pt(m_curr, EM_NHOOD_GS, x, y);

		    if (end)
			break;
		}
	    }
	}
    }

    free_ann_tree(atree);
    free_ann_tree(btree);

    return m;
}

