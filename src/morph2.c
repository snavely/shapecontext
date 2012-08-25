/* morph2.c */
/* Compute the morph of two images */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "anniface.h"
#include "defines.h"
#include "error.h"
#include "image.h"
#include "morph2.h"
#include "pyramid.h"

img_t *img_morph2(img_t *a, img_t *b) {
    ANNkd_tree_t *atree, *btree;  /* NN trees for a, b */
    int arad, brad;
    vec_t *img_pts;
    int x, y, w, h, i, j, count, idx;
    img_t *img_out;
    double eps = 0.0;

    if (a->w != b->w || a->h != b->h) {
        printf("[img_morph2] Size mismatch in images.\n");
        return NULL;
    }

    w = a->w; h = a->h;

    /* Compute neighborhood trees for a and b */
    // set_ann_z_weight(0.35);
    em_set_nhood_size(5);

    atree = img_create_ann_tree(a, EM_NHOOD_GS, &arad);
    btree = img_create_ann_tree(b, EM_NHOOD_GS, &brad);
    
    img_out = img_new(w, h);

    /* Initialize the set of image points */
    img_pts = (vec_t *) malloc(sizeof(vec_t) * w * h);

#if 1
    count = 0;
    for (y = 0; y < h; y++) {
        for (x = 0; x < w; x++) {
            img_pts[count] = vec_new(27); /* 27 dimensions */
            Vx(img_pts[count]) = x;
            Vy(img_pts[count]) = y;
            count++;
        }
    }
#endif

    /* Find a good neighborhood for each pixel */
    idx = 0;

    for (y = 2; y < h - 2; y += 5) {
	// printf(".");
	// fflush(stdout);

        for (x = 2; x < w - 2; x += 5) {
            double min_dist = DBL_MAX;
            int imin = 0;
            double adist, bdist;

	    Vx(img_pts[idx]) = x;
	    Vy(img_pts[idx]) = y;

            /* Start off by coloring all pixels in the neighborhood
             * the same color */
            for (i = 0; i < 256; i++) {
                for (j = 2; j < 27; j++)
                    Vn(img_pts[idx], j) = i;

                adist = query_ann_tree_dist(atree, img_pts[idx], eps, NULL);
                bdist = query_ann_tree_dist(btree, img_pts[idx], eps, NULL);
                
                if (adist + bdist < min_dist) {
                    imin = i;
                    min_dist = adist + bdist;
                }

		printf("%d,%d,%d: %0.3f\n", x, y, i, adist + bdist);
            }

            /* Set the point back to the minimum */
            for (j = 2; j < 27; j++)
                Vn(img_pts[idx], j) = imin;

            /* Set the pixels in the image to the chosen color */
	    for (j = y - 2; j <= y + 2; j++)
		for (i = x - 2; i <= x + 2; i++)
		    img_set_pixel(img_out, i, j, imin, imin, imin);

	    idx++;
        }
    }

    printf("\n");

    return img_out;
}

/* Combine num_parts - 1 parts dest with 1 part src, 
 * storing the result back in dest */
void vec_average(vec_t dest, vec_t src, double old_weight, double new_weight) {
    double ratio0 = old_weight / new_weight;
    double ratio1 = 1.0 - ratio0;
    int i;

    if (dest.d != src.d) {
	printf("[vec_average] Mismatch in vector sizes\n");
	return;
    }
    
    for (i = 0; i < dest.d; i++) {
	Vn(dest, i) = ratio0 * Vn(dest, i) + ratio1 * Vn(src, i);
    }
}

void wtdave(double *d, double val, double old_weight, double new_weight) {
    double ratio0 = old_weight / new_weight;
    double ratio1 = 1.0 - ratio0;
    
    *d = ratio0 * (*d) + ratio1 * val;
}

/* Morph two image by averaging nearby points
 *   a, b : input images
 *   rgb  : true if the result should be rgb, false if grayscale
 *   t    : interpolation parameter -- 0 = first image, 1 = second image */
img_t **img_morph3(img_t *a, img_t *b, int diam,
		   img_dmap_t *amap, img_dmap_t *bmap,
		   double zweight, int rgb, 
		   double tmin, double tmax, int steps) 
{
    // int diam = 9; /* HACK ALERT */
    int w, h, x, y, xidx, yidx, idx;
    img_t **img_out;
    vec_t *img_pts;
    double *mapped_pts;
    double *pixel_map;
    double *pixel_counts;
    int i, k, dim = 2 + diam * diam * (rgb ? 3 : 1);
    int rad, xnn, ynn, xnn_back, ynn_back;
    error_metric_t em = rgb ? EM_NHOOD_RGB : EM_NHOOD_GS;
    double dist; // eps = 0.0
    // double aminvar, amaxvar, bminvar, bmaxvar, var;
    double var;
    double amindist, amaxdist, bmindist, bmaxdist, a_avedist, b_avedist;
    double t, dt = (tmax - tmin) / ((double) steps - 1);
    int step, num_a_samples, num_b_samples;
    vec_t sample, nn, nn_back, pt, *a_samples, *b_samples;
    v2_t back_disp;

    double xdisp, ydisp;

    /* svar0 == 2 :: combine samples from both images
     * svar0 == 1, svar1 == 0 :: take samples from image A
     * svar0 == 1, svar1 == 1 :: take samples from image B */
    int svar0 = 1, svar1 = 0;

    int interpolate_color = 1;

    if (a->w != b->w || a->h != b->h) {
	printf("[img_morph3] Size mismatch in images\n");
	return NULL;
    }

    w = a->w;  h = a->h;

    em_set_nhood_size(diam);

#if 0
    /* Compute variances */
    img_find_min_max_variance(a, 3, &aminvar, &amaxvar);
    img_find_min_max_variance(b, 3, &bminvar, &bmaxvar);
#endif

    amindist = bmindist = DBL_MAX;
    amaxdist = bmaxdist = 0.0;

    /* Initialize the map of points */
    img_pts = (vec_t *) malloc(sizeof(vec_t) * 4 * w * h);
    mapped_pts = (double *) malloc(sizeof(double) * 4 * w * h);

    for (i = 0; i < 4 * w * h; i++)
	mapped_pts[i] = 0.0;

    rad = diam / 2;

    /* Compute min/max and average distances */
    a_avedist = b_avedist = 0.0;
    num_a_samples = num_b_samples = 0;

    for (y = rad; y < h - rad; y++) {
	for (x = rad; x < w - rad; x++) {
	    /* A */
	    dist = amap->dists[y * amap->w + x];

	    if (dist != DBL_MAX) {
		dist = sqrt(dist);

		if (dist < amindist)
		    amindist = dist;
		if (dist > amaxdist)
		    amaxdist = dist;

		a_avedist += dist;
		num_a_samples++;
	    }

	    /* B */
	    dist = bmap->dists[y * bmap->w + x];

	    if (dist != DBL_MAX) {
		dist = sqrt(dist);

		if (dist < bmindist)
		    bmindist = dist;
		if (dist > bmaxdist)
		    bmaxdist = dist;

		b_avedist += dist;

		num_b_samples++;
	    }
	}
    }

    a_avedist /= ((double) num_a_samples);
    b_avedist /= ((double) num_b_samples);

    printf("[img_morph3] {amin, bmin} x {amax, bmax} = {%0.3f, %0.3f} x {%0.3f, %0.3f}\n",
	   amindist, bmindist, amaxdist, bmaxdist);
    printf("[img_morph3] {a_ave, b_ave} = {%0.3f, %0.3f}\n", a_avedist, b_avedist);
    printf("[img_morph3] Tables complete\n");

    img_out = (img_t **) malloc(sizeof(img_t *) * steps);
    
    for (i = 0; i < steps; i++) {
	img_out[i] = img_new(2 * w, 2 * h);
    }
    
    pt = vec_new(dim);
    pixel_counts = (double *) malloc(sizeof(double) * (2 * w) * (2 * h));
    pixel_map = (double *) malloc(sizeof(double) * (rgb ? 3 : 1) * (2 * w) * (2 * h));

    /* Initialize the samples for each image */
    a_samples = (vec_t *) malloc(sizeof(vec_t) * w * h);
    b_samples = (vec_t *) malloc(sizeof(vec_t) * w * h);

    for (y = rad; y < h - rad; y++) {
	for (x = rad; x < w - rad; x++) {
	    a_samples[y * w + x] = em_create_query_pt(a, em, x, y);
	    b_samples[y * w + x] = em_create_query_pt(b, em, x, y);
	}
    }

    printf("[img_morph3] Samples complete\n");

    for (step = 0, t = tmin; t <= tmax; t += dt, step++) {

	/* Reset mapped points */
	for (i = 0; i < 4 * w * h; i++)
	    mapped_pts[i] = 0.0;

	for (i = 0; i < 4 * w * h; i++)
	    pixel_counts[i] = 0.0;

	/* Clear pixels */
	for (i = 0; i < (rgb ? 3 : 1) * (4 * w * h); i++)
	    pixel_map[i] = 0.0;

	k = 0;
	for (y = rad; y < h - rad; y++) {
	    for (x = rad; x < w - rad; x++, k++) {
		double weight;
		double dist, distsq;
		int orig_only = 0;

		if ((k % svar0) == svar1) {
		    /* Take even samples from a */

		    sample = a_samples[y * w + x];

		    if (amap->dists[y * amap->w + x] == DBL_MAX)
			continue;
		    
		    xnn = Vx(amap->nns[y * amap->w + x]);
		    ynn = Vy(amap->nns[y * amap->w + x]);

		    if (xnn < rad || xnn >= w - rad || ynn < rad || ynn >= h - rad) {
			orig_only = 1;
			// printf("(xnn, ynn) ==> (%d, %d)\n", xnn, ynn);
		    } else {
			if (amap->dists[y * amap->w + x] == DBL_MAX) {
			    xnn = ynn = 0;
			    nn = b_samples[0];
			    continue;
			} else {
			    nn = b_samples[ynn * w + xnn];
			}
		    
			// xnn_back = Vx(bmap->nns[ynn * bmap->w + xnn]);
			// ynn_back = Vy(bmap->nns[ynn * bmap->w + xnn]);

			// nn_back = a_samples[ynn_back * w + xnn_back];

			distsq = amap->dists[y * amap->w + x];
			dist = sqrt(distsq);
		    }
		} else {
		    /* Take odd samples from b */

		    sample = b_samples[y * w + x];

		    if (bmap->dists[y * bmap->w + x] == DBL_MAX)
			continue;

		    xnn = Vx(bmap->nns[y * bmap->w + x]);
		    ynn = Vy(bmap->nns[y * bmap->w + x]);

		    if (xnn < rad || xnn >= w - rad || ynn < rad || ynn >= h - rad) {
			orig_only = 1;
			// printf("(xnn, ynn) ==> (%d, %d)\n", xnn, ynn);
		    } else {
			if (bmap->dists[y * bmap->w + x] == DBL_MAX) {
			    xnn = ynn = 0;
			    nn = a_samples[0];
			    continue;
			} else {
			    nn = a_samples[ynn * w + xnn];
			}
		    
			// xnn_back = Vx(amap->nns[ynn * amap->w + xnn]);
			// ynn_back = Vy(amap->nns[ynn * amap->w + xnn]);

			// nn_back = b_samples[ynn_back * w + xnn_back];
			
			distsq = bmap->dists[y * bmap->w + x];
			dist = sqrt(distsq);
		    }
		}

		/* Compute weighted average of the two samples */
		for (i = 0; i < dim; i++) {
		    if ((interpolate_color || i < 2) && !orig_only) {
			if ((k % svar0) == svar1)
			    Vn(pt, i) = (1.0 - t) * Vn(sample, i) + t * Vn(nn, i);
			else
			    Vn(pt, i) = (1.0 - t) * Vn(nn, i) + t * Vn(sample, i);
		    } else {
			if (i == 0) {
			    if ((k % svar0) == svar1)
				Vn(pt, i) = (1.0 - t) * Vn(sample, i) + t * xnn;
			    else
				Vn(pt, i) = (1.0 - t) * xnn + t * Vn(sample, i);

			    Vn(pt, i) = CLAMP(Vn(pt, i), 0, w - 1);
			} else if (i == 1) {
			    if ((k % svar0) == svar1)
				Vn(pt, i) = (1.0 - t) * Vn(sample, i) + t * ynn;
			    else
				Vn(pt, i) = (1.0 - t) * ynn + t * Vn(sample, i);

			    Vn(pt, i) = CLAMP(Vn(pt, i), 0, h - 1);
			} else {

#define BG 40.0
#if 1
			if ((k % svar0) == svar1)
			    Vn(pt, i) = (1.0 - t) * Vn(sample, i) + t * BG;
			else
			    Vn(pt, i) = (1.0 - t) * BG + t * Vn(sample, i);
#else
			Vn(pt, i) = Vn(sample, i);
#endif
			}
		    }
		}

		/* Insert the point into the map */
		xidx = 2 * Vx(pt);
		yidx = 2 * Vy(pt);

#if 0
		if (distsq != DBL_MAX) {
		    if (xidx < 2 * rad || yidx < 2 * rad ||
			xidx > 2 * (w - rad - 1) || yidx > 2 * (h - rad - 1)) {
			printf("error0!\n");
		    }
		}
#endif
		
		idx = yidx * 2 * w + xidx;
	    
#if 0
		/* Average with all previous points */
		if ((k % svar0) == svar1) {
		    var = img_variance_grayscale(a, x, y, 3);
		    var = (var - aminvar) / (amaxvar - aminvar);
		} else {
		    var = img_variance_grayscale(b, x, y, 3);
		    var = (var - bminvar) / (bmaxvar - bminvar);
		}
#endif

#define ALPHA 100.0
#define BETA 1.0

#if 0
		xdisp = Vx(nn) - Vx(sample);
		ydisp = Vy(nn) - Vy(sample);

#define SIGMA 100.0
		// weight = ALPHA * (var + 0.05) + BETA * (1.0 - dist);
		if ((k % svar0) == svar1)
		    weight = exp(-dist / (SIGMA * a_avedist));
		else
		    weight = exp(-dist / (SIGMA * b_avedist));

#if 1
		if ((k % svar0) == svar1 && dist > SIGMA * a_avedist)
		    weight = 0.0;
		else if ((k % svar0) != svar1 && dist > SIGMA * b_avedist)
		    weight = 0.0;

#if 0
		var = img_variance(a, x, y, 5);

		if (var < 30.0)
		    weight = 0.0;
#endif

		/* Set the weight to zero if the correspondence does
		 * not work both ways */

#if 0
		back_disp = v2_new(Vx(nn_back) - Vx(sample),
				   Vy(nn_back) - Vy(sample));

		// #define SIMILARITY_THRESHOLD (10000.0 * 10000.0)
		// #define SIMILARITY_THRESHOLD (8.0 * 8.0)
#define SIMILARITY_THRESHOLD (10000.0 * 10000.0)
		if (v2_norm(back_disp) > SIMILARITY_THRESHOLD)
		    weight = 0.0;
#endif
#endif

		//printf("dist = %0.3f\n", dist);
		//printf("var = %0.3f\n", var);
		//printf("weight = %0.3f\n", weight);

		// if (xdisp * xdisp + ydisp * ydisp > 40 * 40)
		//   weight = 0.0;

#else
#if 0
		back_disp = v2_new(Vx(nn_back) - Vx(sample),
				   Vy(nn_back) - Vy(sample));

#define SIMILARITY_THRESHOLD (4.0 * 4.0)
		if (v2_norm(back_disp) > SIMILARITY_THRESHOLD)
		    weight = 0.0;
		else 
#endif
		    weight = 1.0;
#endif

		if (weight > 0.0) {
		    if (mapped_pts[idx] == 0.0) {
			img_pts[idx] = vec_new(dim);
			vec_copy(img_pts[idx], pt);
		    } else {
			/* Weight this based on the variance and the error */
			vec_average(img_pts[idx], pt, mapped_pts[idx], mapped_pts[idx] + weight);
		    }

		    mapped_pts[idx] += weight;
		}
		
		// vec_free(sample);
		// vec_free(nn);

		if ((k % 2) == 0)
		    x--;
	    }
	}
    
#if 1
	/* Infer the image from the points */
	for (y = 0; y < 2 * h; y++) {
	    for (x = 0; x < 2 * w; x++) {
		idx = y * 2 * w + x;
		if (mapped_pts[idx] > 0) {
		    /* Smooth the neighborhood into the image */
		    // int vp = 2; /* vector position */
		    int nx, ny;
		    for (ny = -rad; ny <= rad; ny++) {
			for (nx = -rad; nx <= rad; nx++) {
			    int vp = (rgb ? 3 : 1) * (((ny + rad) / 2) * diam + ((nx + rad) / 2)) + 2;
			    int pidx = (y + ny) * (2 * w) + (x + nx);

			    if (x + nx < 0 || x + nx >= 2 * w || y + ny < 0 || y + ny >= 2 * h)
				continue;

			    if (pidx >= (2 * w) * (2 * h))
				printf("error!\n");
			    
			    if (!rgb) {
				double val = Vn(img_pts[idx], vp);
				wtdave(pixel_map + pidx, val, 
				       pixel_counts[pidx], 
				       pixel_counts[pidx] + mapped_pts[pidx]);
			    } else {
				double r = Vn(img_pts[idx], vp + 0);
				double g = Vn(img_pts[idx], vp + 1);
				double b = Vn(img_pts[idx], vp + 2);

				wtdave(pixel_map + 3 * pidx + 0, r,
				       pixel_counts[pidx],
				       pixel_counts[pidx] + mapped_pts[idx]);
				wtdave(pixel_map + 3 * pidx + 1, g,
				       pixel_counts[pidx],
				       pixel_counts[pidx] + mapped_pts[idx]);
				wtdave(pixel_map + 3 * pidx + 2, b,
				       pixel_counts[pidx],
				       pixel_counts[pidx] + mapped_pts[idx]);
			    }
			    
			    pixel_counts[pidx] += mapped_pts[idx];
			    // if (abs(nx) % 2 == 1) vp++;
			}
		    }
		}
	    }
	}

	/* Render the image */
	for (y = 0; y < 2 * h; y++) {
	    for (x = 0; x < 2 * w; x++) {
		int pidx = y * 2 * w + x;

		if (!rgb) {
		    int val = rint(pixel_map[pidx]);
		    img_set_pixel(img_out[step], x, y, val, val, val);
		} else {
		    if (pixel_counts[pidx] > 0.0) {
			int r = rint(pixel_map[3 * pidx + 0]);
			int g = rint(pixel_map[3 * pidx + 1]);
			int b = rint(pixel_map[3 * pidx + 2]);

			img_set_pixel(img_out[step], x, y, r, g, b);
		    } else {
			img_set_pixel(img_out[step], x, y, 0, 0, 200);
		    }
		}
	    }
	}

	/* Free the image points */
	for (i = 0; i < 4 * w * h; i++) {
	    if (mapped_pts[i] > 0)
		vec_free(img_pts[i]);
	}
    }
	
#endif

#if 0
    /* Infer the image from the points */
    for (y = 0; y < h; y++) {
	for (x = 0; x < w; x++) {
	    idx = 2 * y * 2 * w + 2 * x;
	    if (mapped_pts[idx] > 0) {
		/* Smooth the neighborhood into the image */
		int vp = 2; /* vector position */
		int nx, ny;
		for (ny = -2; ny <= 2; ny++) {
		    for (nx = -2; nx <= 2; nx++) {
			int pidx = (y + ny) * w + (x + nx);

			pixel_combine(img_out, x + nx, y + ny, Vn(img_pts[idx], vp), pixel_counts[pidx]+1);
			pixel_counts[pidx]++;
			vp++;
		    }
		}
	    }
	}
    }
#endif

#if 0
    /* Infer the image from the points */
    for (y = 0; y < h; y++) {
	for (x = 0; x < w; x++) {
	    idx = 2 * y * 2 * w + 2 * x;
	    if (mapped_pts[idx] > 0) {
		/* Take the point to be the middle of the neighborhood */
		int val = rint(Vn(img_pts[idx], 14));
		
		img_set_pixel(img_out, x, y, val, val, val);
	    }
	}
    }
#endif

    for (idx = 0; idx < w * h; idx++) {
	vec_free(a_samples[idx]);
	vec_free(b_samples[idx]);
    }
    
    free(a_samples);
    free(b_samples);

    vec_free(pt);

    free(img_pts);
    free(mapped_pts);
    free(pixel_counts);
    free(pixel_map);

    return img_out;
}

#if 0
/* Morph two image by averaging nearby points
 * a, b: input images
 * rgb: true if the result should be rgb, false if grayscale
 * t: interpolation parameter -- 0 = first image, 1 = second image */
img_t **img_morph4(img_t *a, img_t *b, int diam,
		   double *adist, double *bdist,
		   v2_t *ann, v2_t *bnn,
		   double zweight, int rgb, 
		   double tmin, double tmax, int steps) 
{
    int w, h, x, y, xidx, yidx, idx;
    img_t **img_out;
    vec_t *img_pts;
    double *mapped_pts;
    double *pixel_map;
    double *pixel_counts;
    int i, k, dim = 2 + diam * diam * (rgb ? 3 : 1);
    int rad, xnn, ynn;
    error_metric_t em = rgb ? EM_NHOOD_RGB : EM_NHOOD_GS;
    double dist; // eps = 0.0
    double aminvar, amaxvar, bminvar, bmaxvar, var;
    double amindist, amaxdist, bmindist, bmaxdist, a_avedist, b_avedist;
    double t, dt = (tmax - tmin) / ((double) steps - 1);
    int step, num_samples;
    vec_t sample, nn, pt;

    double xdisp, ydisp;

    if (a->w != b->w || a->h != b->h) {
	printf("[img_morph4] Size mismatch in images\n");
	return NULL;
    }

    w = a->w;  h = a->h;

    em_set_nhood_size(diam);
    rad = diam / 2;

    /* Compute variances */
    img_find_min_max_variance(a, 3, &aminvar, &amaxvar);
    img_find_min_max_variance(b, 3, &bminvar, &bmaxvar);

    amindist = bmindist = DBL_MAX;
    amaxdist = bmaxdist = 0.0;

    /* We create four maps of points -- 
     *                I_1, I'_1, I_2, I'_2 */

    /* Initialize the map of points */
    i1a_pts = malloc(sizeof(vec_t) * w * h);  /* Points */
    i1a_con = malloc(sizeof(double) * w * h); /* Confidences */
    i1_nns  = malloc(sizeof(vec_t) * w * h);  /* NNs */

    i1b_pts = malloc(sizeof(vec_t) * 4 * w * h);  /* Points */
    i1b_con = malloc(sizeof(double) * 4 * w * h); /* Confidences */

    i2a_pts = malloc(sizeof(vec_t) * w * h);  /* Points */
    i2a_con = malloc(sizeof(double) * w * h); /* Confidences */
    i2_nns  = malloc(sizeof(vec_t) * w * h);  /* NNs */

    i2b_pts = malloc(sizeof(vec_t) * 4 * w * h);  /* Points */
    i2b_con = malloc(sizeof(double) * 4 * w * h); /* Confidences */

    /* i1a_pts and i2a_pts are constant, so we can precompute them, as
     * well as the nearest neighbors */

    for (y = rad; y < h - rad; y++) {
	for (x = rad; x < w - rad; x++) {
	    vec_t sample, nn;
	    int idx = y * w + x;

	    /* A */
	    sample = em_create_query_pt(a, em, x, y);
	    xnn = Vx(ann[idx]);
	    ynn = Vy(ann[idx]);
	    nn = em_create_query_pt(b, em, xnn, ynn);

	    i1a_pts[idx] = sample;
	    i1_nns[idx] = nn;

	    /* B */
	    sample = em_create_query_pt(b, em, x, y);
	    xnn = Vx(bnn[idx]);
	    ynn = Vy(bnn[idx]);
	    nn = em_create_query_pt(a, em, xnn, ynn);

	    i2a_pts[idx] = sample;
	    i2_nns[idx] = nn;

	    /* Confidence always 1.0 for both images */
	    i1a_con[idx] = 1.0;
	    i2a_con[idx] = 1.0;
	}
    }

    /* Compute the maximum and average distances */
    a_avedist = b_avedist = 0.0;
    num_samples = 0;

    for (y = rad; y < h - rad; y++) {
	for (x = rad; x < w - rad; x++) {
	    /* A */
	    dist = adist[y * w + x]; /* Squared distance */
	    dist = sqrt(dist);

	    if (dist < amindist)
		amindist = dist;
	    if (dist > amaxdist)
		amaxdist = dist;

	    a_avedist += dist;

	    /* B */
	    dist = bdist[y * w + x]; /* Squared distance */
	    dist = sqrt(dist);

	    if (dist < bmindist)
		bmindist = dist;
	    if (dist > bmaxdist)
		bmaxdist = dist;

	    b_avedist += dist;

	    num_samples++;
	}
    }

    a_avedist /= ((double)num_samples);
    b_avedist /= ((double)num_samples);

    printf("[img_morph4] {amin, bmin} x {amax, bmax}= {%0.3f, %0.3f} x {%0.3f, %0.3f}\n",
	   amindist, bmindist, amaxdist, bmaxdist);
    printf("[img_morph4] {a_ave, b_ave} = {%0.3f, %0.3f}\n", a_avedist, b_avedist);
    printf("[img_morph4] Tables complete\n");

    img_out = malloc(sizeof(img_t *) * steps);
    
    for (i = 0; i < steps; i++) {
	img_out[i] = img_new(2 * w, 2 * h);
    }
    
    pt = vec_new(dim);
    pixel_counts = malloc(sizeof(double) * (2 * w) * (2 * h));
    pixel_map = malloc(sizeof(double) * (rgb ? 3 : 1) * (2 * w) * (2 * h));

    for (step = 0, t = tmin; t <= tmax; t += dt, step++) {

	/* Reset the confidence maps */
	for (i = 0; i < 4 * w * h; i++) {
	    i1b_con[i] = 0.0;
	    i2b_con[i] = 0.0;
	}

	for (i = 0; i < 4 * w * h; i++)
	    pixel_counts[i] = 0.0;

	/* Clear pixels */
	for (i = 0; i < (rgb ? 3 : 1) * (4 * w * h); i++)
	    pixel_map[i] = 0.0;

	/* Compute the maps and confidences for i1b and i2b */
	for (y = rad; y < h - rad; y++) {
	    for (x = rad; x < w - rad; x++) {
		int idx = y * w + x;

		vec_t sample = i1a_pts[idx];
		vec_t nn = i1a_nns[idx];

		double dist = sqrt(adist[idx]);

		/* Compute weighted average the two samples */
		for (i = 0; i < dim; i++) {
		    Vn(pt, i) = (1.0 - t) * Vn(sample, i) + t * Vn(nn, i);
		}

		/* Compute confidence */
		confidence = exp(-dist / sigma);

		/* Insert point into the map based on the confidence */
		xidx = 2 * Vx(pt);
		yidx = 2 * Vy(pt);

		
	    }
	}

	k = 0;
	for (y = rad; y < h - rad; y++) {
	    for (x = rad; x < w - rad; x++, k++) {
		double weight;
		double dist;

		if ((k % 2) == 0) {
		    /* Take even samples from a */

		    //if (adist[y * w + x] > (amindist + zweight * 128.0))
		    //	continue;

		    sample = em_create_query_pt(a, em, x, y);
		    xnn = Vx(ann[y * w + x]);
		    ynn = Vy(ann[y * w + x]);
		    nn = em_create_query_pt(b, em, xnn, ynn);

		    dist = adist[y * w + x];
		    /* Scale between 0 and 1 */
		    dist = (dist - amindist) / (amaxdist - amindist);
		} else {
		    /* Take odd samples from b */
		    // if (bdist[y * w + x] > (bmindist + zweight * 128.0))
		    //   continue;

		    sample = em_create_query_pt(b, em, x, y);
		    xnn = Vx(bnn[y * w + x]);
		    ynn = Vy(bnn[y * w + x]);
		    nn = em_create_query_pt(a, em, xnn, ynn);

		    dist = bdist[y * w + x];
		    /* Scale between 0 and 1 */
		    dist = (dist - bmindist) / (bmaxdist - bmindist);
		}

		/* Compute weighted average the two samples */
		for (i = 0; i < dim; i++) {
		    if ((k % 2) == 0)
			Vn(pt, i) = (1.0 - t) * Vn(sample, i) + t * Vn(nn, i);
		    else
			Vn(pt, i) = (1.0 - t) * Vn(nn, i) + t * Vn(sample, i);
		}

		/* Insert the point into the map */
		xidx = 2 * Vx(pt);
		yidx = 2 * Vy(pt);

#if 0
		if (xidx < 2 * rad || yidx < 2 * rad ||
		    xidx > 2 * (w - rad - 1) || yidx > 2 * (h - rad - 1))
		    printf("error0!\n");
#endif

		idx = yidx * 2 * w + xidx;
	    
		/* Average with all previous points */
		if ((k % 2) == 0) {
		    var = img_variance_grayscale(a, x, y, 3);
		    var = (var - aminvar) / (amaxvar - aminvar);
		} else {
		    var = img_variance_grayscale(b, x, y, 3);
		    var = (var - bminvar) / (bmaxvar - bminvar);
		}

#define ALPHA 100.0
#define BETA 1.0

#if 1
		xdisp = Vx(nn) - Vx(sample);
		ydisp = Vy(nn) - Vy(sample);

		weight = ALPHA * (var + 0.05) + BETA * (1.0 - dist);
		//printf("dist = %0.3f\n", dist);
		//printf("var = %0.3f\n", var);
		//printf("weight = %0.3f\n", weight);

		// if (xdisp * xdisp + ydisp * ydisp > 40 * 40)
		//   weight = 0.0;

#else
		weight = 1.0;
#endif

		if (weight > 0.0) {
		    if (mapped_pts[idx] == 0.0) {
			img_pts[idx] = vec_new(dim);
			vec_copy(img_pts[idx], pt);
		    } else {
			/* Weight this based on the variance and the error */
			vec_average(img_pts[idx], pt, mapped_pts[idx], mapped_pts[idx] + weight);
		    }

		    mapped_pts[idx] += weight;
		}
		
		vec_free(sample);
		vec_free(nn);
	    }
	}
    
#if 1
	/* Infer the image from the points */
	for (y = 0; y < 2 * h; y++) {
	    for (x = 0; x < 2 * w; x++) {
		idx = y * 2 * w + x;
		if (mapped_pts[idx] > 0) {
		    /* Smooth the neighborhood into the image */
		    // int vp = 2; /* vector position */
		    int nx, ny;
		    for (ny = -rad; ny <= rad; ny++) {
			for (nx = -rad; nx <= rad; nx++) {
			    int vp = (rgb ? 3 : 1) * (((ny + rad) / 2) * diam + ((nx + rad) / 2)) + 2;
			    int pidx = (y + ny) * (2 * w) + (x + nx);
			    
			    if (pidx >= (2 * w) * (2 * h))
				printf("error!\n");
			    
			    if (!rgb) {
				double val = Vn(img_pts[idx], vp);
				wtdave(pixel_map + pidx, val, 
				       pixel_counts[pidx], 
				       pixel_counts[pidx] + mapped_pts[pidx]);
			    } else {
				double r = Vn(img_pts[idx], vp + 0);
				double g = Vn(img_pts[idx], vp + 1);
				double b = Vn(img_pts[idx], vp + 2);

				wtdave(pixel_map + 3 * pidx + 0, r,
				       pixel_counts[pidx],
				       pixel_counts[pidx] + mapped_pts[idx]);
				wtdave(pixel_map + 3 * pidx + 1, g,
				       pixel_counts[pidx],
				       pixel_counts[pidx] + mapped_pts[idx]);
				wtdave(pixel_map + 3 * pidx + 2, b,
				       pixel_counts[pidx],
				       pixel_counts[pidx] + mapped_pts[idx]);
			    }
			    
			    pixel_counts[pidx] += mapped_pts[idx];
			    // if (abs(nx) % 2 == 1) vp++;
			}
		    }
		}
	    }
	}

	/* Render the image */
	for (y = 0; y < 2 * h; y++) {
	    for (x = 0; x < 2 * w; x++) {
		int pidx = y * 2 * w + x;

		if (!rgb) {
		    int val = rint(pixel_map[pidx]);
		    img_set_pixel(img_out[step], x, y, val, val, val);
		} else {
		    if (pixel_counts[pidx] > 0.0) {
			int r = rint(pixel_map[3 * pidx + 0]);
			int g = rint(pixel_map[3 * pidx + 1]);
			int b = rint(pixel_map[3 * pidx + 2]);

			img_set_pixel(img_out[step], x, y, r, g, b);
		    } else {
			img_set_pixel(img_out[step], x, y, 0, 0, 200);
		    }
		}
	    }
	}

	/* Free the image points */
	for (i = 0; i < 4 * w * h; i++) {
	    if (mapped_pts[i] > 0)
		vec_free(img_pts[i]);
	}
    }
	
#endif

#if 0
    /* Infer the image from the points */
    for (y = 0; y < h; y++) {
	for (x = 0; x < w; x++) {
	    idx = 2 * y * 2 * w + 2 * x;
	    if (mapped_pts[idx] > 0) {
		/* Smooth the neighborhood into the image */
		int vp = 2; /* vector position */
		int nx, ny;
		for (ny = -2; ny <= 2; ny++) {
		    for (nx = -2; nx <= 2; nx++) {
			int pidx = (y + ny) * w + (x + nx);

			pixel_combine(img_out, x + nx, y + ny, Vn(img_pts[idx], vp), pixel_counts[pidx]+1);
			pixel_counts[pidx]++;
			vp++;
		    }
		}
	    }
	}
    }
#endif

#if 0
    /* Infer the image from the points */
    for (y = 0; y < h; y++) {
	for (x = 0; x < w; x++) {
	    idx = 2 * y * 2 * w + 2 * x;
	    if (mapped_pts[idx] > 0) {
		/* Take the point to be the middle of the neighborhood */
		int val = rint(Vn(img_pts[idx], 14));
		
		img_set_pixel(img_out, x, y, val, val, val);
	    }
	}
    }
#endif

    vec_free(pt);

    free(img_pts);
    free(mapped_pts);
    free(pixel_counts);
    free(pixel_map);

    return img_out;
}
#endif


/* Draw A using neighborhoods from B
 * a, b: input images
 * rgb: true if the result should be rgb, false if grayscale
 * t: interpolation parameter -- 0 = first image, 1 = second image */
img_t *img_remap(img_t *a, img_t *b, 
		  double *bdist, vec_t *bnn, int rgb) 
{
    int rad = 2; /* HACK ALERT: This needs to be changed when
		  * nhood_radius changes */
    int x, y, w = b->w, h = b->h, i;
    img_t *img_out = img_new(w, h);
    int *pixel_counts = (int *) malloc(sizeof(int) * w * h);

    for (i = 0; i < w * h; i++)
	pixel_counts[i] = 0;

    /* Infer the image from the points */
    for (y = rad; y < h - rad; y += 1) { //2 * rad + 1) {
	for (x = rad; x < w - rad; x += 1) { // * rad + 1) {
	    int idx = y * w + x;

	    /* Smooth the neighborhood into the image */
	    int nx, ny;
	    for (ny = -rad; ny <= rad; ny++) {
		for (nx = -rad; nx <= rad; nx++) {
		    int vp = (rgb ? 3 : 1) * (((ny + rad) / 2) * (2 * rad + 1) + ((nx + rad) / 2)) + 2;
		    int pidx = (y + ny) * w + (x + nx);
		    double r = Vn(bnn[idx], vp + 0);
		    double g = Vn(bnn[idx], vp + 1);
		    double b = Vn(bnn[idx], vp + 2);

		    pixel_combine(img_out, x + nx, y + ny, r, g, b, pixel_counts[pidx]+1);
		    pixel_counts[pidx]++;
		}
	    }
	}
    }
    
    free(pixel_counts);

    return img_out;
}
