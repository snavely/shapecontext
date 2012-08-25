/* match.c */
/* Match features between two sets of points */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "dmap-io.h"

#include "emd.h"
#include "match.h"

img_dmap_t *match_shape_contexts_euclidean(img_t *img_edge1, 
					   img_t *img_edge2, 
					   vec_t *v1, vec_t *v2) 
{
    int d = v1[0].d;
    int i;
    int x, y;
    int w1 = img_edge1->w, h1 = img_edge1->h;
    int w2 = img_edge2->w, h2 = img_edge2->h;
    int idx, count;

    /* Create a search tree */
    ANNkd_tree_t *tree;
    vec_t axis_weights;
    img_dmap_t *dmap;

    double x_mean = 0.0, y_mean = 0.0;

    dmap = img_dmap_new(w1, h1);

    axis_weights = vec_new(d);
    
    for (i = 0; i < d; i++) {
	Vn(axis_weights,i) = 1.0;
    }

    set_ann_points_visit(200);
    tree = create_ann_tree(w2 * h2, d, v2, axis_weights);

    /* Now do the search */
    idx = count = 0;
    for (y = 0; y < h1; y++) {
	for (x = 0; x < w1; x++) {
	    double dist;
	    int nn_idx = query_pri_ann_tree_idx(tree, v1[idx], 0.0, &dist);

	    int y2 = nn_idx / w2;
	    int x2 = nn_idx % w2;

	    // printf("(%d,%d) ==> (%d,%d,%d) = %0.3f\n", x, y, nn_idx, x2, y2, dist);
	    if (x == 0 && y == 0) {
		dmap->dists[idx] = DBL_MAX;
	    } else {
		dmap->nns[idx] = v2_new(x2, y2);
		dmap->dists[idx] = dist;

		x_mean += x2 - x;
		y_mean += y2 - y;

		count++;
	    }
	    
	    idx++;
	}
    }

    x_mean /= count;
    y_mean /= count;

    printf("Mean motion: %0.3f, %0.3f\n", x_mean, y_mean);

    return dmap;    
}

static float emd_compare_features(feature_t *f1, feature_t *f2) 
{
    float dtheta = (float) (f2->theta - f1->theta);
    float drho = (float) f2->rho - f1->rho;

    /* Take care of weirdness */
    while (dtheta >= M_PI)
	dtheta -= M_PI;
    
    while (dtheta <= -M_PI)
	dtheta += M_PI;

    dtheta /= M_PI;

    /* Or could just take |cos| of difference? */

    return sqrt(dtheta * dtheta + drho * drho);
}

img_dmap_t *match_shape_contexts_emd(img_t *img_edge1, 
				     img_t *img_edge2, 
				     vec_t *v1, vec_t *v2,
				     int num_rings, int num_wedges) 
{
    int i, j, k;
    int x1, y1, x2, y2;
    int w1 = img_edge1->w, h1 = img_edge1->h;
    int w2 = img_edge2->w, h2 = img_edge2->h;
    int idx1, idx2;

    double x_mean = 0.0, y_mean = 0.0;

    feature_t *f1, *f2;
    float *weights1, *weights2;
    signature_t s1, s2;

    img_dmap_t *dmap;

    f1 = (feature_t *) malloc(sizeof(feature_t) * num_rings * num_wedges);
    f2 = (feature_t *) malloc(sizeof(feature_t) * num_rings * num_wedges);

    weights1 = (float *) malloc(sizeof(float) * num_rings * num_wedges);
    weights2 = (float *) malloc(sizeof(float) * num_rings * num_wedges);

    s1.n = num_rings * num_wedges;
    s1.Features = f1;
    s1.Weights = weights1;
    
    s2.n = num_rings * num_wedges;
    s2.Features = f2;
    s2.Weights = weights2;

    /* Fill features for signature 1 */
    for (i = 0; i < num_rings; i++) {
	for (j = 0; j < num_wedges; j++) {
	    f1[i * num_wedges + j].rho = i;
	    f1[i * num_wedges + j].theta = j * 2.0 * M_PI / num_wedges;

	    f2[i * num_wedges + j].rho = i;
	    f2[i * num_wedges + j].theta = j * 2.0 * M_PI / num_wedges;
	}
    }

    dmap = img_dmap_new(w1, h1);

    idx1 = 0;
    for (y1 = 0; y1 < h1; y1++) {
	for (x1 = 0; x1 < w1; x1++, idx1++) {
	    float min_dist = FLT_MAX;
	    int min_idx = -1;

	    printf(".");
	    fflush(stdout);

	    /* Fill the weights and weights for signature 1 */
	    for (k = 0; k < num_rings * num_wedges; k++) {
		weights1[k] = Vn(v1[idx1], k);
	    }

	    /* Match to all other features */
	    for (idx2 = 0; idx2 < w2 * h2; idx2++) {
		float dist;

		for (k = 0; k < num_rings * num_wedges; k++) {
		    weights2[k] = Vn(v2[idx2], k);
		}

		dist = emd(&s1, &s2, emd_compare_features, NULL, NULL);

		if (dist < min_dist) {
		    min_dist = dist;
		    min_idx = idx2;
		}
	    }
	    
	    x2 = min_idx % w2;
	    y2 = min_idx / w2;

	    dmap->nns[idx1] = v2_new(x2, y2);
	    dmap->dists[idx1] = min_dist;

	    x_mean += x2 - x1;
	    y_mean += y2 - y1;
	}
    }

    x_mean /= idx1;
    y_mean /= idx1;

    printf("Mean motion: %0.3f, %0.3f\n", x_mean, y_mean);

    free(f1);
    free(f2);

    free(weights1);
    free(weights2);

    return dmap;
}

img_dmap_t *match_shape_contexts(img_t *img_edge1, 
				 img_t *img_edge2, 
				 vec_t *v1, vec_t *v2,
				 int num_rings, int num_wedges,
				 distance_t dtype)
{
    switch (dtype) {
    case DISTANCE_EUCLIDEAN:
	return match_shape_contexts_euclidean(img_edge1, img_edge2, v1, v2);
    case DISTANCE_CHI_SQUARED:
	printf("Error: Chi-Squared metric has not yet been implemented\n");
	return NULL;
    case DISTANCE_EMD:
	return match_shape_contexts_emd(img_edge1, img_edge2, v1, v2, 
					num_rings, num_wedges);
    }

    return NULL;
}
