/* shape-context.c */
/* Shape context code */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <omp.h>

#include "image.h"
#include "util.h"
#include "vector.h"

#include "shapecontext.h"

/* Compute shape contexts for every point on an image */
void compute_shape_context(img_t *img_edge, int num_rings, int num_wedges,
			   double factor, double sigma, 
                           int normalize, int take_sqrt,
			   vec_t *descriptors)
{
    int w = img_edge->w;
    int h = img_edge->h;
    
    int x, y;

    for (y = 0; y < h; y++) {
	printf(".");
	fflush(stdout);

#pragma omp parallel for
	for (x = 0; x < w; x++) {
            int idx = y * w + x;
	    compute_shape_context_pt(img_edge, x, y, num_rings, num_wedges,
				     factor, sigma, normalize, take_sqrt,
                                     descriptors[idx].p);
	}
    }
    printf("\n");
}

/* Compute a shape context given an edge confidence image and a point */
void compute_shape_context_pt(img_t *img_edge, int x_p, int y_p, 
			      int num_rings, int num_wedges,
			      double factor, double sigma, 
                              int normalize, int take_sqrt,
			      double *descriptor)
{
    /* Compute the radius of the log-polar circle */
    double radius = pow(factor, num_rings) - 1;
    int rad_int = iround(ceil(radius));
    int feature_size = num_rings * num_wedges;

    double *bins = (double *) malloc(sizeof(double) * feature_size);
    double *weights = (double *) malloc(sizeof(double) * feature_size);

    double inv_log_factor = 1.0 / log(factor);
    double inv_2pi = 1.0 / (2.0 * M_PI);
    
    int dx, dy, i;

    for (i = 0; i < feature_size; i++) {
	bins[i] = 0.0;
	weights[i] = 0.0;
    }

    for (dy = -rad_int; dy <= rad_int; dy++) {
	for (dx = -rad_int; dx <= rad_int; dx++) {
	    int xi = x_p + dx;
	    int yi = y_p + dy;

	    /* Compute the bin for this pixel */
	    double rho = sqrt(dx * dx + dy * dy);
	    double theta = atan2(dy, dx);

	    /* Get the edge confidence */
	    int conf;

	    /* Bounds checks */
	    if (xi < 0 || xi >= img_edge->w || yi < 0 || yi >= img_edge->h)
		conf = 0;
	    else 
		conf = img_get_pixel(img_edge, xi, yi).r;

	    if (theta < 0.0)
		theta += 2 * M_PI;

	    if (rho >= radius)
		continue;

	    /* Compute the floor of the bucket number */
	    double R = log(rho + 1.0) * inv_log_factor; // / log(factor);
	    double W = num_wedges * theta * inv_2pi;
	    int Rf = iround(floor(R));
	    int Wf = iround(floor(W));
		
	    double t = R - Rf;
	    double u = W - Wf;

	    double wt0 = (1.0 - t) * (1.0 - u);
	    double wt1 = (1.0 - t) * u;
	    double wt2 = t * (1.0 - u);
	    double wt3 = t * u;

	    // printf("(R, W) ==> (%0.3f, %0.3f)\n", R, W);

	    if (Rf < 0 || Rf >= num_rings || Wf < 0 || Wf >= num_wedges)
		printf("error: %d %d\n", Rf, Wf);

	    /* Use bilinear interpolation to compute bins */
	    bins[Rf * num_wedges + Wf] += conf * wt0;
	    bins[Rf * num_wedges + (Wf+1) % num_wedges] += conf * wt1;

	    weights[Rf * num_wedges + Wf] += wt0;
	    weights[Rf * num_wedges + (Wf+1) % num_wedges] += wt1;

	    if (Rf + 1 < num_rings) {
		bins[(Rf+1) * num_wedges + Wf] += conf * wt2;
		bins[(Rf+1) * num_wedges + (Wf+1) % num_wedges] += conf * wt3;

		weights[(Rf+1) * num_wedges + Wf] += wt2;
		weights[(Rf+1) * num_wedges + (Wf+1) % num_wedges] += wt3;
	    }

#if 0
	    if (weights[2] > 0.0) {
		printf("%0.3e, %0.3e\n", R, W);
		printf("%d, %d\n", Rf, Wf);
	    }
#endif
	}
    }

    /* Normalize bins */
    for (i = 0; i < feature_size; i++) {
	if (weights[i] > 0.0)
	    descriptor[i] = bins[i] / weights[i];
	else
	    descriptor[i] = 0.0;
	
	// printf("bns[%d] %0.5e\n", i, bins[i]);
	// printf("wts[%d] %0.5e\n", i, weights[i]);
	// printf("dsc[%d] %0.5e\n", i, descriptor[i]);
    }

    if (normalize) {
	/* Normalize histogram */

	double mass = 0.0;
	for (i = 0; i < feature_size; i++) {
	    mass += descriptor[i];
	}

	if (mass > 0.0) {
	    for (i = 0; i < feature_size; i++) {
		descriptor[i] /= mass;

                if (take_sqrt)
                    descriptor[i] = sqrt(descriptor[i]);
		// descriptor[i] *= 2048.0;
	    }
	}
    }
    
    /* Copy to the descriptor */
    // memcpy(descriptor, bins, sizeof(double) * feature_size);

    free(bins);
    free(weights);
}
