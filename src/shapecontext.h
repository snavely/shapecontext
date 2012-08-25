/* shape-context.h */

#ifndef __shape_context_h__
#define __shape_context_h__

#ifdef __cplusplus
extern "C" {
#endif

#include "image.h"

/* Compute shape contexts for every point on an image */
void compute_shape_context(img_t *img_edge, int num_rings, int num_wedges,
			   double factor, double sigma, int normalize, 
			   vec_t *descriptors);

/* Compute a shape context given an edge confidence image and a point */
void compute_shape_context_pt(img_t *img_edge, int x_p, int y_p, 
			      int num_rings, int num_wedges,
			      double factor, double sigma, int normalize, 
			      double *descriptor);

#ifdef __cplusplus
}
#endif

#endif /* __shape_context_h__ */
