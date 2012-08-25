/* match.h */
/* Match features between two sets of points */

#ifndef __match_h__
#define __match_h__

// #include "dmap.h"
#include "image.h"
#include "vector.h"

typedef enum {
    DISTANCE_EUCLIDEAN,
    DISTANCE_CHI_SQUARED,
    DISTANCE_EMD,
} distance_t;

#if 1
img_dmap_t *match_shape_contexts(img_t *img_edge1, 
				 img_t *img_edge2, 
				 vec_t *v1, vec_t *v2,
				 int num_rings, int num_wedges,
				 distance_t dtype);
#else
int match_shape_contexts(img_t *img_edge1, 
                         img_t *img_edge2, 
                         vec_t *v1, vec_t *v2,
                         int num_rings, int num_wedges,
                         distance_t dtype);
#endif

#endif /* __match_h__ */
