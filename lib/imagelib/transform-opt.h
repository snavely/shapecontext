/* transform-opt.h */
/* Methods for finding optimal transforms */

#ifndef __transform_opt_h__
#define __transform_opt_h__

/* List of known transformation classes */
typedef enum {
    TRANSFORM_TRANSLATE,
    TRANSFORM_TRANSLATE_ROTATE,
    TRANSFORM_AFFINE,
    TRANSFORM_HOMOGRAPHY
} transform_class_t;

/* Find the optimal transform of a given class, given point
 * correspondences */
void find_optimal_transform(transform_class_t tclass, int n,
			    v3_t *r_pts, v3_t *l_pts, double *Tout);

#endif /* __transform_opt_h__ */
