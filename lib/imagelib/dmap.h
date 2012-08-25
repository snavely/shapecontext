/* dmap.h */

#ifndef __dmap_h__
#define __dmap_h__

#ifdef __cplusplus
extern "C" {
#endif

#include <sys/types.h>

#include "error_metric.h"
#include "image.h"
#include "vector.h"

typedef struct {
    u_int16_t w, h;  /* Width, height of the distance map */
    double *dists;   /* Matrix of distances */
    v2_t *nns;       /* Matrix of nearest neighbors */
    iv2_t *uppers;   /* Pointer to coordinates on upper level */
} img_dmap_t;

/* Create a map of distances between b and a */
img_dmap_t *img_dmap_create(img_t *a, img_t *b, error_metric_t em, double eps, int full, int use_sqrt);

/* Allocate uppers for a dmap */
void img_dmap_allocate_uppers(img_dmap_t *map);

/* Return the sub map indicated by the given bounds */
img_dmap_t *img_dmap_sub_map(img_dmap_t *map, int xmin, int ymin, int w, int h);

/* Shift the nearest neighbors in a dmap by the amoung specified */
img_dmap_t *img_dmap_shift(img_dmap_t *map, int x_shift, int y_shift);

/* Shrink wrap the given dmap, making it as small as possible */
void img_dmap_shrink_wrap(img_dmap_t *map, img_dmap_t *map_out);

#ifdef __cplusplus
}
#endif

#endif /* __dmap_h__ */
