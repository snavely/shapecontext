#ifndef __metric_h__
#define __metric_h__

#include <sys/types.h>
#include "binimage.h"
#include "transform.h"
#include "transform3D.h"

/* An image is an array of intensities with a width and a height */
typedef struct {
    int w, h;
    u_int32_t *data;
} image_t;

/* Compute the L2 norm of two images */
double L2Norm(bimg_t *a, trans2D_t *T, bimg_t *b);
// u_int32_t L2Norm(image_t *a, trans2D_t *T, image_t *b);

/* Compute the chamfer distance between two images a and b,
 * given the distance transform d for a */
double chamfer_distance(bimg_t *a, bimg_t *b, dtrans_t *d);
double chamfer_distance_3D(bvol_t *a, bvol_t *b, dtrans3D_t *d);

/* Compute the chamfer distance between two images a and the image of
 * b under T, given the distance transform d for a */
double transform_chamfer_distance(bimg_t *a, bimg_t *b, trans2D_t *T, dtrans_t *d);
double transform_chamfer_distance_3D(bvol_t *a, bvol_t *b, trans3D_t *T, dtrans3D_t *d);

/* Computes the chamfer distance of a 3D volume given a 2D transform */
double transform_chamfer_distance_2D_3D(bvol_t *a, bvol_t *b, trans2D_t *T, dtrans3D_t *d);

/* Approximate the chamfer distance for point 'in' by interpolating
 * from the eight nearest gridpoints */
double point_distance_3D(double in[3], dtrans3D_t *d);

#endif /* __metric_h__ */
