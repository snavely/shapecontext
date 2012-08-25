#ifndef __dist_h__
#define __dist_h__

#ifdef __cplusplus
extern "C" {
#endif

#define INF ((u_int32_t) 0xfeedface)

#include "binimage.h"
#include "binvolume.h"
#include "vector.h"

/* Enumeration for DT channels */
typedef enum {
    DT_RED,
    DT_GREEN,
    DT_BLUE,
    DT_INTENSITY,
} dt_channel_t;

/* Distance transform */
typedef struct {
    int xmin, ymin;   /* Min x/y coords stored in the transform */
    int w, h;         /* Width, height */
    float *dist;      /* 2D array of distances */
    iv2_t *disp;      /* 3D array of displacement vectors */
} dtrans_t;

/* 3D distance transform */
typedef struct {
    int xmin, ymin, zmin;  /* Min x/y/z coords stored in the transform */
    int w, h, d;           /* Width, height, depth */
    float axis_weights[3]; /* Axis weights */
    float *dist;           /* 3D array of (squared) distances */
    iv3_t *disp;           /* 3D array of displacement vectors */
} dtrans3D_t;

/* 5D distance transform */
typedef struct {
    int xmin, ymin, rmin, gmin, bmin;  /* Min x/y/z coords stored in the transform */
    int w, h, r_w, g_w, b_w;           /* Width, height, depth */
    float axis_weights[5]; /* Axis weights */
    float *dist;           /* 5D array of (squared) distances */
    iv5_t *disp;           /* 5D array of displacement vectors */
} dtrans5D_t;

/* Compute the distance transform of a binary image */
dtrans_t *compute_DT(bimg_t *b, double d1, double d2);

/* Compute the manhattan distance transform of a binary image */
dtrans_t *compute_manhattan_DT(bimg_t *b);

/* Compute the chessboard distance transform of a binary image */
dtrans_t *compute_chessboard_DT(bimg_t *b);

/* Compute the Euclidean distance transform of a binary image */
dtrans_t *compute_euclidean_DT(bimg_t *b);

/* Compute the 3D Euclidean distance transform of a binary volume */
dtrans3D_t *compute_euclidean_DT3D(bvol_t *b);

/* Compute a 3D distance transform using Huttenlocher's technique */
dtrans3D_t *compute_euclidean_DT3D_huttenlocher(img_t *img, 
						dt_channel_t channel);

/* Compute a 3D distance transform using Huttenlocher's technique */
dtrans3D_t *compute_euclidean_DT3D_huttenlocher_z(img_t *img,
						  dt_channel_t channel, int sub);

/* Compute a 3D distance transform using Huttenlocher's technique */
dtrans3D_t *compute_euclidean_DT3D_huttenlocher_zf(float *patch, int w, int h,
						   dt_channel_t channel, int sub);

/* Compute a 5D distance transform using Huttenlocher's technique */
dtrans5D_t *compute_euclidean_DT5D_huttenlocher(img_t *img, int rsub, int gsub, int bsub);
    
/* Use the brute force method to compute the exact DT.  Used for 
 * verifying the results of faster algorithms */
dtrans_t *compute_exact_EDT(bimg_t *b);
dtrans3D_t *compute_exact_EDT3D(bvol_t *b);

/* Set the distance between each grid point in the intensity direction.
 * This can be used to adjust how how much variation in color between
 * two images effects the image alignment */
void set_color_step_size(double step);

/* Check that the distances stored in the given Euclidean DT are 
 * consistent with the displacements */
void verify_DT(dtrans_t *d);
void verify_DT3D(dtrans3D_t *d);

#ifdef __cplusplus
}
#endif

#endif
