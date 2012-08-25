#ifndef __distutil_h__
#define __distutil_h__

#ifdef __cplusplus
extern "C" {
#endif

#include <sys/types.h>

#include "dist.h"

/* Create new DTs in which every pixel has distance infinity */
dtrans_t *new_empty_DT(int xmin, int ymin, int w, int h);
dtrans3D_t *new_empty_DT3D(int xmin, int ymin, int zmin, int w, int h, int d);
dtrans3D_t *new_empty_DT3D_uninit(int xmin, int ymin, int zmin, int w, int h, int d);
dtrans5D_t *new_empty_DT5D(int xmin, int ymin, int rmin, int gmin, int bmin, int w, int h, int rs, int bs, int gs);
    
/* Create a new DT in which every object pixel has distance 0 and every 
 * other pixel has distance infinity */
dtrans_t *new_DT(bimg_t *b);
dtrans3D_t *new_DT3D(bvol_t *b);

/* Get the value of the distance (or displacement) transform at 
 * the given coordinates */
float getdist(dtrans_t *d, int x, int y);
float getdist3D(dtrans3D_t *d, int x, int y, int z);
#define getdist3D_unsafe(d, x, y, z) ((d)->dist[(((z) - (d)->zmin) * (d)->h + ((y) - (d)->ymin)) * (d)->w + ((x) - (d)->xmin)])
float getdist5D(dtrans5D_t *d, int x, int y, int r, int g, int b);
#define getdist5D_unsafe(d, x, y, r, g, b) \
  ((d)->dist[(((y) - (d)->ymin) * (d)->w + ((x) - (d)->xmin)) * (d)->r_w * (d)->g_w * (d)->b_w + \
	     (((r) - (d)->rmin) * (d)->g_w + ((g) - (d)->gmin)) * (d)->b_w + ((b) - (d)->bmin)])

iv2_t getdisp(dtrans_t *d, int x, int y);
iv3_t getdisp3D(dtrans3D_t *d, int x, int y, int z);
iv5_t getdisp5D(dtrans5D_t *d, int x, int y, int r, int g, int b);

/* Set the value of the distance transform at the given 
 * coordinates */
void setdist(dtrans_t *d, int x, int y, float dist);
void setdist3D(dtrans3D_t *d, int x, int y, int z, float dist);
#define setdist3D_unsafe(d, x, y, z, s) (d)->dist[(((z) - (d)->zmin) * (d)->h + ((y) - (d)->ymin)) * (d)->w + ((x) - (d)->xmin)] = (s)
void setdist5D(dtrans5D_t *d, int x, int y, int r, int g, int b, float dist);
#define setdist5D_unsafe(d, x, y, r, g, b, s) \
  ((d)->dist[(((y) - (d)->ymin) * (d)->w + ((x) - (d)->xmin)) * (d)->r_w * (d)->g_w * (d)->b_w + \
	     (((r) - (d)->rmin) * (d)->g_w + ((g) - (d)->gmin)) * (d)->b_w + ((b) - (d)->bmin)] = (s))

void setdisp(dtrans_t *d, int x, int y, int16_t dx, int16_t dy);
void setdisp3D(dtrans3D_t *d, int x, int y, int z, int16_t dx, int16_t dy, int16_t dz);
void setdisp5D(dtrans5D_t *d, int x, int y, int r, int g, int b, 
	       int16_t dx, int16_t dy, int16_t dr, int16_t dg, int16_t db);

/* Returns true is the given value represents infinity */
int is_inf_disp(int16_t d);

/* Print an ascii representation of the given distance transform */
void print_DT(dtrans_t *d);

/* Compute the difference (squared) of two transforms */
dtrans_t *diff_DT(dtrans_t *d1, dtrans_t *d2);
dtrans3D_t *diff_DT3D(dtrans3D_t *d1, dtrans3D_t *d2);

/* Add up all the distances in the transform */
float sum_DT(dtrans_t *d);
float sum_DT3D(dtrans3D_t *d);

/* Free the given 2D/3D distance transform */
void free_DT(dtrans_t *d);
void free_DT3D(dtrans3D_t *d);
void free_DT5D(dtrans5D_t *d);
    
/* Read a 3D distance transform from a file */
dtrans3D_t *DT3D_read(char *filename);
/* Write a 3D distance transform to a file */
void DT3D_write(dtrans3D_t *dt, char *filename);

/* Smooth a given distance transform */
void DT3D_smooth(dtrans3D_t *dt, float sigma);

/* Render a slice of a distance transform as an image */
img_t *DT3D_render(dtrans3D_t *d, int z);

/* Compute the minimum energy image given a 3D DT */
img_t *DT3D_min_image(dtrans3D_t *dt);

#ifdef __cplusplus
}
#endif

#endif /* __distutil_h__ */
