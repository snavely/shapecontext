#ifndef __morph2_h__
#define __morph2_h__

#include "image.h"
#include "pyramid.h"

/* Morph two images */
img_t *img_morph2(img_t *a, img_t *b);

img_t **img_morph3(img_t *a, img_t *b, int diam,
		   img_dmap_t *apyr, img_dmap_t *bpyr,
		   double zweight, int rgb, 
		   double tmin, double tmax, int steps);

img_t *img_remap(img_t *a, img_t *b, 
		 double *bdist, vec_t *bnn, int rgb);

#endif /* __morph2_h__ */
