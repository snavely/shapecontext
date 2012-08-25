/* alignDT.h */
/* Align two images using ICP and distance transforms */


#ifndef __alignDT_h__
#define __alignDT_h__

#include "image.h"
#include "transform-opt.h"

/* Find a linear transformation that aligns image B to image A.
 * Uses ICP and resamples the image after each iteration.  The
 * transform that describes the initial guess is given. */
void align_Timage_ICP_DT(img_t *a, img_t *b, double color_step_size, 
			 transform_class_t tclass, double *Tin, double *Tout,
			 int symmetric, int verbose);

void align_Timage_ICP_DT2(img_t *a, img_t *b, double color_step_size, 
			  transform_class_t tclass, double *Tin, double *Tout, 
			  int symmetric, int verbose);

#endif /* __alignDT_h__ */
