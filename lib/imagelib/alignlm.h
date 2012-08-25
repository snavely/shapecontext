/* alignlm.h */
/* Routines for aligning pairs of images */

#ifndef __alignlm_h__
#define __alignlm_h__

#include "transform.h"
#include "transform3D.h"

/* Find a linear transformation that aligns image B to image A */
trans2D_t *align_bin_image(bimg_t *a, bimg_t *b);

/* Find a linear transformation that aligns volume B to volume A */
trans2D_t *align_bin_volume_LM(bvol_t *a, bvol_t *b);
trans2D_t *align_bin_volume_LM2(bvol_t *a, bvol_t *b);

#endif /* __alignlm_h__ */
