/* morph.h */
/* Routines for morphing two images to create a third */

#ifndef __morph_h__
#define __morph_h__

/* Morph two images */
img_t *img_morph(img_t *a, img_t *b);

/* Incrementally compute the morph between img_a and img_b */
img_t *img_morph_incremental(img_t *img_a, img_t *img_b);

#endif /* __morph_h__ */
