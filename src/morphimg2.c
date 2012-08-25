/* morphimg.c */

#include <stdlib.h>

#include "anniface.h"
#include "bmp.h"
#include "error.h"
#include "fileio.h"
#include "image.h"
#include "morph.h"
#include "morph2.h"
#include "pyramid.h"
#include "pyramid-io.h"

int main(int argc, char **argv) {
    // bmp_t *b1, *b2, *bout;
    img_t *i1, *i2, **m;
    FILE *f;
    double zweight;
    int diameter;

    img_pyr_t *a_imgpyr, *b_imgpyr;

#if 1
    double tmin = 0, tmax = 1.0;
    int i, steps = 32;
#else
    double tmin = 0.5, tmax = 0.5;
    int i, steps = 1;
#endif

    img_dist_pyr_t *apyr, *bpyr;

    if (argc != 8) {
        printf("Usage: morphimg2 <zweight> <diameter> <in1.bmp> <in2.bmp> <in1.pyr> <in2.pyr> <out.bmp>\n");
        return 1;
    }

    zweight = atof(argv[1]);
    diameter = atoi(argv[2]);

#if 0
    /* Read the first bitmap */
    f = fopen(argv[3], "r");
    if (f == NULL) {
        printf("Could not open file %s for reading\n", argv[3]);
        return 1;
    }
    b1 = read_bmp(f);
    fclose(f);
    
    /* Read the second bitmap */
    f = fopen(argv[4], "r");
    if (f == NULL) {
        printf("Could not open file %s for reading\n", argv[4]);
        return 1;
    }
    b2 = read_bmp(f);
    fclose(f);

    /* Convert the bitmaps to images */
    if (b1 == NULL || b2 == NULL) {
        printf("Error reading bitmaps\n");
        return 1;
    }

    i1 = bmp2img(b1), i2 = bmp2img(b2);
    
    if (i1 == NULL || i2 == NULL) {
        printf("Error in bmp2img conversion\n");
        return 1;
    }
#endif

    i1 = img_read_bmp_file(argv[3]);
    i2 = img_read_bmp_file(argv[4]);

    /* Read the first map */
    f = open_file(argv[5], "r");
    apyr = img_read_distance_pyramid(f);
    fclose(f);
    
    /* Read the second map */
    f = open_file(argv[6], "r");
    bpyr = img_read_distance_pyramid(f);
    fclose(f);

    set_ann_z_weight(zweight);

    a_imgpyr = img_create_gaussian_pyramid(i1, 0);
    b_imgpyr = img_create_gaussian_pyramid(i2, 0);

#if 1
    m = img_morph3(i1, i2, diameter, &(apyr->dmaps[0]), &(bpyr->dmaps[0]), zweight, 1, tmin, tmax, steps);

    for (i = 0; i < steps; i++) {
	char outfile[64];

#if 0
	bout = img2bmp(m[i]);
    
	if (bout == NULL) {
	    printf("Error in img2bmp conversion\n");
	    return 1;
	}
#endif

	sprintf(outfile, "%s%03d.bmp", argv[7], i);

	f = fopen(outfile, "w");
	if (f == NULL) {
	    printf("Error opening %s for writing\n", outfile);
	    return 1;
	}
	    
	// write_bmp(f, bout);
	// fclose(f);
        img_write_bmp_file(m[i], outfile);

	// free_bmp(bout);
    }
#else
    iout = img_morph(i1, i2);
    bout = img2bmp(iout);
    
    if (bout == NULL) {
	printf("Error in img2bmp conversion\n");
	return 1;
    }

    f = fopen(argv[7], "w");
    if (f == NULL) {
	printf("Error opening %s for writing\n", argv[7]);
	return 1;
    }

    write_bmp(f, bout);
    fclose(f);
    
    free_bmp(bout);
#endif

    // free_bmp(b1);
    // free_bmp(b2);
    img_free(i1);
    img_free(i2);

    return 0;
}
