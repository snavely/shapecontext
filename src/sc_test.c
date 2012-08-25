/* sc_test.c */
/* Test driver for shape context app */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "dmap.h"
#include "dmap-io.h"
#include "image.h"
#include "pyramid.h"
#include "pyramid-io.h"
#include "util.h"
#include "vector.h"

#include "shapecontext.h"

int main(int argc, char **argv) 
{
    char *in_image = NULL;
    img_t *img_edge;
    int x, y, r, w;
    char *out_image;
    img_t *img_out;
    double *descriptor;

    if (argc != 5) {
	printf("Usage: %s <edge-image.bmp> <x> <y> <out.bmp>\n", argv[0]);
	return -1;
    }

    in_image = argv[1];
    x = atoi(argv[2]);
    y = atoi(argv[3]);
    out_image = argv[4];

    img_edge = img_read_bmp_file(in_image);
    
#define NUM_RINGS 5 // 4 // 6 // 14
#define NUM_WEDGES 12
#define FACTOR 2.4 // 3.0 // 1.4
#define SIGMA 1.0 /* For now, sigma is not used */
#define NORMALIZE 1

    /* Compute shape context for the given point */
    descriptor = (double *) malloc(sizeof(double) * NUM_RINGS * NUM_WEDGES);

    compute_shape_context_pt(img_edge, x, y, 
			     NUM_RINGS, NUM_WEDGES, FACTOR, SIGMA, NORMALIZE, 
			     descriptor);

#define CELL_SIZE 32
    img_out = img_new(NUM_WEDGES * CELL_SIZE, NUM_RINGS * CELL_SIZE);

    for (r = 0; r < NUM_RINGS; r++) {
	for (w = 0; w < NUM_WEDGES; w++) {
	    double c = descriptor[r * NUM_WEDGES + w];
	    int ci = iround(c);
	    for (y = r * CELL_SIZE; y < (r+1) * CELL_SIZE; y++) {
		for (x = w * CELL_SIZE; x < (w+1) * CELL_SIZE; x++) {
		    img_set_pixel(img_out, x, y, ci, ci, ci);
		}
	    }
	}
    }

    img_write_bmp_file(img_out, out_image);

    return 0;
}
