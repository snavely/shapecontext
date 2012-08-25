/* pts2dmap.c */
/* Take a set of input points and generate a dmap */

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "defines.h"
#include "dmap.h"
#include "dmap-io.h"
#include "image.h"
#include "pyramid.h"
#include "pyramid-io.h"
#include "util.h"

#include "point.h"

#define NUM_FEATURES 13

#define IMAGE_WIDTH 512
#define IMAGE_HEIGHT 768

img_t *readPoints(char *in_pts, char *in_img) 
{

    /* Read the first file */
    FILE *f = fopen(in_pts, "r");

    /* Read boundaries */
    int xMin1, yMin1, xMax1, yMax1;
    fscanf(f, "%d %d %d %d", &xMin1, &yMin1, &xMax1, &yMax1);
    xMax1 = xMin1 + xMax1;
    yMax1 = yMin1 + yMax1;

    /* Read the feature positions */
    std::vector<point_t> pts1;
    for (int j = 0; j < NUM_FEATURES; j++) {
	int x, y;
	point_t pt;
	fscanf(f, "%d %d", &x, &y);

	pt.x = x;
	pt.y = y;
	    
	pts1.push_back(pt);
    }

    fclose(f);

    img_t *img = img_read_bmp_file(in_img);

    /* Draw images with selected points */
    img_t *img_out = img_copy(img);

    for (int i = 0; i < NUM_FEATURES; i++) {
	if (i == 0 || i == 12)
	    continue;

	point_t pt1 = pts1[i];
	int x = iround(pt1.x);
	int y = iround(pt1.y);

#define RAD 6
#if 0
	img_draw_line(img_out, x - RAD, y - RAD, x - RAD, y + RAD, 
		      0xff, 0xff, 0xff);
	img_draw_line(img_out, x - RAD, y + RAD, x + RAD, y + RAD, 
		      0xff, 0xff, 0xff);
	img_draw_line(img_out, x + RAD, y + RAD, x + RAD, y - RAD, 
		      0xff, 0xff, 0xff);
	img_draw_line(img_out, x + RAD, y - RAD, x - RAD, y - RAD, 
		      0xff, 0xff, 0xff);
#endif

	img_draw_pt(img_out, x, y, 2 * RAD + 1, 0xff, 0x0, 0x0);
    }

    return img_out;
}


void printUsage(char *name) 
{
    printf("Usage: %s <img.in> <pts.in> <img.out>\n", name);
}

int main(int argc, char **argv) 
{
    char *in_pts;
    char *in_img;
    char *out_img;

    if (argc < 4) {
	printUsage(argv[0]);
	return -1;
    }

    in_img = argv[1];
    in_pts = argv[2];
    out_img = argv[3];

    img_t *img_out = readPoints(in_pts, in_img);
    img_write_bmp_file(img_out, out_img);

    return 0;
}
