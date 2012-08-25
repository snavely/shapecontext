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

void readPoints(char *fname1, char *fname2, img_dmap_t **out1to2, img_dmap_t **out2to1) 
{

    /* Read the first file */
    FILE *f = fopen(fname1, "r");

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


    /* Read the second file */
    f = fopen(fname2, "r");

    /* Read boundaries */
    int xMin2, yMin2, xMax2, yMax2;
    fscanf(f, "%d %d %d %d", &xMin2, &yMin2, &xMax2, &yMax2);
    xMax2 = xMin2 + xMax2;
    yMax2 = yMin2 + yMax2;

    // 139 379 304 354

    /* Read the feature positions */
    std::vector<point_t> pts2;
    for (int j = 0; j < NUM_FEATURES; j++) {
	int x, y;
	point_t pt;
	fscanf(f, "%d %d", &x, &y);

	pt.x = x;
	pt.y = y;
	    
	pts2.push_back(pt);
    }

    fclose(f);


    int w1 = xMax1 - xMin1 + 1;
    int h1 = yMax1 - yMin1 + 1;

    int w2 = xMax2 - xMin2 + 1;
    int h2 = yMax2 - yMin2 + 1;

    printf("(%d, %d) (%d, %d)\n", xMin1, yMin1, xMax1, yMax1);
    printf("(%d, %d) (%d, %d)\n", xMin2, yMin2, xMax2, yMax2);
    printf("(%d x %d), (%d x %d)\n", w1, h1, w2, h2);

    /* Select the image with the larger dimensions */
    int w = MAX(w1, w2);
    int h = MAX(h1, h2);

    double wratio1 = ((double) w) / ((double) w1);
    double wratio2 = ((double) w) / ((double) w2);
    double hratio1 = ((double) h) / ((double) h1);
    double hratio2 = ((double) h) / ((double) h2);

    img_dmap_t *d1to2;
    img_dmap_t *d2to1;

    d1to2 = img_dmap_new(w, h);
    d2to1 = img_dmap_new(w, h);

    /* Fill in dmaps */
    for (int y = 0; y < h; y++) {
	for (int x = 0; x < w; x++) {
	    d1to2->dists[y * w + x] = DBL_MAX;
	}
    }

    for (int y = 0; y < h; y++) {
	for (int x = 0; x < w; x++) {
	    d2to1->dists[y * w + x] = DBL_MAX;
	}
    }

    for (int i = 0; i < NUM_FEATURES; i++) {
	if (i == 0 || i == 12)
	    continue;

	int x1 = iround((pts1[i].x - xMin1) * wratio1);
	int y1 = iround((pts1[i].y - yMin1) * hratio1);
	int x2 = iround((pts2[i].x - xMin2) * wratio2);
	int y2 = iround((pts2[i].y - yMin2) * hratio2);
	
	d1to2->dists[y1 * w + x1] = 0.0;
	d1to2->nns[y1 * w + x1] = v2_new(x2, y2);

	d2to1->dists[y2 * w + x2] = 0.0;
	d2to1->nns[y2 * w + x2] = v2_new(x1, y1);
    }

    *out1to2 = d1to2;
    *out2to1 = d2to1;

    printf("convert -crop %dx%d+%d+%d x.tga tmp.tga\n",
	   w1, h1, xMin1, IMAGE_HEIGHT - yMax1 - 1);
    printf("convert -geometry %dx%d! tmp.tga out.tga \n",
	   w, h);

    printf("convert -crop %dx%d+%d+%d x.tga tmp.tga\n",
	   w2, h2, xMin2, IMAGE_HEIGHT - yMax2 - 1);
    printf("convert -geometry %dx%d! tmp.tga out.tga \n",
	   w, h);
}


void printUsage(char *name) 
{
    printf("Usage: %s <pts1.in> <pts2.in> <dmap1.out> <dmap2.out>\n", name);
}

int main(int argc, char **argv) 
{
    char *in_pts1, *in_pts2;
    char *out_dmap1, *out_dmap2;

    img_dmap_t *dmap1, *dmap2;
    img_dist_pyr_t *dpyr1, *dpyr2;

    if (argc < 5) {
	printUsage(argv[0]);
	return -1;
    }

    in_pts1 = argv[1];
    in_pts2 = argv[2];
    out_dmap1 = argv[3];
    out_dmap2 = argv[4];

    readPoints(in_pts1, in_pts2, &dmap1, &dmap2);

    dpyr1 = dmap2dpyr(dmap1);
    dpyr2 = dmap2dpyr(dmap2);
    img_write_distance_pyramid_file(dpyr1, out_dmap1);
    img_write_distance_pyramid_file(dpyr2, out_dmap2);

    return 0;
}
