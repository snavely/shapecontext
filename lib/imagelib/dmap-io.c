/* dmap-io.c */
/* Routines for reading / writing / rendering distance maps */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>

#include "defines.h"
#include "dmap.h"
#include "dmap-io.h"
#include "fileio.h"
#include "image.h"
#include "vector.h"

/* Create a clean dmap */
img_dmap_t *img_dmap_new(int w, int h) {
    img_dmap_t *dmap = (img_dmap_t *)malloc(sizeof(img_dmap_t));

    dmap->w = w;
    dmap->h = h;
    dmap->dists = (double *)malloc(sizeof(double) * w * h);
    dmap->nns = (v2_t *)malloc(sizeof(v2_t) * w * h);
    dmap->uppers = NULL;

    return dmap;
}

/* Produces an image of the chamfer metric between two images */
img_t *img_dmap_render(img_dmap_t *dmap) {
    int w = dmap->w, h = dmap->h;
    double *dists = dmap->dists;

    img_t *img_out = img_new(w, h);
    int idx, num_pixels = w * h;

    double max_dist = 0.0;
    double min_dist = DBL_MAX;

    for (idx = 0; idx < num_pixels; idx++) {
	if (dists[idx] == DBL_MAX)
	    continue;

	if (dists[idx] > max_dist)
	    max_dist = dists[idx];

	if (dists[idx] < min_dist)
	    min_dist = dists[idx];
    }

    /* Fill in the pixels */
    for (idx = 0; idx < num_pixels; idx++) {
	if (dists[idx] == DBL_MAX) {
	    img_out->pixels[idx].r = img_out->pixels[idx].g = 0;
	    img_out->pixels[idx].b = 255;
	} else {
	    int c = (int) (256 * (dists[idx] - min_dist) / (max_dist - min_dist));
	    img_out->pixels[idx].r = img_out->pixels[idx].g = img_out->pixels[idx].b = c;
	}
    }

    return img_out;
}

/* Produces an image of the disparity between two images */
img_t *img_dmap_render_disparity(img_dmap_t *dmap) {
    int w = dmap->w, h = dmap->h;
    double *dists = dmap->dists;
    int x, y;

    img_t *img_out = img_new(w, h);
    int idx; // , num_pixels = w * h;

    double max_dist = 0.0;
    double min_dist = DBL_MAX;

    idx = 0;
    for (y = 0; y < dmap->h; y++) {
	for (x = 0; x < dmap->w; x++) {
	    double dist;
	    double dx, dy;
	    v2_t nn;

	    if (dists[idx] == DBL_MAX) {
		idx++;
		continue;
	    }
	    
	    nn = dmap->nns[idx];

	    dx = Vx(nn) - x;
	    dy = Vy(nn) - y;

	    // dist = sqrt(dx * dx + dy * dy);
	    dist = dx;

	    if (dist < min_dist)
		min_dist = dist;
	    
	    if (dist > max_dist)
		max_dist = dist;

	    idx++;
	}
    }

    printf("mindist = %0.3f\n", min_dist);
    printf("maxdist = %0.3f\n", max_dist);

    /* Fill in the pixels */
    idx = 0;
    for (y = 0; y < dmap->h; y++) {
	for (x = 0; x < dmap->w; x++) {
	    if (dists[idx] == DBL_MAX) {
		img_out->pixels[idx].r = img_out->pixels[idx].g = 0;
		img_out->pixels[idx].b = 255;
	    } else {
		v2_t nn = dmap->nns[idx];
		double dx = Vx(nn) - x;
		double dy = Vy(nn) - y;
		// double dist = sqrt(dx * dx + dy * dy);
		double dist = dx;
		// int c = (int) (256 * sqrt(dist / max_dist));
		int c = CLAMP((int) rint(1 * (-dist + 15)), 0, 255);
		img_out->pixels[idx].r = img_out->pixels[idx].g = img_out->pixels[idx].b = c;
	    }

	    idx++;
	}
    }

    return img_out;
}

/* Produces an image of the flow map */
img_t *img_dmap_render_flow(img_dmap_t *dmap) {
    
    int w = dmap->w, h = dmap->h;
    double *dists = dmap->dists;
    v2_t *nns = dmap->nns;

    img_t *img_out = img_new(w, h);
    int idx, x, y;

    double max_dist = 0.0;

    /* Fill in the pixels */
    for (y = 0; y < h; y++) {
	for (x = 0; x < w; x++) {
	    idx = y * dmap->w + x;
	    v2_t nn = nns[idx];

#define FACTOR 6.0
	    double dx = CLAMP(FACTOR * (Vx(nn) - x) + 127.0, 0.0, 255.0);
	    double dy = CLAMP(FACTOR * (Vy(nn) - y) + 127.0, 0.0, 255.0);
#undef FACTOR

	    if (dists[idx] == DBL_MAX) {
		img_out->pixels[idx].r = img_out->pixels[idx].g = 0;
		img_out->pixels[idx].b = 0;
	    } else {
		img_out->pixels[idx].r = (int) rint(dx);
		img_out->pixels[idx].g = (int) rint(dy);
		img_out->pixels[idx].b = 0.0;
	    }
	}
    }

    return img_out;    
}

void img_dmap_write(FILE *f, img_dmap_t *dmap) {
    int count;
    short int has_uppers;

    /* Write the identifier */
    write_word((u_int32_t *)"DMAP", f);

    /* Write the width and height */
    write_short(&dmap->w, f);
    write_short(&dmap->h, f);

    if (dmap->uppers == NULL)
	has_uppers = 0;
    else
	has_uppers = 1;
    
    write_short(&has_uppers, f);

    /* Write the distances */
    for (count = 0; count < dmap->w * dmap->h; count++)
	write_double(&(dmap->dists[count]), f);

    /* Write each vector */
    for (count = 0; count < dmap->w * dmap->h; count++) {
	write_double(&(Vx(dmap->nns[count])), f);
	write_double(&(Vy(dmap->nns[count])), f);
    }

    if (dmap->uppers != NULL) {
	for (count = 0; count < dmap->w * dmap->h; count++) {
	    write_short(&Vx(dmap->uppers[count]), f);
	    write_short(&Vy(dmap->uppers[count]), f);
	}
    }
}

img_dmap_t *img_dmap_read_file(char *fname) {    
    FILE *f = open_file(fname, "r");
    img_dmap_t *dmap = img_dmap_read(f);

    if (dmap == NULL) {
	printf("[img_dmap_read_file] Error reading dmap file %s.\n", fname);
    }

    fclose(f);

    return dmap;
}

img_dmap_t *img_dmap_read(FILE *f) {
    int count;
    u_int16_t w, h;
    char id[5];
    img_dmap_t *dmap;
    short int has_uppers;

    /* Read the identifier */
    read_word((u_int32_t *)id, f);
    id[4] = 0;

    if (strcmp(id, "DMAP") != 0) {
	printf("[img_read_distance_map] Invalid distance map file\n");
	return NULL;
    }

    /* Read the width and height */
    read_short(&w, f);
    read_short(&h, f);
    read_short(&has_uppers, f);

    /* Initialize the map */
    dmap = img_dmap_new(w, h);

    /* Read the distances */
    for (count = 0; count < w * h; count++) {
	read_double(&(dmap->dists[count]), f);
    }
    
    /* Read the nearest neighbors */
    for (count = 0; count < w * h; count++) {
	read_double(&Vx(dmap->nns[count]), f);
	read_double(&Vy(dmap->nns[count]), f);
    }

    if (has_uppers) {
	dmap->uppers = (iv2_t *) malloc(w * h * sizeof(iv2_t));

	/* Read the uppers */
	for (count = 0; count < w * h; count++) {
	    read_short(&(Vx(dmap->uppers[count])), f);
	    read_short(&(Vy(dmap->uppers[count])), f);
	}
    }

    return dmap;
}

void img_dmap_free(img_dmap_t *dmap) 
{
    free(dmap->dists);
    free(dmap->nns);
    free(dmap);
}
