/* dmap.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "dmap.h"
#include "dmap-io.h"
#include "error.h"
#include "image.h"

/* Create a map of distances between b and a */
img_dmap_t *img_dmap_create(img_t *a, img_t *b, error_metric_t em, double eps, int full, int use_sqrt) {
    int rad;
    ANNkd_tree_t *atree = img_create_ann_tree(a, em, &rad);
    img_dmap_t *dmap = img_dmap_new(b->w, b->h);

    int x, y;

    for (y = 0; y < b->h; y++) {
	for (x = 0; x < b->w; x++) {
	    if (img_pixel_is_valid(b, x, y)) {
		double dist;
		vec_t sample = em_create_query_pt(b, em, x, y);
		vec_t nn = query_ann_tree_pt(atree, sample, eps, &dist, full, use_sqrt);

		dmap->nns[y * b->w + x] = v2_new(Vx(nn), Vy(nn));
		dmap->dists[y * b->w + x] = dist;

		vec_free(sample);
		vec_free(nn);
	    } else {
		dmap->dists[y * b->w + x] = DBL_MAX;
		dmap->nns[y * b->w + x] = v2_new(0.0, 0.0);
	    }
	}
    }

    free_ann_tree(atree);

    return dmap;
}

/* Allocate uppers for a dmap */
void img_dmap_allocate_uppers(img_dmap_t *map) {
    map->uppers = (iv2_t *)malloc(sizeof(iv2_t) * map->w * map->h);
}

/* Return the sub map indicated by the given bounds */
img_dmap_t *img_dmap_sub_map(img_dmap_t *map, int xmin, int ymin, int w, int h) {
    img_dmap_t *map_out = img_dmap_new(w, h);
    
    int xmax = xmin + w - 1;
    int ymax = ymin + h - 1;

    int x, y;

    for (y = ymin; y <= ymax; y++) {
	int y_off = y - ymin;
	for (x = xmin; x <= xmax; x++) {
	    int x_off = x - xmin;
	    int idx = y * map->w + x;
	    int idx_off = y_off * w + x_off;

	    map_out->dists[idx_off] = map->dists[idx];
	    map_out->nns[idx_off] = map->nns[idx];
	}
    }

    return map_out;
}

/* Shift the nearest neighbors in a dmap by the amoung specified */
img_dmap_t *img_dmap_shift(img_dmap_t *map, int x_shift, int y_shift) {
    img_dmap_t *map_out = img_dmap_new(map->w, map->h);
    int x, y;

    for (y = 0; y < map->h; y++) {
	for (x = 0; x < map->w; x++) {
	    int idx = y * map->w + x;
	    v2_t nn;

	    map_out->dists[idx] = map->dists[idx];
	    nn = map->nns[idx];
	    Vx(nn) += x_shift;
	    Vy(nn) += y_shift;
	    map_out->nns[idx] = nn;
	}
    }

    return map_out;    
}

/* Shrink wrap the given dmap, making it as small as possible */
void img_dmap_shrink_wrap(img_dmap_t *map, img_dmap_t *map_out) {
    /* Find the smallest enclosing width and height */
    int x_max = 0, y_max = 0;
    int x, y, idx;
    int w_new, h_new;

    idx = 0;
    for (y = 0; y < map->h; y++) {
	for (x = 0; x < map->w; x++, idx++) {
	    if (map->dists[idx] != DBL_MAX) {
		if (x > x_max)
		    x_max = x;
		if (y > y_max)
		    y_max = y;
	    }
	}
    }

    w_new = x_max + 1;
    h_new = y_max + 1;

    map_out->w = w_new;
    map_out->h = h_new;
    map_out->dists = (double *)malloc(sizeof(double) * w_new * h_new);
    map_out->nns = (v2_t *)malloc(sizeof(v2_t) * w_new * h_new);
    map_out->uppers = NULL;

    if (map->uppers != NULL) {
	map_out->uppers = (iv2_t *) malloc(sizeof(iv2_t) * w_new * h_new);
    }

    /* Copy over the map */
    for (y = 0; y < h_new; y++) {
	memcpy(map_out->dists + y * w_new, map->dists + y * map->w, sizeof(double) * w_new);
	memcpy(map_out->nns + y * w_new, map->nns + y * map->w, sizeof(v2_t) * w_new);

	if (map->uppers != NULL) {
	    memcpy(map_out->uppers + y * w_new, map->uppers + y * map->w, sizeof(iv2_t) * w_new);
	}
    }
}
