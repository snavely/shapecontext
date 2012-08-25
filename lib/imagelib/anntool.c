/* ann.c */
/* Contains routines that use the ANN interface defined in 
 * anniface.h */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "anniface.h"
#include "binvolume.h"
#include "bmp.h"
#include "circle.h"
#include "histogram.h"
#include "image.h"
#include "resample.h"
#include "vector.h"

static double ann_gradient_weight = 3.0;

#if 0
/* Create a kd tree of the points stored in a binary volume */
ANNkd_tree_t *bvol_create_kd_tree(bvol_t *b) {
    ANNkd_tree_t *tree;
    int n = bvol_num_points(b);
    v3_t *pts = malloc(sizeof(v3_t) * n);
    int idx = 0;

    if (b->type == BVOL_BITS) {
        int x, y, z;
        for (z = 0; z < b->d; z++)
            for (y = 0; y < b->h; y++)
                for (x = 0; x < b->w; x++)
                    if (bvol_getbit(b, x, y, z))
                        pts[idx++] = v3_new((double) x, 
                                            (double) y, 
                                            (double) z);

    } else { /* b->type == BVOL_HEIGHTS */
        int x, y;
        
        for (y = 0; y < b->h; y++) {
            for (x = 0; x < b->w; x++) {
                int z = BVOL_HEIGHT(b, x, y);
                pts[idx++] = v3_new((double) x,
                                    (double) y,
                                    (double) z);
            }
        }
    }

    tree = create_ann_tree_3D(n, pts, v3_new(1.0, 1.0, get_ann_z_weight()));
    free(pts);

    return tree;
}
#endif

/* Create a 3D kd-tree from a (grayscale) image */
ANNkd_tree_t *img_create_kd_tree_grayscale(img_t *img) {
    int w = img->w, h = img->h, n = w * h;
    int x, y, idx = 0;
    v3_t *pts = malloc(sizeof(v3_t) * n);
    ANNkd_tree_t *tree;
    clock_t start, end;
    
    for (y = 0; y < h; y++) {
        for (x = 0; x < w; x++) {
            // int z = IMG_GET_PIXEL(img, x, y).r;
	    color_t c = img_get_pixel(img, x, y);
	    double intensity = (double) (c.r + c.g + c.b) / 3.0;
            pts[idx++] = v3_new((double) x, (double) y, intensity);
        }
    }

    start = clock();
    tree = create_ann_tree_3D(n, pts, v3_new(1.0, 1.0, get_ann_z_weight()));
    end = clock();

    // printf("kd-tree creation took %0.3fs\n", ((double) (end - start)) / CLOCKS_PER_SEC);

    free(pts);

    return tree;
}

/* Create a 3D kd-tree from a (grayscale) image, filling in some
 * interior points using linear interpolation */
ANNkd_tree_t *img_create_kd_tree_grayscale_fill(img_t *img, int fill_ratio) {
    int w = img->w, h = img->h, n = (fill_ratio * (w - 1) + 1) * (fill_ratio * (h - 1) + 1);
    double x, y;
    int idx = 0;
    v3_t *pts = malloc(sizeof(v3_t) * n);
    ANNkd_tree_t *tree;
    double d = 1.0 / (double) fill_ratio;
    clock_t start, end;

    for (y = 0.0; y <= h - 1; y += d) {
        for (x = 0.0; x <= w - 1; x += d) {
            fcolor_t c = pixel_lerp(img, x, y);
	    double intensity = (double) (c.r + c.g + c.b) / 3.0;
            pts[idx++] = v3_new(x, y, intensity);
        }
    }

    start = clock();
    tree = create_ann_tree_3D(n, pts, v3_new(1.0, 1.0, get_ann_z_weight()));
    end = clock();
    
    // printf("kd-tree creation took %0.3fs\n", ((double) (end - start)) / CLOCKS_PER_SEC);

    free(pts);

    return tree;
}

/* Create a 4D kd-tree from a grayscale image.  The fourth dimension
 * is the gradient */
ANNkd_tree_t *img_create_kd_tree_gradient(img_t *img) {
    int w = img->w, h = img->h, n = w * h;
    ANNkd_tree_t *tree;
    vec_t *pts = malloc(sizeof(vec_t) * n), axis_weights;
    int x, y, idx = 0, i;
    
    for (y = 0; y < h; y++) {
	for (x = 0; x < w; x++) {
	    color_t col = img_get_pixel(img, x, y);
	    pts[idx] = vec_new(4);

	    Vn(pts[idx], 0) = (double) x;
	    Vn(pts[idx], 1) = (double) y;
	    Vn(pts[idx], 2) = (double) col.r;
	    Vn(pts[idx], 3) = img_gradient(img, x, y);

	    idx++;
	}
    }

    axis_weights = vec_new(4);

    Vn(axis_weights, 0) = 1.0;
    Vn(axis_weights, 1) = 1.0;
    Vn(axis_weights, 2) = get_ann_z_weight();
    Vn(axis_weights, 3) = ann_gradient_weight;

    tree = create_ann_tree(n, 4, pts, axis_weights);

    /* Cleanup */
    for (i = 0; i < n; i++)
        vec_free(pts[i]);

    vec_free(axis_weights);
    free(pts);

    return tree;
}

/* Create a 2+n^2 dimensional kd-tree from a grayscale image, 
 * with dimensions for each pixel value in the neighborhood */
ANNkd_tree_t *img_create_kd_tree_nhood(img_t *img, int n) {
    int w = img->w, h = img->h, num_pts = (w - n + 1) * (h - n + 1);
    int d = 2 + n * n; /* 2 spatial dimensions, n*n color dimensions */

    vec_t *pts = malloc(sizeof(vec_t) * num_pts), axis_weights;
    ANNkd_tree_t *tree;
    int x, y, idx = 0, i, j;

    clock_t start, end;

    for (y = n / 2; y <= h - 1 - n / 2; y++) {
	for (x = n / 2; x <= w - 1 - n / 2; x++) {
	    pts[idx] = vec_new(d);

	    /* First two dimensions are spatial (lower left corner of
	     * neighborhood */
	    Vn(pts[idx], 0) = (double) x;
	    Vn(pts[idx], 1) = (double) y;
	    
	    /* Remaining dimensions are colors */
	    for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
		    color_t c = img_get_pixel(img, x + j - n / 2, y + i - n / 2);
		    Vn(pts[idx], 2 + i * n + j) = (double) (c.r + c.g + c.b) / 3.0;
		}
	    }
	    
	    idx++;
	}
    }

    /* Fill the axis weights */
    axis_weights = vec_new(d);

    Vn(axis_weights, 0) = 1.0;
    Vn(axis_weights, 1) = 1.0;

    for (i = 2; i < d; i++)
	Vn(axis_weights, i) = get_ann_z_weight();

    /* Create the tree! */
    start = clock();
    tree = create_ann_tree(num_pts, d, pts, axis_weights);
    end = clock();

    // printf("kd-tree creation took %0.3fs\n", ((double) (end - start)) / CLOCKS_PER_SEC);

    /* Cleanup */
    for (i = 0; i < num_pts; i++)
	vec_free(pts[i]);

    vec_free(axis_weights);
    free(pts);

    return tree;
}

/* Create a 2+3*n^2 dimensional kd-tree from a rgb image, 
 * with dimensions for each pixel value in the neighborhood */
ANNkd_tree_t *img_create_kd_tree_nhood_rgb(img_t *img, int n, int fill_ratio, 
					   int xmin, int xmax, int ymin, int ymax) {
    int w = xmax - xmin + 1, h = ymax - ymin + 1, 
	num_pts = (fill_ratio * (w - n) + 1) * (fill_ratio * (h - n) + 1);
    int d = 2 + 3 * n * n; /* 2 spatial dimensions, 3*n*n color dimensions */
    double dx = 1.0 / (double) fill_ratio;
    double dy = dx;
    vec_t *pts, axis_weights;
    ANNkd_tree_t *tree;
    int i;

    /* The computation above is not accurate when the image is too
     * small */
    if (w < n || h < n)
	num_pts = 0;

    if (num_pts == 0) {

	/* The image is too small -- create a tree with a single point */
	num_pts = 1;

	pts = malloc(sizeof(vec_t));
	pts[0] = vec_new(d);

	for (i = 0; i < d; i++)
	    Vn(pts[0], i) = 0.0;

	printf("[img_create_kd_tree_nhood_rgb] Error: no points\n");
    } else {
	int idx = 0, j;
	double x, y;
	int rad = n / 2;

	pts = malloc(sizeof(vec_t) * num_pts);

	for (y = ymin + rad; y <= ymin + h - 1 - rad; y += dy) {
	    for (x = xmin + rad; x <= xmin + w - 1 - rad; x += dx) {
		pts[idx] = vec_new(d);
		
		/* First two dimensions are spatial (lower left corner of
		 * neighborhood */
		Vn(pts[idx], 0) = x;
		Vn(pts[idx], 1) = y;
	    
		/* Remaining dimensions are colors */
		for (i = 0; i < n; i++) {
		    for (j = 0; j < n; j++) {
			/* color_t c = img_get_pixel(img, x + j - n / 2, y + i - n / 2); */
			fcolor_t c = pixel_lerp(img, x + j - rad, y + i - rad);
			Vn(pts[idx], 2 + (i * n + j) * 3 + 0) = (double) c.r;
			Vn(pts[idx], 2 + (i * n + j) * 3 + 1) = (double) c.g;
			Vn(pts[idx], 2 + (i * n + j) * 3 + 2) = (double) c.b;
		    }
		}
	    
		idx++;
	    }
	}
    }

    /* Fill the axis weights */
    axis_weights = vec_new(d);

    Vn(axis_weights, 0) = 1.0;
    Vn(axis_weights, 1) = 1.0;

    for (i = 2; i < d; i++)
	Vn(axis_weights, i) = get_ann_z_weight();
    
    /* Create the tree! */
    tree = create_ann_tree(num_pts, d, pts, axis_weights);

    /* Cleanup */
    for (i = 0; i < num_pts; i++)
	vec_free(pts[i]);

    vec_free(axis_weights);
    free(pts);

    return tree;
}

/* Create a 5D kd-tree from a (rgb) image, with separate dimensions
 * for each of the tree color channels */
ANNkd_tree_t *img_create_kd_tree_rgb(img_t *img, int fill_ratio) {
    int w = img->w, h = img->h, n = (fill_ratio * (w - 1) + 1) * (fill_ratio * (h - 1) + 1);
    int d = 5;
    double dx = 1.0 / (double) fill_ratio;
    double dy = dx;

    int idx = 0, i;
    double x, y;
    vec_t *pts = malloc(sizeof(vec_t) * n), axis_weights;
    ANNkd_tree_t *tree;
    
    for (y = 0.0; y <= h - 1; y += dy) {
        for (x = 0.0; x <= w - 1; x += dx) {
            fcolor_t col = pixel_lerp(img, x, y);
            pts[idx] = vec_new(d);

            Vn(pts[idx], 0) = (double) x;
            Vn(pts[idx], 1) = (double) y;
            Vn(pts[idx], 2) = (double) col.r;
            Vn(pts[idx], 3) = (double) col.g;
            Vn(pts[idx], 4) = (double) col.b;
            
            idx++;
        }
    }

    axis_weights = vec_new(d);

    Vn(axis_weights, 0) = 1.0;
    Vn(axis_weights, 1) = 1.0;
    Vn(axis_weights, 2) = get_ann_z_weight();
    Vn(axis_weights, 3) = get_ann_z_weight();
    Vn(axis_weights, 4) = get_ann_z_weight();

    tree = create_ann_tree(n, d, pts, axis_weights);

    /* Cleanup */
    for (i = 0; i < n; i++)
        vec_free(pts[i]);

    vec_free(axis_weights);
    free(pts);

    return tree;
}

/* Return the approximate nearest neighbor to q out of the points in
 * the (weighted) tree */
v3_t query_weighted_ann_tree_3D_pt(ANNkd_tree_t *tree, v3_t q, double eps) {
    v3_t pt;

    Vz(q) *= get_ann_z_weight();
    pt = query_ann_tree_3D_pt(tree, q, eps);
    Vz(pt) /= get_ann_z_weight();
    return pt;
}
