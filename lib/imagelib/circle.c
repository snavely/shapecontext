/* circle.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "circle.h"

/* The antiderivative of sqrt(r * r - x * x) */
static double integrate(double r, double x) {
    double subexpr = x * sqrt(r * r - x * x);
    return 0.5 * subexpr - 0.5 * r * r * atan(subexpr / (- r * r + x * x));
}

/* Compute the area of the intersection of a circle with a unit square
 * (pixel).
 * r:     the radius of the circle
 * xc,yc: the coordinates of the center of the circle
 * x,y:   the coordinates of the lower left corner of the pixel */
double circle_intersect_pixel_old(int r, int xc, int yc, int x, int y) {
    int xmin, ymin, xmax, ymax, tmp;

    x -= xc;
    y -= yc;

    /* Check if the pixel is completely outside or completely inside
     * the circle by finding the corners of the pixel with least and
     * greatest distance from the center */
    
    xmin = (x >= 0) ? x : x + 1;
    ymin = (y >= 0) ? y : y + 1;

    xmax = (x >= 0) ? (x + 1) : x;
    ymax = (y >= 0) ? (y + 1) : y;

    if (xmin * xmin + ymin * ymin >= r * r)
	return 0.0;
    if (xmax * xmax + ymax * ymax <= r * r)
	return 1.0;
    
    /* We know that the circle intersects the pixel, but not
     * completely.  Now find the pixel in the first quadrant with the
     * same area as (x, y) */

    if (x < 0) 
	x = -(x+1);
    if (y < 0)
	y = -(y+1);
    
    if (x != y) {
	/* Normal case.  Find the pixel in the second octant with the
	 * same area as (x, y) */
	
	if (x > y) {
	    tmp = x;
	    x = y;
	    y = tmp;
	}

	/* Compute the area under the curve */
	return integrate((double) r, x + 1) - integrate((double) r, x) - y;
    } else {
	/* Special case.  We can't take the entire integral, we must
	 * find a bound on it */
	double bound = sqrt(r * r - y * y);

	return integrate((double) r, bound) - integrate((double) r, x) - y * (bound - x);
    }
}

/* Compute the area of the intersection of a circle with a unit square
 * (pixel).
 * r:     the radius of the circle
 * xc,yc: the coordinates of the center of the circle
 * x,y:   the coordinates of the center of the pixel */
double circle_intersect_pixel(double r, int xc, int yc, int xp, int yp) {
    double xmin, ymin, xmax, ymax, tmp;
    double x = (double) xp - 0.5;
    double y = (double) yp - 0.5;
    double ubound, lbound;

    x -= xc;
    y -= yc;


    /* If the radius is less than half a pixel width, we consider this
     * a special case */
    if (r <= 0.5) {
	if (x == -0.5 && y == -0.5) {
	    return M_PI * r * r;
	} else {
	    return 0.0;
	}
    }

    /* Check if the pixel is completely outside or completely inside
     * the circle by finding the corners of the pixel with least and
     * greatest distance from the center */

    if (x == -0.5) {
	xmin = 0.0; /* Special case when the pixel straddles the y-axis */
    } else {
	xmin = (x >= 0) ? x : (x + 1);
    }
    
    if (y == -0.5) {
	ymin = 0.0; /* Special case when the pixel straddles the x-axis */
    } else {
	ymin = (y >= 0) ? y : (y + 1);
    }
    
    xmax = (x >= 0) ? (x + 1) : x;
    ymax = (y >= 0) ? (y + 1) : y;

    if (xmin * xmin + ymin * ymin >= r * r)
	return 0.0;
    if (xmax * xmax + ymax * ymax <= r * r)
	return 1.0;
    
    /* We know that the circle intersects the pixel, but not
     * completely.  Now find the pixel in the first quadrant with the
     * same area as (x, y) */

    if (x < -0.5) 
	x = -(x+1);
    if (y < -0.5)
	y = -(y+1);

    /* Find the pixel in the second octant with the same area as (x, y) */
    if (x > y) {
	tmp = x;
	x = y;
	y = tmp;
    }

    ubound = sqrt(r * r - y * y);

    if ((y + 1) * (y + 1) > r * r) {
	lbound = x;
    } else {
	lbound = sqrt(r * r - (y + 1) * (y + 1));
    }
    
#if 1
    /* Four cases */
    if (lbound <= x && ubound >= x + 1) {
	/* Easy case... just take the integral */
	return integrate(r, x + 1) - integrate(r, x) - y;
    } else if (lbound <= x && ubound < x + 1) {
	return integrate(r, ubound) - integrate(r, x) - y * (ubound - x);
    } else if (lbound > x && ubound >= x + 1) {
	return integrate(r, x + 1) - integrate(r, lbound) - y * (x + 1 - lbound) + (lbound - x);
    } else if (lbound > x && ubound < x + 1) {
	return integrate(r, ubound) - integrate(r, lbound) - y * (ubound - lbound) + (lbound - x);
    } else {
	printf("error in circle()!\n");
    }
    
    return 0.0;

#else
    if (x != y && (ubound - x) >= 1.0) {
	/* Compute the area under the curve */
	return integrate((double) r, x + 1) - integrate((double) r, x) - y;
    } else {
	/* Special case.  We can't take the entire integral, we must
	 * find a bound on it */

	if (ubound - x < 1.0) {
	    return integrate(r, ubound) - integrate(r, x) - y * (ubound - x);
	} else { /* x == y */
	    return 2 * integrate(r, x + 1.0) - integrate(r, x) - integrate(r, ubound) - y * (ubound - x) + 2 * y * (ubound - x - 1.0);
	}
    }
#endif
}

/* Compute the length of one side of the square grid needed to fully
 * contain a circle of radius r */
double circle_grid_side_length(double r) {
    int rmin = -ceil(r+0.5) + 1, rmax = ceil(r+0.5) - 1;
    return (rmax - rmin + 1);
}

/* Compute the size of the grid needed to fully contain a circle of
 * radius r */
double circle_grid_size(double r) {
    int rmin = -ceil(r+0.5) + 1, rmax = ceil(r+0.5) - 1;
    return (rmax - rmin + 1) * (rmax - rmin + 1);
}


/* Compute the weights to assign grid points that intersect with a
 * circle (a grid point is weighted proportionally to the area of the
 * intersection of the grid point and the circle).
 * r:       the radius of the circle
 * weights: an output array of weights of size (2 * r) * (2 * r) */
void compute_circle_weights(double rgrid, double r, double *weights) {
    int x, y;
    int rmin = -ceil(rgrid+0.5) + 1, rmax = ceil(rgrid+0.5) - 1;
    int w = rmax - rmin + 1;

    for (y = rmin; y <= rmax; y++) {
	for (x = rmin; x <= rmax; x++) {
	    double isect = circle_intersect_pixel(r, 0, 0, x, y);
	    
	    if (isect < 0.0)
		printf("circle parity\n");

	    weights[(y - rmin) * w + x - rmin] = isect;
	}
    }
}

/* Compute the weights to assign grid points that intersect with a
 * ring with outer radius router and inner radius rinner (rinner <=
 * router) , filling in the given array of weights */
void compute_ring_weights(double rgrid, double router, double rinner, double *weights) {
    int x, y;
    int rmin = -ceil(rgrid+0.5) + 1, rmax = ceil(rgrid+0.5) - 1;
    int w = rmax - rmin + 1;

    for (y = rmin; y <= rmax; y++) {
	for (x = rmin; x <= rmax; x++) {
	    double outer = circle_intersect_pixel(router, 0, 0, x, y);
	    double inner = circle_intersect_pixel(rinner, 0, 0, x, y);
	    
	    if (outer < 0.0 || inner > outer) {
		printf("[compute_ring_weights] Ring parity: outer = %0.3f (%0.3f), inner = %0.3f (%0.3f)\n",
		       outer, router, inner, rinner);
	    }

	    weights[(y - rmin) * w + x - rmin] = outer - inner;
	}
    }
}

/* Compute the weights to assign grid points based on a ring-shaped
 * gaussian with radius r and standard deviation sigma.  Returns the
 * total area */
double compute_gaussian_ring_weights(double rgrid, double r, double sigma, double *weights, int normalize) {
    int x, y;
    int rmin = -ceil(rgrid+0.5) + 1, rmax = ceil(rgrid+0.5) - 1;
    int w = rmax - rmin + 1;
    int idx;
    double area = 0.0;

    idx = 0;
    for (y = rmin; y <= rmax; y++) {
	for (x = rmin; x <= rmax; x++, idx++) {
	    double dist = sqrt(x * x + y * y);
	    double val = exp(-(dist - r) * (dist - r) / (2 * sigma * sigma));

	    // printf("%d, %d, %0.3f, %0.3f\n", x, y, dist, val);

	    weights[idx] = val;
	    area += val;
	}
    }

    if (normalize) {
	double inv_area = 1.0 / area;

	area = 0.0;
	for (idx = 0; idx < w * w; idx++) {
	    weights[idx] *= inv_area;
	    area += weights[idx];
	}
    }

    return area;
}
