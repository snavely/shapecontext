/* circle.h */
/* Functions dealing with circles */

#ifndef __circle_h__
#define __circle_h__

/* Compute the area of the intersection of a circle with a unit square
 * (pixel).
 * r:     the radius of the circle
 * xc,yc: the coordinates of the center of the circle
 * x,y:   the coordinates of the lower left corner of the pixel */
double circle_intersect_pixel(double r, int xc, int yc, int xp, int yp);

/* Compute the weights to assign grid points that intersect with a
 * circle (a grid point is weighted proportionally to the area of the
 * intersection of the grid point and the circle).
 * r:       the radius of the circle
 * weights: an output array of weights of size (2 * r) * (2 * r) */
void compute_circle_weights(double rgrid, double r, double *weights);

/* Compute the length of one side of the square grid needed to fully
 * contain a circle of radius r */
double circle_grid_side_length(double r);

/* Compute the size of the grid needed to contain a circle of radius r */
double circle_grid_size(double r);

/* Compute the weights to assign grid points that intersect with a
 * ring with outer radius router and inner radius rinner (rinner <=
 * router) , filling in the given array of weights */
void compute_ring_weights(double rgrid, double router, double rinner, double *weights);

/* Compute the weights to assign grid points based on a ring-shaped
 * gaussian with radius r and standard deviation sigma */
double compute_gaussian_ring_weights(double rgrid, double r, double sigma, double *weights, int normalize);

#endif /* __circle_h__ */
