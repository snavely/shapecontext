/* alignDT.h */
/* Align two images using ICP and distance transforms */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "dist.h"
#include "distutil.h"
#include "image.h"
#include "matrix.h"
#include "resample.h"
#include "transform.h"
#include "transform-opt.h"
#include "util.h"

// #define UPSAMPLE_TARGET

static int number_of_calls = 0;

/* Find a linear transformation that aligns image B to image A.
 * Uses ICP and resamples the image after each iteration.  The
 * transform that describes the initial guess is given. */
void align_Timage_ICP_DT(img_t *a, img_t *b, double color_step_size, 
			 transform_class_t tclass, double *Tin, double *Tout, 
			 int symmetric, int verbose) 
{
    dtrans3D_t *a_dt;
    dtrans3D_t *b_dt = NULL;

#ifdef UPSAMPLE_TARGET
    trans2D_t *Tgrow;
    img_t *a_up;
#endif

    double T[9], Ttmp[9], Tnext[9];
    v3_t *r_pts, *l_pts;

    int round;
    double error = DBL_MAX;
    double error_last;

    int w_a, h_a, w_b, h_b;
    int num_pts_a, num_pts_b;

    clock_t start, end;

    w_a = a->w;  h_a = a->h;
    w_b = b->w;  h_b = b->h;

    num_pts_a = w_a * h_a;
    num_pts_b = w_b * h_b;

    memcpy(T, Tin, sizeof(double) * 9);

#define MAX_ROUNDS 512
#define ERROR_EPSILON 1.0e-8
#define EPSILON 1.0e-8

    r_pts = (v3_t *) malloc(sizeof(v3_t) * (num_pts_a + num_pts_b));
    l_pts = (v3_t *) malloc(sizeof(v3_t) * (num_pts_a + num_pts_b));

    set_color_step_size(color_step_size);

    /* Compute a distance transform for A */
    start = clock();
    
    printf("[align_Timage_ICP_DT] Computing DT for A...\n");
#ifdef UPSAMPLE_TARGET
    Tgrow = new_scaling_transform(2.0, 2.0);
    a_up = img_resample_bbox(a, Tgrow);
    a_dt = compute_euclidean_DT3D_huttenlocher(a_up, DT_INTENSITY);
#else
    a_dt = compute_euclidean_DT3D_huttenlocher(a, DT_INTENSITY);
#endif

    if (symmetric) {
	printf("[align_Timage_ICP_DT] Computing DT for B...\n");
	b_dt = compute_euclidean_DT3D_huttenlocher(b, DT_INTENSITY);
    }

    end = clock();

    printf("DT computation took %0.5e s\n", (end - start) / ((double) CLOCKS_PER_SEC));

    start = clock();

    round = 0;
    do {
	int x, y;
	int n = 0;
	char fname[256];
	double Tinv[9];

	matrix_invert(3, T, Tinv);

	error_last = error;
	error = 0.0;

	if (verbose) {
	    img_t *img_tmp, *img_tmp2;
	    trans2D_t *tr = new_transform_vector(T);
	    img_tmp = img_resample_bbox(b, tr);

	    img_tmp2 = img_fuse(2, 1, a, img_tmp);

	    sprintf(fname, "alignDT%02d-%04d.bmp", number_of_calls, round);
	    img_write_bmp_file(img_tmp2, fname);

	    transform_free(tr);
	    img_free(img_tmp);
	    img_free(img_tmp2);
	}
	

	/* Add the B -> A matches */
	for (y = 0; y < h_b; y++) {
	    for (x = 0; x < w_b; x++) {
		/* Transform the point */
		double v[3] = { (double) x, (double) y, 1.0 };
		double Tv[3];
		double x_norm, y_norm;
		color_t c;
		int x_i, y_i, z_i;
		iv3_t disp;
		double dx, dy, dz;

		matrix_product(3, 3, 3, 1, T, v, Tv);

		x_norm = Tv[0] / Tv[2];
		y_norm = Tv[1] / Tv[2];

		/* Bounds check */
		if (x_norm < 0.0 || x_norm >= w_a - 1 ||
		    y_norm < 0.0 || y_norm >= h_a - 1)
		    continue;

		/* Find the displacement vector */
		c = img_get_pixel(b, x, y);

#ifdef UPSAMPLE_TARGET
		x_i = iround(2.0 * x_norm);
		y_i = iround(2.0 * y_norm);
#else
		x_i = iround(x_norm);
		y_i = iround(y_norm);
#endif

		z_i = iround(color_intensity(c));

		disp = getdisp3D(a_dt, x_i, y_i, z_i);

#if 0
		printf("(%0.3f, %0.3f, %d) ==> (%d, %d, %d)\n", 
		       x_norm, y_norm, z_i, 
		       x_i + Vx(disp), y_i + Vy(disp), z_i + Vz(disp));
#endif

		l_pts[n] = v3_new(x_norm, y_norm, 1.0);
#ifdef UPSAMPLE_TARGET
		r_pts[n] = v3_new(0.5 * (x_i + Vx(disp)), 0.5 * (y_i + Vy(disp)), 1.0);
#else
		r_pts[n] = v3_new(x_i + Vx(disp), y_i + Vy(disp), 1.0);
#endif

		dx = Vx(r_pts[n]) - Vx(l_pts[n]);
		dy = Vy(r_pts[n]) - Vy(l_pts[n]);
		dz = Vz(r_pts[n]) - Vz(l_pts[n]);

		error += ((dx * dx + dy * dy) / color_step_size + Vz(disp) * Vz(disp));

		n++;

		if (n > (num_pts_a + num_pts_b)) {
		    printf("[align_Timage_ICP_DT] Error: n > num_pts\n");
		}
	    }
	}
	
	if (symmetric) {
	    /* Now add the A -> B matches */ 
	    for (y = 0; y < h_a; y++) {
		for (x = 0; x < w_a; x++) {
		    /* Transform the point */
		    double v[3] = { (double) x, (double) y, 1.0 };
		    double Tinvv[3];
		    double x_norm, y_norm;
		    color_t c;
		    int x_i, y_i, z_i;
		    iv3_t disp;
		    double dx, dy, dz;
		    double match[3], Tmatch[3];

		    matrix_product(3, 3, 3, 1, Tinv, v, Tinvv);

		    x_norm = Tinvv[0] / Tinvv[2];
		    y_norm = Tinvv[1] / Tinvv[2];

		    /* Bounds check */
		    if (x_norm < 0.0 || x_norm >= w_b - 1 ||
			y_norm < 0.0 || y_norm >= h_b - 1)
			continue;

		    // if (symmetric && (rand() % 5) != 0)
		    //   continue;

		    /* Find the displacement vector */
		    c = img_get_pixel(a, x, y);

		    x_i = iround(x_norm);
		    y_i = iround(y_norm);
		    z_i = iround(color_intensity(c));

		    disp = getdisp3D(b_dt, x_i, y_i, z_i);

		    match[0] = x_i + Vx(disp);
		    match[1] = y_i + Vy(disp);
		    match[2] = 1.0;

		    matrix_product(3, 3, 3, 1, T, match, Tmatch);

		    l_pts[n] = v3_new(Tmatch[0] / Tmatch[2], 
				      Tmatch[1] / Tmatch[2],
				      1.0);
		    r_pts[n] = v3_new(x, y, 1.0);

		    dx = Vx(r_pts[n]) - Vx(l_pts[n]);
		    dy = Vy(r_pts[n]) - Vy(l_pts[n]);
		    dz = Vz(r_pts[n]) - Vz(l_pts[n]);

		    error += sqrt(dx * dx + dy * dy + 
				  color_step_size * color_step_size * dz * dz);

		    n++;

		    if (n > (num_pts_a + num_pts_b)) {
			printf("[align_Timage_ICP_DT] Error: n > num_pts\n");
		    }
		}
	    }
	}

	/* Now solve for the correct transformation */
	find_optimal_transform(tclass, n, r_pts, l_pts, Ttmp);
	matrix_print(3, 3, T);

	matrix_product(3, 3, 3, 3, Ttmp, T, Tnext);
	memcpy(T, Tnext, sizeof(double) * 9);

	error = error / n;

	printf("[align_Timage_ICP_DT] Round %d: E = %0.5e (n = %d)\n", 
	       round, error, n);

	round++;
    } while (fabs(error_last - error) >= EPSILON && error > ERROR_EPSILON && round < MAX_ROUNDS);

    end = clock();
    
    printf("[align_Timage_ICP_DT] Alignment took %0.5e s\n", (end - start) / ((double) CLOCKS_PER_SEC));
    printf("[align_Timage_ICP_DT] ICP converged after %d rounds\n", round);
    printf("[align_Timage_ICP_DT] with a final error of %0.5e\n", error);

    if (error < ERROR_EPSILON) {
	printf("[align_Timage_ICP_DT] Exited early with error of %0.5e\n",
	       error);
    } else {
	printf("[align_Timage_ICP_DT] Exited after maximum number of rounds "
	       "exceeded.\n");
    }

    printf("Final transform:\n");
    matrix_print(3, 3, T);

    memcpy(Tout, T, sizeof(double) * 9);

    free(r_pts);
    free(l_pts);

    free_DT3D(a_dt);

    if (symmetric) 
	free_DT3D(b_dt);

    number_of_calls++;
}


/* Find a linear transformation that aligns image B to image A.
 * Uses ICP and resamples the image after each iteration.  The
 * transform that describes the initial guess is given. */
void align_Timage_ICP_DT2(img_t *a, img_t *b, double color_step_size, 
			  transform_class_t tclass, double *Tin, double *Tout, 
			  int symmetric, int verbose) 
{
    dtrans3D_t *a_dt;
    dtrans3D_t *b_dt = NULL;

#ifdef UPSAMPLE_TARGET
    trans2D_t *Tgrow;
    img_t *a_up;
#endif

    double T[9], Ttmp[9], Tnext[9];
    v3_t *r_pts, *l_pts;

    int round;
    double error = DBL_MAX;
    double error_last;

    int w_a, h_a, w_b, h_b;
    int num_pts_a, num_pts_b;

    clock_t start, end;

    img_t *img_curr;

    int done;

    w_a = a->w;  h_a = a->h;
    w_b = b->w;  h_b = b->h;

    num_pts_a = w_a * h_a;
    num_pts_b = w_b * h_b;

    memcpy(T, Tin, sizeof(double) * 9);

#define MAX_ROUNDS 512
#define ERROR_EPSILON 1.0e-8
#define EPSILON 1.0e-8

    r_pts = (v3_t *) malloc(sizeof(v3_t) * (num_pts_a + num_pts_b));
    l_pts = (v3_t *) malloc(sizeof(v3_t) * (num_pts_a + num_pts_b));

    set_color_step_size(color_step_size);

    /* Compute a distance transform for A */
    start = clock();
    
    printf("[align_Timage_ICP_DT] Computing DT for A...\n");
#ifdef UPSAMPLE_TARGET
    Tgrow = new_scaling_transform(2.0, 2.0);
    a_up = img_resample_bbox(a, Tgrow);
    a_dt = compute_euclidean_DT3D_huttenlocher(a_up, DT_INTENSITY);
#else
    a_dt = compute_euclidean_DT3D_huttenlocher(a, DT_INTENSITY);
#endif

    if (symmetric) {
	printf("[align_Timage_ICP_DT] Computing DT for B...\n");
	b_dt = compute_euclidean_DT3D_huttenlocher(b, DT_INTENSITY);
    }

    end = clock();

    printf("DT computation took %0.5e s\n", (end - start) / ((double) CLOCKS_PER_SEC));

    start = clock();

    round = 0;
    done = 0;

    do {
	int x, y;
	int n = 0;
	char fname[256];
	double Tinv[9];

	trans2D_t *tr = new_transform_vector(T);
	img_curr = img_resample_bbox(b, tr);

	printf("origin -> (%0.3f, %0.3f)\n", 
	       Vx(img_curr->origin), Vy(img_curr->origin));

	matrix_invert(3, T, Tinv);

	error_last = error;
	error = 0.0;

	if (verbose) {
	    img_t *img_tmp, *img_tmp2;
	    trans2D_t *tr = new_transform_vector(T);
	    img_tmp = img_resample_bbox(b, tr);

	    img_tmp2 = img_fuse(2, 1, a, img_tmp);

	    sprintf(fname, "alignDT%02d-%04d.bmp", number_of_calls, round);
	    img_write_bmp_file(img_tmp2, fname);

	    transform_free(tr);
	    img_free(img_tmp);
	    img_free(img_tmp2);
	}
	

	/* Add the B -> A matches */
	for (y = 0; y < img_curr->h; y++) {
	    for (x = 0; x < img_curr->w; x++) {
		/* Transform the point */
		color_t c;
		int x_i, y_i, z_i;
		iv3_t disp;
		double dx, dy, dz;

		if (!img_pixel_is_valid(img_curr, x, y))
		    continue;

		x_i = iround(x + Vx(img_curr->origin));
		y_i = iround(y + Vy(img_curr->origin));

		/* Bounds check */
		if (x_i < 20 || x_i >= w_a - 30 ||
		    y_i < 20 || y_i >= h_a - 30)
		    continue;

		/* Find the displacement vector */
		c = img_get_pixel(img_curr, x, y);
		z_i = iround(color_intensity(c));

		disp = getdisp3D(a_dt, x_i, y_i, z_i);

		l_pts[n] = v3_new(x_i, y_i, 1.0);
		r_pts[n] = v3_new(x_i + Vx(disp), y_i + Vy(disp), 1.0);

		dx = Vx(r_pts[n]) - Vx(l_pts[n]);
		dy = Vy(r_pts[n]) - Vy(l_pts[n]);
		dz = Vz(r_pts[n]) - Vz(l_pts[n]);

		error += ((dx * dx + dy * dy) / color_step_size + Vz(disp) * Vz(disp));

		n++;

		if (n > (num_pts_a + num_pts_b)) {
		    printf("[align_Timage_ICP_DT] Error: n > num_pts\n");
		}
	    }
	}
	
	if (symmetric) {
	    /* Now add the A -> B matches */ 
	    for (y = 0; y < h_a; y++) {
		for (x = 0; x < w_a; x++) {
		    /* Transform the point */
		    double v[3] = { (double) x, (double) y, 1.0 };
		    double Tinvv[3];
		    double x_norm, y_norm;
		    color_t c;
		    int x_i, y_i, z_i;
		    iv3_t disp;
		    double dx, dy, dz;
		    double match[3], Tmatch[3];

		    matrix_product(3, 3, 3, 1, Tinv, v, Tinvv);

		    x_norm = Tinvv[0] / Tinvv[2];
		    y_norm = Tinvv[1] / Tinvv[2];

		    /* Bounds check */
		    if (x_norm < 0.0 || x_norm >= w_b - 1 ||
			y_norm < 0.0 || y_norm >= h_b - 1)
			continue;

		    // if (symmetric && (rand() % 5) != 0)
		    //   continue;

		    /* Find the displacement vector */
		    c = img_get_pixel(a, x, y);

		    x_i = iround(x_norm);
		    y_i = iround(y_norm);
		    z_i = iround(color_intensity(c));

		    disp = getdisp3D(b_dt, x_i, y_i, z_i);

		    match[0] = x_i + Vx(disp);
		    match[1] = y_i + Vy(disp);
		    match[2] = 1.0;

		    matrix_product(3, 3, 3, 1, T, match, Tmatch);

		    l_pts[n] = v3_new(Tmatch[0] / Tmatch[2], 
				      Tmatch[1] / Tmatch[2],
				      1.0);
		    r_pts[n] = v3_new(x, y, 1.0);

		    dx = Vx(r_pts[n]) - Vx(l_pts[n]);
		    dy = Vy(r_pts[n]) - Vy(l_pts[n]);
		    dz = Vz(r_pts[n]) - Vz(l_pts[n]);

		    error += sqrt(dx * dx + dy * dy + 
				  color_step_size * color_step_size * dz * dz);

		    n++;

		    if (n > (num_pts_a + num_pts_b)) {
			printf("[align_Timage_ICP_DT] Error: n > num_pts\n");
		    }
		}
	    }
	}

	/* Now solve for the correct transformation */
	find_optimal_transform(tclass, n, r_pts, l_pts, Ttmp);
	matrix_print(3, 3, T);

	matrix_product(3, 3, 3, 3, Ttmp, T, Tnext);
	memcpy(T, Tnext, sizeof(double) * 9);

	error = error / n;

	printf("[align_Timage_ICP_DT] Round %d: E = %0.5e (n = %d)\n", 
	       round, error, n);

	img_free(img_curr);
	transform_free(tr);

	if (tclass == TRANSFORM_TRANSLATE && 
	    Ttmp[2] * Ttmp[2] + Ttmp[5] * Ttmp[5] < 0.01) 
	    done = 1;

	round++;
    } while (fabs(error_last - error) >= EPSILON && error > ERROR_EPSILON && round < MAX_ROUNDS && !done);

    end = clock();
    
    printf("[align_Timage_ICP_DT] Alignment took %0.5e s\n", (end - start) / ((double) CLOCKS_PER_SEC));
    printf("[align_Timage_ICP_DT] ICP converged after %d rounds\n", round);
    printf("[align_Timage_ICP_DT] with a final error of %0.5e\n", error);

    if (error < ERROR_EPSILON) {
	printf("[align_Timage_ICP_DT] Exited early with error of %0.5e\n",
	       error);
    } else {
	printf("[align_Timage_ICP_DT] Exited after maximum number of rounds "
	       "exceeded.\n");
    }

    printf("Final transform:\n");
    matrix_print(3, 3, T);

    memcpy(Tout, T, sizeof(double) * 9);

    free(r_pts);
    free(l_pts);

    free_DT3D(a_dt);

    if (symmetric) 
	free_DT3D(b_dt);

    number_of_calls++;
}
