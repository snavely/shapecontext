/* alignlm.c */
/* Routines for aligning pairs of images using Levenberg-Marquardt
 * for optimization */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>

#include "minpack.h"

#include "align.h"
#include "binimage.h"
#include "binvolume.h"
#include "defines.h"
#include "dist.h"
#include "distutil.h"
#include "horn.h"
#include "image.h"
#include "lerp.h"
#include "matrix.h"
#include "metric.h"
#include "resample.h"
#include "transform.h"
#include "transform3D.h"

// #define LM_OUTPUT_INTERMEDIATE
// #define LM_FREE
// #define LM_TRANSLATE_ONLY
#define LM_TRANSLATE_ROTATE

static bimg_t *img1, *img2;
static dtrans_t *img1_dt;

static bvol_t *vol1, *vol2;
static dtrans3D_t *vol1_dt;

static int lm_iter = 0;

static void eval2D(const int *m, const int *n, double *x, double *fvec, int *iflag) {
#ifdef LM_FREE
    trans2D_t T = { { { x[0], x[1], x[2] }, 
                      { x[3], x[4], x[5] },
                      { x[6], x[7], 1.0 } } };
#endif

#ifdef LM_TRANSLATE_ONLY
    trans2D_t T = { { { 1.0, 0.0, x[0] }, 
		      { 0.0, 1.0, x[1] },
		      { 0.0, 0.0, 1.0 } } };
#endif

#ifdef LM_TRANSLATE_ROTATE
    double phi = x[2];
    /* This matrix corresponds to a rotation by phi followed by a
     * translation by x */
    trans2D_t T = { { { cos(phi), -sin(phi), x[0] }, 
		      { sin(phi),  cos(phi), x[1] },
		      {      0.0,       0.0, 1.0 } } };
#endif

    /* Compute the chamfer distance */
    double dist = transform_chamfer_distance(img1, img2, &T, img1_dt);
  
    fvec[0] = dist;
    fvec[1] = fvec[2] = fvec[3] = fvec[4] = fvec[5] = fvec[6] = fvec[7] = 0;
}

trans2D_t *align_bin_image(bimg_t *a, bimg_t *b) {
    /* Find a linear transformation that aligns image B to image A */
    trans2D_t *T = new_identity_transform();

#ifdef LM_FREE
    double xvec[8] = { 1.0, 0.0, 0.0, 
                       0.0, 1.0, 0.0,
                       0.0, 0.0 };
    const int m = 8, n = 8;
#else
    double xvec[3] = { 0.0, 0.0, 0.0 };
    const int m = 8, n = 3;
#endif

    double fvec[8];
    double tol = 1.0e-10;
    int info;
    int iwa[n];
    const int lwa = m * n + 5 * n + m;
    double wa[lwa];
    double min = DBL_MAX, old_min;
    int round = 0;

    img1 = a;
    img2 = b;

    /* Compute the distance transform for a */
    img1_dt = compute_euclidean_DT(a);

    // lmdif1_((void *)eval, &m, &n, xvec, fvec, &tol, &info, iwa, wa, &lwa);

#if 1
    do {
        old_min = min;
        lmdif1_((void *)eval2D, &m, &n, xvec, fvec, &tol, &info, iwa, wa, &lwa);
        eval2D(&m, &n, xvec, fvec, NULL);
        min = fvec[0];
        xvec[0] = rint(xvec[0]);
        xvec[1] = rint(xvec[1]);
        printf("Min after round %d: %f\n", round, min);
        round++;
    } while (min < old_min);
#endif

#ifdef LM_FREE
    for (i = 0; i < 8; i++)
        T->T[i / 3][i % 3] = xvec[i];
#else
    {
        trans2D_t Ttrans = { { { 1.0, 0.0, xvec[0] }, 
                               { 0.0, 1.0, xvec[1] },
                               { 0.0, 0.0, 1.0 } } };
        
        double phi = xvec[2];
    
        trans2D_t Trot = { { { cos(phi), -sin(phi), 0.0 }, 
                             { sin(phi),  cos(phi), 0.0 },
                             {      0.0,       0.0, 1.0 } } };

        printf("angle: %0.3f\n", phi);
        printf("offset: %0.3f, %0.3f\n", xvec[0], xvec[1]);

        T = transform_product(&Ttrans, &Trot);
    }
#endif    
    
    switch (info) {
        case 0:
            printf("Improper input parameters\n");
            break;
        case 1:
            printf("Sum of squares tolerance reached\n");
            break;
        case 2:
            printf("x is within tolerance\n");
            break;
        case 3:
            printf("Sum of squares and x are within tolerance\n");
            break;
        case 4:
            printf("fvec orthogonal\n");
            break;
        case 5:
            printf("max function calls made\n");
            break;
        case 6:
            printf("tolerance is too small (squares)\n");
            break;
        case 7:
            printf("tolerance is too small (x)\n");
            break;
    }

    return T;
}

static void eval3D(const int *m, const int *n, double *x, double *fvec, int *iflag) {
    int i;

#ifdef LM_OUTPUT_INTERMEDIATE
    char fname[64];
    FILE *f;
    img_t *img, *Timg;
    bmp_t *bmp;
#endif

#ifdef LM_FREE
    trans2D_t T = { { { x[0], x[1], x[2] }, 
                      { x[3], x[4], x[5]},
                      { x[6], x[7], 1.0 } } };
		    
#endif

#ifdef LM_TRANSLATE_ONLY
    trans2D_t T = { { { 1.0, 0.0, x[0] }, 
		      { 0.0, 1.0, x[1] },
		      { 0.0, 0.0, 1.0 } } };
#endif

#ifdef LM_TRANSLATE_ROTATE
    double phi = x[2];
    trans2D_t T = { { { cos(phi), -sin(phi), x[0] }, 
		      { sin(phi),  cos(phi), x[1] },
		      {      0.0,       0.0,  1.0 } } };
#endif

    /* Compute the chamfer distance */
    double dist = transform_chamfer_distance_2D_3D(vol1, vol2, &T, vol1_dt);

#if 0
    if (iflag && *iflag == 0) {
#ifdef LM_TRANSLATE_ONLY
        printf("x = (%0.3f, %0.3f)\n", x[0], x[1]);
#endif

#ifdef LM_TRANSLATE_ROTATE
	printf("x = (%0.3f, %0.3f), phi = %0.3f\n", x[0], x[1], x[2]);
#endif

#ifdef LM_OUTPUT_INTERMEDIATE
	img = bvm2img(vol2);
	Timg = resample_img(img, &T);

        sprintf(fname, "sloth%03d.bmp", lm_iter);
        f = fopen(fname, "w");
        if (f == NULL) {
            printf("Cannot open %s for writing\n", fname);
            exit(1);
        }
        bmp = img2bmp(Timg);
        write_bmp(f, bmp);
        fclose(f);

	img_free(img);
	img_free(Timg);
        free_bmp(bmp);
#endif     

	printf("error = %0.3f\n", dist);

	lm_iter++;
    }
#endif

#if 0
    printf("x = (%0.3e, %0.3e)\n", x[0], x[1]);
    printf("error = %0.3f\n", dist);
#endif

    if (isnan(dist))
        printf("Result of eval is NaN\n");

    fvec[0] = dist;
    for (i = 1; i < *m; i++)
        fvec[i] = 0;
}

/* Find a linear transformation that aligns image B to image A.  Uses 
 * Levenberg-Marquardt to do the non-linear optimization. */
trans2D_t *align_bin_volume_LM(bvol_t *a, bvol_t *b) {
    trans2D_t *T = new_identity_transform();

#ifdef LM_FREE
    double xvec[7] = { 1.0, 0.0, 0.0,
		       0.0, 1.0, 0.0,
		       0.0, 0.0 /* 1.0 */ };
    const int m = 7, n = 7;
#endif

#ifdef LM_TRANSLATE_ONLY
    double xvec[2] = { 0.0e-6, 0.0e-6 };
    const int m = 2, n = 2;
#endif

#ifdef LM_TRANSLATE_ROTATE
    double xvec[3] = { 0.0, 0.0, 0.0 };
    const int m = 3, n = 3;
    double phi;
#endif    

    double fvec[m];
    double tol = 1.0e-8;
    double min = DBL_MAX, old_min;
    int round = 0;
    clock_t start, end;

    vol1 = a;
    vol2 = b;

    /* Compute the distance transform for a */
    start = clock();
    vol1_dt = compute_euclidean_DT3D(a);
    end = clock();
    
#if 1
    printf("DT computation time: %0.3fs\n", (float)(end - start) / CLOCKS_PER_SEC);
    /* printf("Initial distance: %0.3f\n", chamfer_distance_3D(a, b, vol1_dt)); */
#endif

    start = clock();
#if 1
    lm_iter = 0;

    do {
        old_min = min;
        lmdif_driver2((void *)eval3D, m, n, xvec, tol);
        eval3D(&m, &n, xvec, fvec, NULL);
        min = fvec[0];

        if (isnan(xvec[0]) || isnan(xvec[1]))
            printf("Result is NaN\n");

#if 0
        xvec[0] = rint(xvec[0]);
        xvec[1] = rint(xvec[1]);
#endif

        /* printf("Min after round %d: %0.3f\n", round, min); */
        round++;
    } while (min < old_min);
#endif
    end = clock();
    
#if 1
    printf("Alignment time: %0.3fs\n", (float)(end - start) / CLOCKS_PER_SEC);
#endif

#ifdef LM_FREE
    for (i = 0; i < 15; i++)
        T->T[i / 4][i % 4] = xvec[i];
#endif

#ifdef LM_TRANSLATE_ONLY
    T->T[0][2] = xvec[0];
    T->T[1][2] = xvec[1];
    /* printf("offset: %0.3f, %0.3f, %0.3f\n", xvec[0], xvec[1], 0.0); */
#endif

#ifdef LM_TRANSLATE_ROTATE
    phi = xvec[2];
    T->T[0][0] = cos(phi);
    T->T[0][1] = -sin(phi);
    T->T[1][0] = sin(phi);
    T->T[1][1] = cos(phi);

    T->T[0][2] = xvec[0];
    T->T[1][2] = xvec[1];
#endif

    /* Cleanup */
    free_DT3D(vol1_dt);

    return T;
}

static int in_range_3D(double *p, double *xrange, double *yrange, double *zrange) {
    if (p[0] < xrange[0] || p[0] > xrange[1])
        return 0;
    if (p[1] < yrange[0] || p[1] > yrange[1])
        return 0;
    if (p[2] < zrange[0] || p[2] > zrange[1])
        return 0;
    
    return 1;
}

static void eval3D_2(const int *m, const int *n, double *x, double *fvec, int *iflag) {
    trans2D_t T = { { { 1.0, 0.0, x[0] }, 
                      { 0.0, 1.0, x[1] },
                      { 0.0, 0.0, 1.0 } } };

    double xrange[2] = { 0.0, (double) vol2->w - 1};
    double yrange[2] = { 0.0, (double) vol2->h - 1};
    double zrange[2] = { 0.0, (double) vol2->d - 1};
    double sum = 0.0;
    int num_hits = 0, num_misses = 0, i;
    double ratio;

    if (iflag && *iflag == 0) {
        printf("x = (%0.3f, %0.3f)\n", x[0], x[1]);
    }
    
#define OUT_OF_BOUNDS_CONSTANT 0.0

    /* Compute the chamfer distance for each point in vol2 */
    if (vol2->type == BVOL_HEIGHTS) {
        int x, y, idx = 0;
        for (y = 0; y < vol2->h; y++) {
            for (x = 0; x < vol2->w; x++) {
                int z = vol2->data.heights[y * vol2->w + x];
                double pnew[3];

                transform_point(&T, (double) x, (double) y, pnew + 0, pnew + 1);
		pnew[2] = (double) z;

                if (!in_range_3D(pnew, xrange, yrange, zrange)) {
                    /* Fudge things by making the contribution constant */
                    fvec[idx] = OUT_OF_BOUNDS_CONSTANT;
                    num_misses++;
                } else {
                    fvec[idx] = sqrt(point_distance_3D(pnew, vol1_dt));
                    sum += fvec[idx];
                    num_hits++;
                }

                idx++;
            }
        }
    } else {
        printf("Cannot yet evaluate bitmapped binary volume\n");
    }

    if (num_hits > 0)
        ratio = sqrt(1.0 / (double) num_hits);
    else 
        ratio = 1.0;

    for (i = 0; i < num_hits + num_misses; i++)
        fvec[i] *= ratio;
    
    /* printf("sum = %0.3f\n", sum); */
}

/* Find a linear transformation that aligns image B to image A.  Uses 
 * Levenberg-Marquardt to do the non-linear optimization. */
trans2D_t *align_bin_volume_LM2(bvol_t *a, bvol_t *b) {
    trans2D_t *T = new_identity_transform();

    double xvec[2] = { 0.0, 0.0 };
    int m = bvol_num_points(b), n = 2;
    double *fvec = malloc(sizeof(double) * m);
    double tol = 1.0e-8;
    double min = DBL_MAX, old_min;
    int round = 0;

    vol1 = a;
    vol2 = b;

    /* Compute the distance transform for a */
    vol1_dt = compute_euclidean_DT3D(a);

    do {
        int i;

        old_min = min;
        lmdif_driver2((void *)eval3D_2, m, n, xvec, tol);
        eval3D_2(&m, &n, xvec, fvec, NULL);

        min = 0.0;
        for (i = 0; i < m; i++) {
            min += fvec[i] * fvec[i];
        }

        if (isnan(xvec[0]) || isnan(xvec[1]))
            printf("Result is NaN\n");

        xvec[0] = rint(xvec[0]);
        xvec[1] = rint(xvec[1]);
        round++;
    } while (min < old_min);
    
    T->T[0][2] = xvec[0];
    T->T[1][2] = xvec[1];
    
    /* Cleanup */
    free_DT3D(vol1_dt);

    return T;
}
