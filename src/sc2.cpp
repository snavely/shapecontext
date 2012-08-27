/* sc2.c */
/* Main driver for shape context app */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cspnd.h"
// #include "dmap.h"
// #include "dmap-io.h"
#include "image.h"
#include "matrix.h"
#include "pyramid.h"
#include "pyramid-io.h"
#include "vector.h"

#include "match.h"
#include "point.h"
#include "register.h"
#include "shapecontext.h"

typedef struct {
    float x, y;  // feature position
    vec_t v;     // descriptor vector
} feature_t;

void write_features(const std::vector<feature_t> &features,
                    const char *features_file)
{
    FILE *f = fopen(features_file, "w");
    
    if (f == NULL) {
        printf("[write_features] Error opening file %s for writing\n", 
               features_file);
        return;
    }

    int num_features = (int) features.size();
    int dim = 0;
    if (num_features > 0) {
        dim = features[0].v.d;
    }

    fprintf(f, "%d %d\n", num_features, dim);

    for (int i = 0; i < num_features; i++) {
        fprintf(f, "%0.3f %0.3f\n", features[i].x, features[i].y);

        for (int j = 0; j < dim; j++) {
            fprintf(f, "%0.6lf ", features[i].v.p[j]);
        }

        fprintf(f, "\n");
    }

    fclose(f);
}

std::vector<feature_t> 
    compute_shape_context_driver(img_t *img_edge, 
                                 int num_rings, int num_wedges, 
                                 double factor, double sigma, int normalize) 
{
    int i;
    int w = img_edge->w;
    int h = img_edge->h;
    int num_pixels = w * h;
    int descriptor_size = num_rings * num_wedges;
    // vec_t *v = (vec_t *) malloc(sizeof(vec_t) * num_pixels);

    std::vector<vec_t> v;
    v.resize(num_pixels);

    for (i = 0; i < num_pixels; i++) {
	v[i] = vec_new(descriptor_size);
    }

    compute_shape_context(img_edge, num_rings, num_wedges, factor, sigma, 
			  normalize, &(v[0]));

    std::vector<feature_t> features;
    features.resize(num_pixels);
    
    /* Fill in the features data structure */
    int idx = 0;
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++, idx++) {
            features[idx].x = x;
            features[idx].y = y;
            features[idx].v = v[idx];
        }
    }

    return features;
}

#if 0
int main_smooth_maps(int argc, char **argv) 
{
    MotionParams params;

    char *in_image1 = NULL, *in_image2 = NULL;
    img_t *img_edge1 = NULL, *img_edge2 = NULL;
    char *in_dpyr1to2, *in_dpyr2to1;
    img_dist_pyr_t *dpyr1to2, *dpyr2to1;
    img_dist_pyr_t *tmp1, *tmp2;

    char *out_dpyr1to2, *out_dpyr2to1;
    img_dist_pyr_t *dpyr1to2_out, *dpyr2to1_out;
    img_dmap_t *dmap1to2_out, *dmap2to1_out;
    img_dmap_t *dmap1to2_out2, *dmap2to1_out2;
    img_dmap_t *dmap1to2_out3, *dmap2to1_out3;
    img_dmap_t *dmap1to2_out4, *dmap2to1_out4;

    char *out_params;

    if (argc != 9) {
	printf("Usage: %s smoothMaps <edge-image1.bmp> <edge-image2.bmp> "
	       "<dpyr12> <dpyr21> <dpyr12-out> <dpyr21-out> <params.txt>\n", 
	       argv[0]);
	return -1;
    }

    in_image1 = argv[2];
    in_image2 = argv[3];
    in_dpyr1to2 = argv[4];
    in_dpyr2to1 = argv[5];
    out_dpyr1to2 = argv[6];
    out_dpyr2to1 = argv[7];
    out_params = argv[8];

    img_edge1 = img_read_bmp_file(in_image1);
    img_edge2 = img_read_bmp_file(in_image2);

    dpyr1to2 = img_read_distance_pyramid_file(in_dpyr1to2);
    dpyr2to1 = img_read_distance_pyramid_file(in_dpyr2to1);

    dmap1to2_out = img_dmap_new(img_edge1->w, img_edge1->h);
    dmap2to1_out = img_dmap_new(img_edge2->w, img_edge2->h);

    dmap1to2_out2 = img_dmap_new(img_edge1->w, img_edge1->h);
    dmap2to1_out2 = img_dmap_new(img_edge2->w, img_edge2->h);

    dmap1to2_out3 = img_dmap_new(img_edge1->w, img_edge1->h);
    dmap2to1_out3 = img_dmap_new(img_edge2->w, img_edge2->h);

    dmap1to2_out4 = img_dmap_new(img_edge1->w, img_edge1->h);
    dmap2to1_out4 = img_dmap_new(img_edge2->w, img_edge2->h);

    /* Round 1 */
    std::vector<int> inliers;


#if 1
    inliers = EstimateTransform(&(dpyr1to2->dmaps[0]), 
				&(dpyr2to1->dmaps[0]), 
				MotionAffine, 400, 20.0 /*50.0*/, &params,
				dmap1to2_out, dmap2to1_out, false);

    matrix_print(3, 3, params.M);
    printf("[round1] num_inliers = %d\n", inliers.size());
    fflush(stdout);

    tmp1 = dmap2dpyr(dmap1to2_out);
    tmp2 = dmap2dpyr(dmap2to1_out);

    img_write_distance_pyramid_file(tmp1, "tmp0_12.pyr");
    img_write_distance_pyramid_file(tmp2, "tmp0_21.pyr");
#else
    dmap1to2_out = &(dpyr1to2->dmaps[0]);
    dmap2to1_out = &(dpyr2to1->dmaps[0]);
#endif


    /* ------------- */


#if 0

#if 1
    params.num_basis_pts = 64;
    params.lambda = 5.0e3;
    inliers = EstimateTransform(dmap1to2_out, dmap2to1_out,
				MotionThinPlateSpline, 256, 25.0, &params,
				dmap1to2_out2, dmap2to1_out2, false);

    printf("[round2] num_inliers = %d\n", inliers.size());
    fflush(stdout);

    tmp1 = dmap2dpyr(dmap1to2_out2);
    tmp2 = dmap2dpyr(dmap2to1_out2);

    img_write_distance_pyramid_file(tmp1, "tmp1_12.pyr");
    img_write_distance_pyramid_file(tmp2, "tmp1_21.pyr");
#else
    dmap1to2_out2 = dmap1to2_out;
    dmap2to1_out2 = dmap2to1_out;
#endif


#if 1
    params.num_basis_pts = 512; // 256; // 512;
    params.lambda = 1.0e1;
    inliers = EstimateTransform(dmap1to2_out2, dmap2to1_out2,
				MotionThinPlateSpline, 64, 5.0, &params,
				dmap1to2_out3, dmap2to1_out3, true);

    printf("[round3] num_inliers = %d\n", inliers.size());
    fflush(stdout);
#endif

#if 1
    params.num_basis_pts = 2048;
    params.lambda = 5.0e1;
    inliers = EstimateTransform(dmap1to2_out3, dmap2to1_out3,
				MotionThinPlateSpline, 4, 1.0, &params,
				dmap1to2_out4, dmap2to1_out4, true);

    printf("[round4] num_inliers = %d\n", inliers.size());
    fflush(stdout);
#endif

#endif


    /* -------------- */


#if 1
    /* Round 2 */

#if 1
    params.num_basis_pts = 64;
    params.lambda = 5.0e3;
    inliers = EstimateTransform(dmap1to2_out, dmap2to1_out,
				MotionThinPlateSpline, 256, 15.0, &params,
				dmap1to2_out2, dmap2to1_out2, false);

    printf("[round2] num_inliers = %d\n", inliers.size());
    fflush(stdout);
#endif

#if 1
    params.num_basis_pts = 512;
    params.lambda = 1.0e3;
    inliers = EstimateTransform(dmap1to2_out2, dmap2to1_out2,
				MotionThinPlateSpline, 64, 3.0, &params,
				dmap1to2_out3, dmap2to1_out3, true);

    printf("[round3] num_inliers = %d\n", inliers.size());
    fflush(stdout);
#endif

#if 0
    params.num_basis_pts = 2048;
    params.lambda = 5.0e1;
    inliers = EstimateTransform(dmap1to2_out3, dmap2to1_out3,
				MotionThinPlateSpline, 4, 1.0, &params,
				dmap1to2_out4, dmap2to1_out4, true);

    printf("[round4] num_inliers = %d\n", inliers.size());
    fflush(stdout);
#else
    dmap1to2_out4 = dmap1to2_out3;
    dmap2to1_out4 = dmap2to1_out3;
#endif

#endif


    /* ------------ */


#if 0
#if 1
    params.num_basis_pts = 64;
    params.lambda = 1.0e4;
    inliers = EstimateTransform(dmap1to2_out, dmap2to1_out,
				MotionThinPlateSpline, 256, 32.0, &params,
				dmap1to2_out2, dmap2to1_out2, false);

    printf("[round2] num_inliers = %d\n", inliers.size());
#endif

#if 1
    params.num_basis_pts = 512;
    params.lambda = 1.0e3;
    inliers = EstimateTransform(dmap1to2_out2, dmap2to1_out2,
				MotionThinPlateSpline, 64, 8.0, &params,
				dmap1to2_out3, dmap2to1_out3, true);

    printf("[round3] num_inliers = %d\n", inliers.size());
#endif

#if 0
    params.num_basis_pts = 2048;
    params.lambda = 5.0e1;
    inliers = EstimateTransform(dmap1to2_out3, dmap2to1_out3,
				MotionThinPlateSpline, 4, 1.0, &params,
				dmap1to2_out4, dmap2to1_out4, true);

    printf("[round4] num_inliers = %d\n", inliers.size());
    fflush(stdout);
#else
    dmap1to2_out4 = dmap1to2_out3;
    dmap2to1_out4 = dmap2to1_out3;
#endif

#endif

    dpyr1to2_out = dmap2dpyr(dmap1to2_out4);
    dpyr2to1_out = dmap2dpyr(dmap2to1_out4);

    img_write_distance_pyramid_file(dpyr1to2_out, out_dpyr1to2);
    img_write_distance_pyramid_file(dpyr2to1_out, out_dpyr2to1);

    /* Write the parameters */
    FILE *f = fopen(out_params, "w");
    
    fprintf(f, "%d\n", (int) params.basis_pts.size());
    
    for (int i = 0; i < (int) params.basis_pts.size(); i++) {
	fprintf(f, "%0.3f %0.3f\n", 
		params.basis_pts[i].x, params.basis_pts[i].y);
    }

    /* x, y affine */
    fprintf(f, "%0.3f %0.3f %0.3f\n", 
	    params.x_affine[0], params.x_affine[1], params.x_affine[2]);
    fprintf(f, "%0.3f %0.3f %0.3f\n", 
	    params.y_affine[0], params.y_affine[1], params.y_affine[2]);

    /* x, y weights */
    for (int i = 0; i < (int) params.basis_pts.size(); i++) {
	fprintf(f, "%0.3f %0.3f", params.x_weights[i], params.y_weights[i]);
    }

    fclose(f);

    return 0;
}
#endif

int main_initialize_maps(int argc, char **argv) 
{
    // img_t *img_edge1 = NULL, *img_edge2 = NULL;
    img_dmap_t *dmap1to2, *dmap2to1, *dmap1to2_sym, *dmap2to1_sym;
    img_dist_pyr_t *dpyr1to2, *dpyr2to1;

    // vec_t *v1 = NULL, *v2 = NULL;
    // std::vector<feature_t> v1, v2;
    std::vector<feature_t> f1, f2;

    if (argc != 6 && argc != 8) {
	printf("Usage: %s initMaps <edgeimg1.bmp> <edgeimg2.bmp> "
	       "<features1.out> <features2.out> [dpyr12out.pyr] [dpyr21out.pyr]\n", 
	       argv[0]);
	return -1;
    }

    const char *in_image1 = argv[2];
    const char *in_image2 = argv[3];

    // out_dpyr1to2 = argv[4];
    // out_dpyr2to1 = argv[5];
    const char *features1_out = argv[4];
    const char *features2_out = argv[5];

    bool compute_matches = false;
    const char *out_dpyr1to2 = NULL;
    const char *out_dpyr2to1 = NULL;

    if (argc == 8) {
        compute_matches = true;
        out_dpyr1to2 = argv[6];
        out_dpyr2to1 = argv[7];
    }

    img_t *img_edge1 = img_read_bmp_file(in_image1);
    img_t *img_edge2 = img_read_bmp_file(in_image2);
    
#define NUM_RINGS 5 // 6 // 4 // 6 // 14
#define NUM_WEDGES 12 // 8
#define FACTOR 2.4 // 2.0 // 3.0 // 2.0
#define SIGMA 1.0 /* For now, sigma is not used */
#define NORMALIZE 1

    /* First shape context */
    printf("Computing first shape context...\n");
    f1 = compute_shape_context_driver(img_edge1, 
				      NUM_RINGS, NUM_WEDGES, 
                                      FACTOR, SIGMA, NORMALIZE);

    /* Second shape context */
    printf("Computing second shape context...\n");
    f2 = compute_shape_context_driver(img_edge2, 
				      NUM_RINGS, NUM_WEDGES, 
                                      FACTOR, SIGMA, NORMALIZE);

    write_features(f1, features1_out);
    write_features(f2, features2_out);

    if (compute_matches) {
        std::vector<vec_t> v1, v2;
        int f1_size = (int) f1.size();
        int f2_size = (int) f2.size();

        for (int i = 0; i < f1_size; i++) {
            v1[i] = f1[i].v;
        }

        for (int i = 0; i < f2_size; i++) {
            v2[i] = f2[i].v;
        }

        printf("Matching shape contexts [1==>2]...\n");
        dmap1to2 = match_shape_contexts(img_edge1, img_edge2, 
                                        &(v1[0]), &(v2[0]), 
                                        NUM_RINGS, NUM_WEDGES, 
                                        DISTANCE_EUCLIDEAN);

        printf("Matching shape contexts [2==>1]...\n");
        dmap2to1 = match_shape_contexts(img_edge2, img_edge1, 
                                        &(v2[0]), &(v1[0]),
                                        NUM_RINGS, NUM_WEDGES, 
                                        DISTANCE_EUCLIDEAN);

        printf("Pruning matches...\n");

        cspnd_set_similarity_threshold(4.0);
    
        dmap1to2_sym = 
            img_estimate_correspondence(img_edge2, img_edge1, 
                                        dmap2to1, dmap1to2);

        dmap2to1_sym = 
            img_estimate_correspondence(img_edge1, img_edge2, 
                                        dmap1to2, dmap2to1);

        dpyr1to2 = dmap2dpyr(dmap1to2_sym);
        dpyr2to1 = dmap2dpyr(dmap2to1_sym);

        img_write_distance_pyramid_file(dpyr1to2, out_dpyr1to2);
        img_write_distance_pyramid_file(dpyr2to1, out_dpyr2to1);
    }
    
    return 0;
}

#define NUM_FEATURE_POINTS 13

#if 0
int main_ada_maps(int argc, char **argv) 
{
    if (argc != 10) {
	printf("        %s adaMaps <edge-image1.bmp> <edge-image2.bmp> "
	       "<dpyr12-ada> <dpyr21-ada> <dpyr12-in> <dpyr21-in> "
	       "<dpyr12-out> <dpyr21-out>\n", 
	       argv[0]);

	return -1;
    }
    
    char *in_image1 = argv[2];
    char *in_image2 = argv[3];

    char *ada_dpyr1to2 = argv[4];
    char *ada_dpyr2to1 = argv[5];

    char *in_dpyr1to2 = argv[6];
    char *in_dpyr2to1 = argv[7];

    char *out_dpyr1to2 = argv[8];
    char *out_dpyr2to1 = argv[9];
    
    img_t *img_edge1 = img_read_bmp_file(in_image1);
    img_t *img_edge2 = img_read_bmp_file(in_image2);
    
    img_dist_pyr_t *ada1to2 = img_read_distance_pyramid_file(ada_dpyr1to2);
    img_dist_pyr_t *ada2to1 = img_read_distance_pyramid_file(ada_dpyr2to1);
    
    img_dist_pyr_t *in1to2 = img_read_distance_pyramid_file(in_dpyr1to2);
    img_dist_pyr_t *in2to1 = img_read_distance_pyramid_file(in_dpyr2to1);

    int w = img_edge1->w;
    int h = img_edge1->h;

    /* First, compute a TPS warp using the ADA points */
    MotionParams params;

    params.num_basis_pts = NUM_FEATURE_POINTS;
    params.lambda = 1.0e1;

    img_dmap_t *tmp1to2a = img_dmap_new(w, h);
    img_dmap_t *tmp2to1a = img_dmap_new(w, h);

    std::vector<int> inliers;

    inliers = EstimateTransform(&(ada1to2->dmaps[0]), &(ada2to1->dmaps[0]),
				MotionThinPlateSpline, 1, 100.0, &params,
				tmp1to2a, tmp2to1a, true);

    img_dist_pyr_t *tmp1to2a_pyr = dmap2dpyr(tmp1to2a);
    img_dist_pyr_t *tmp2to1a_pyr = dmap2dpyr(tmp2to1a);

    img_write_distance_pyramid_file(tmp1to2a_pyr, "tmp12.pyr");
    img_write_distance_pyramid_file(tmp2to1a_pyr, "tmp21.pyr");

    printf("[AfterInit] num_inliers = %d\n", inliers.size());
    fflush(stdout);
    
    /* Now find the inliers for the shape context points */
    std::vector<point_t> p1, p2;
    std::vector<point_match_t> matches;
    VectorizeDMAP(&(in1to2->dmaps[0]), &(in2to1->dmaps[0]), 
		  p1, p2, matches);

    CountInliers(p1, p2, matches,
		 params.basis_pts,
		 params.x_affine, params.x_weights,
		 params.y_affine, params.y_weights,
		 25.0, inliers);

    printf("[AfterFit] num_inliers = %d\n", inliers.size());

    img_dmap_t *tmp1to2b = img_dmap_new(w, h);
    img_dmap_t *tmp2to1b = img_dmap_new(w, h);

    PruneDMAP(&(in1to2->dmaps[0]), &(in2to1->dmaps[0]),
	      p1, p2, matches, inliers, tmp1to2b, tmp2to1b);

    /* Finally, relax the fit */
    img_dmap_t *tmp1to2c = img_dmap_new(w, h);
    img_dmap_t *tmp2to1c = img_dmap_new(w, h);

    params.num_basis_pts = 512;
    params.lambda = 1.0e3;
    inliers = EstimateTransform(tmp1to2b, tmp2to1b,
				MotionThinPlateSpline, 64, 6.0, &params,
				tmp1to2c, tmp2to1c, true);

    printf("[AfterRelax] num_inliers = %d\n", inliers.size());

    img_dist_pyr_t *out1to2 = dmap2dpyr(tmp1to2c);
    img_dist_pyr_t *out2to1 = dmap2dpyr(tmp2to1c);

    img_write_distance_pyramid_file(out1to2, out_dpyr1to2);
    img_write_distance_pyramid_file(out2to1, out_dpyr2to1);

    return 0;
}
#endif

void printUsage(char *name) 
{
    printf("Usage: %s initMaps <edgeimg1.bmp> <edgeimg2.bmp> "
	   "[dpyr12-out] [dpyr21-out]\n", name);
    printf("       %s smoothMaps <edge-image1.bmp> <edge-image2.bmp> "
	   "<dpyr12-in> <dpyr21-in> <dpyr12-out> <dpyr21-out>"
	   "<params.txt>\n", 
	   name);
    printf("       %s adaMaps <edge-image1.bmp> <edge-image2.bmp "
	   "<dpyr12-ada> <dpyr21-ada> <dpyr12-in> <dpyr21-in> "
	   "<dpyr12-out> <dpyr21-out>\n", 
	   name);
}

int main(int argc, char **argv) 
{
    if (argc < 2) {
	printUsage(argv[0]);
	return -1;
    }
    
    if (strcmp(argv[1], "initMaps") == 0) {
	return main_initialize_maps(argc, argv);
#if 0
    } else if (strcmp(argv[1], "smoothMaps") == 0) {
	return main_smooth_maps(argc, argv);
    } else if (strcmp(argv[1], "adaMaps") == 0) {
	return main_ada_maps(argc, argv);
#endif
    } else {
	printUsage(argv[0]);
    }

    return -1;
}
