/* anniface.h */
/* Interface to the ANN (Approximate Nearest Neighbor) library */

#ifndef __anniface_h__
#define __anniface_h__


#ifdef __cplusplus
extern "C" {
#endif

#include "vector.h"

typedef struct ANNkd_tree_t ANNkd_tree_t;

/* *********** 3D ANNkd_tree functions *********** */

/* Create a point search structure from the given array of 3D points */
ANNkd_tree_t *create_ann_tree_3D(int n, v3_t *pts, v3_t axis_weights);

/* Return the approximate nearest neighbor to q out of the points in
 * the tree */
v3_t query_ann_tree_3D_pt(ANNkd_tree_t *tree, v3_t q, double eps);

/* Return the distance to the approximate nearest neighbor to q out of 
 * the points in the tree */
double query_ann_tree_3D_dist(ANNkd_tree_t *tree, v3_t q, double eps);

/* Return a number of indices for the k nearest neighbors */
void query_ann_tree_3D_idx(ANNkd_tree_t *tree, v3_t q, double eps, 
			   int num_idxs, int *idxs, double *dists);

/* Run the deconstructor on the given ANN tree */
void free_ann_tree_3D(ANNkd_tree_t *tree);


/* *********** General ANNkd_tree functions *********** */

/* Create a point search structure from the given array of points 
 * of dimension `d' */
ANNkd_tree_t *create_ann_tree(int n, int d, vec_t *pts, vec_t axis_weights);

/* Return the point inside the ANN tree at the given index */
vec_t extract_ann_tree_pt(ANNkd_tree_t *tree, int idx);

/* Return the approximate nearest neighbor to q out of the points in
 * the tree */
vec_t query_ann_tree_pt_brute_force(ANNkd_tree_t *tree, vec_t q, double *distout);

/* Return the approximate nearest neighbor to q out of the points in
 * the tree */
vec_t query_ann_tree_pt_chamfer(ANNkd_tree_t *tree, vec_t q, double *distout);
    
/* Return the approximate nearest neighbor to q out of the points in
 * the tree */
vec_t query_ann_tree_pt(ANNkd_tree_t *tree, vec_t q, double eps, double *distout, int full, int use_sqrt);

/* Return the distance to the approximate nearest neighbor to q out of 
 * the points in the tree */
double query_ann_tree_dist(ANNkd_tree_t *tree, vec_t q, double eps, int *idx);

int query_ann_tree_idx(ANNkd_tree_t *tree, vec_t q, 
		       double eps, double *dist_out);
int query_pri_ann_tree_idx(ANNkd_tree_t *tree, vec_t q, 
			   double eps, double *dist_out);
void query_ann_tree_idx2(ANNkd_tree_t *tree, vec_t q, double eps, 
			 int *idx1, int *idx2, 
			 double *dist_out1, double *dist_out2);

/* Return the squared distance to the approximate nearest neighbor to q out of 
 * the points in the tree */
double query_ann_tree_distsq(ANNkd_tree_t *tree, vec_t q, double eps);

/* Return the distance to the approximate nearest neighbor to q out of 
 * the points in the tree, removing any spatial (x,y) distance */
double query_ann_tree_color_dist(ANNkd_tree_t *tree, vec_t q, double eps);
    
/* Run the deconstructor on the given ANN tree */
void free_ann_tree(ANNkd_tree_t *tree);

/* Set the weight of each step along the z (color) axis */
void set_ann_z_weight(double d);

/* Set the maximum number of points ANN will visit */
void set_ann_points_visit(int n);
    
/* Return the current weight on the z axis */
double get_ann_z_weight();

/* Print the points in the tree */
void print_ann_tree(ANNkd_tree_t *tree);

#ifdef __cplusplus
}
#endif

#endif /* __anniface_h__ */
