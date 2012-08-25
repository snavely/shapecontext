/* anniface.c */
/* Interface to the ANN (Approximate Nearest Neighbor) library */

#include "anniface.h"
#include "vector.h"
#include "ANN/ANN.h"

#include <stdio.h>

static double ann_z_weight = 1.0;

/* Set the weight of each step along the z (color) axis */
void set_ann_z_weight(double d) {
    if (ann_z_weight < 0.0) {
        printf("Z-weight must be non-negative\n");
        return;
    }

    ann_z_weight = d;
}

/* Return the current weight on the z axis */
double get_ann_z_weight() {
    return ann_z_weight;
}

/* Create a point search structure from the given array of points */
ANNkd_tree_t *create_ann_tree_3D(int n, v3_t *pts, v3_t axis_weights) {
    int i;
    ANNkd_tree *tree;

    /* Create a new array of points */
    ANNpointArray ann_pts = annAllocPts(n, 3);

    /* Create default axis weights */
    ANNpoint aw = annAllocPt(3, 1.0);

    aw[0] = Vx(axis_weights);
    aw[1] = Vy(axis_weights);
    aw[2] = Vz(axis_weights);

    for (i = 0; i < n; i++) {
        ann_pts[i][0] = Vx(pts[i]);
        ann_pts[i][1] = Vy(pts[i]);
        ann_pts[i][2] = Vz(pts[i]);
    }

    tree = new ANNkd_tree(ann_pts, n, 3, aw);

    return (ANNkd_tree_t *)tree;
}

/* Return the approximate nearest neighbor to q out of the points in
 * the tree */
v3_t query_ann_tree_3D_pt(ANNkd_tree_t *tree, v3_t q, double eps) {
    ANNkd_tree *t = (ANNkd_tree *)tree;
    ANNpoint query = annAllocPt(3);
    ANNidx nn_idx;
    ANNdist dist;
    v3_t ret;

    query[0] = Vx(q) * t->axis_weights[0];
    query[1] = Vy(q) * t->axis_weights[1];
    query[2] = Vz(q) * t->axis_weights[2];

    t->annkSearch(query, 1, &nn_idx, &dist, eps);

    Vx(ret) = t->pts[nn_idx][0] / t->axis_weights[0];
    Vy(ret) = t->pts[nn_idx][1] / t->axis_weights[1];
    Vz(ret) = t->pts[nn_idx][2] / t->axis_weights[2];

    annDeallocPt(query);

    return ret;
}

/* Return the distance to the approximate nearest neighbor to q out of 
 * the points in the tree */
double query_ann_tree_3D_dist(ANNkd_tree_t *tree, v3_t q, double eps) {
    ANNkd_tree *t = (ANNkd_tree *)tree;
    ANNpoint query = annAllocPt(3);
    ANNidx nn_idx;
    ANNdist dist;

    query[0] = Vx(q);
    query[1] = Vy(q);
    query[2] = Vz(q);

    t->annkSearch(query, 1, &nn_idx, &dist, eps);

    annDeallocPt(query);

    return dist;
}

/* Return a number of indices for the k nearest neighbors */
void query_ann_tree_3D_idx(ANNkd_tree_t *tree, v3_t q, double eps, 
			   int num_idxs, int *idxs, double *dists) 
{
    ANNkd_tree *t = (ANNkd_tree *)tree;
    ANNpoint query = annAllocPt(3);
    ANNidx *nn_idx = (ANNidx *) malloc(sizeof(ANNidx) * num_idxs);
    ANNdist *dist = (ANNdist *) malloc(sizeof(ANNdist) * num_idxs);

    query[0] = Vx(q);
    query[1] = Vy(q);
    query[2] = Vz(q);

    t->annkSearch(query, num_idxs, nn_idx, dist, eps);

    for (int i = 0; i < num_idxs; i++) {
	idxs[i] = nn_idx[i];
	dists[i] = dist[i];
    }

    annDeallocPt(query);
    free(dist);
    free(nn_idx);
}

/* Run the deconstructor on the given ANN tree */
void free_ann_tree_3D(ANNkd_tree_t *tree) {
    ANNkd_tree *t = (ANNkd_tree *)tree;

    /* Dealloc the points given to the tree */
    free(t->pts);
    delete t;
}



/* Create a point search structure from the given array of points */
ANNkd_tree_t *create_ann_tree(int n, int d, vec_t *pts, vec_t axis_weights) {
    int i, j;
    ANNkd_tree *tree;

    if (n == 0) {
	printf("[create_ann_tree] Invalid input: n == 0\n");
    }

    if (d == 0) {
	printf("[create_ann_tree] Invalid input: d == 0\n");
    }

    /* Create a new array of points */
    ANNpointArray ann_pts = annAllocPts(n, d);

    /* Create the default axis weights */
    ANNpoint aw = annAllocPt(d, 1.0);

    // printf("[create_ann_tree] Creating tree with %d points of size %d\n", n, d);

    for (i = 0; i < d; i++)
        aw[i] = Vn(axis_weights, i);

    for (i = 0; i < n; i++) {
        for (j = 0; j < d; j++) {
            ann_pts[i][j] = (ANNcoord) pts[i].p[j];
        }
    }

    tree = new ANNkd_tree(ann_pts, n, d, aw);

    return (ANNkd_tree_t *)tree;
}

/* Return the point inside the ANN tree at the given index */
vec_t extract_ann_tree_pt(ANNkd_tree_t *tree, int idx) {
    ANNkd_tree *t = (ANNkd_tree *)tree;
    vec_t ret;
    int i;

    if (idx >= t->n_pts) {
	printf("[extract_ann_tree_pt] ANN point index out of bounds (%d <= %d)\n", t->n_pts, idx);
	return ret;
    }
    
    ret = vec_new(t->dim);
    
    for (i = 0; i < t->dim; i++)
	ret.p[i] = t->pts[idx][i] / t->axis_weights[i];

    return ret;
}

/* Return the approximate nearest neighbor to q out of the points in
 * the tree */
vec_t query_ann_tree_pt_brute_force(ANNkd_tree_t *tree, vec_t q, double *distout) {
    ANNkd_tree *t = (ANNkd_tree *)tree;
    ANNpoint query = annAllocPt(t->dim);
    ANNpoint min_pt = NULL;

    double min_dist = DBL_MAX;
    vec_t ret;
    int i,j;

    if (q.d != t->dim) {
        printf("[query_ann_tree_pt_brute_force] Invalid query\n");
        return ret;
    }
    
    for (i = 0; i < q.d; i++)
        query[i] = q.p[i] * t->axis_weights[i];

    for (i = 0; i < t->n_pts; i++) {
	double dist = 0.0;

	for (j = 0; j < t->dim; j++)
	    dist += (query[j] - t->pts[i][j]) * (query[j] - t->pts[i][j]);

	if (dist < min_dist) {
	    min_dist = dist;
	    min_pt = t->pts[i];
	}
    }

    ret = vec_new(t->dim);

    for (i = 0; i < q.d; i++)
        ret.p[i] = min_pt[i] / t->axis_weights[i];

    annDeallocPt(query);

    if (distout != NULL)
	*distout = min_dist;

    return ret;
}

/* Return the approximate nearest neighbor to q out of the points in
 * the tree */
vec_t query_ann_tree_pt_chamfer(ANNkd_tree_t *tree, vec_t q, double *distout) {
    ANNkd_tree *t = (ANNkd_tree *)tree;
    ANNpoint query = annAllocPt(t->dim);
    ANNpoint min_pt = NULL;

    double min_dist = DBL_MAX;
    vec_t ret;
    int i,xj,yj,xk,yk;
    double wt;

    if (q.d != t->dim) {
        printf("[query_ann_tree_pt_chamfer] Invalid query\n");
        return ret;
    }
    
    for (i = 0; i < q.d; i++)
        query[i] = q.p[i] * t->axis_weights[i];

#define ZWEIGHT (1.0 / (443.405006737633 * 0.10))
    wt = ZWEIGHT * ZWEIGHT;

    for (i = 0; i < t->n_pts; i++) {
	double dist = 0.0;

	for (yj = 0; yj < 5; yj++) {
	    for (xj = 0; xj < 5; xj++) {
		int jidx = 2 + yj * 15 + xj * 3;
		int xmin = -1, ymin = -1;

		double chamfer_min = DBL_MAX;
		double rdiff_min = DBL_MAX, gdiff_min = DBL_MAX, bdiff_min = DBL_MAX;

		for (yk = 0; yk < 5; yk++) {
		    for (xk = 0; xk < 5; xk ++) {
			int kidx = 2 + yk * 15 + xk * 3;

			double rdiff = (query[jidx+0] / t->axis_weights[jidx+0]) - (t->pts[i][kidx+0] / t->axis_weights[kidx+0]);
			double gdiff = (query[jidx+1] / t->axis_weights[jidx+1]) - (t->pts[i][kidx+1] / t->axis_weights[kidx+1]);
			double bdiff = (query[jidx+2] / t->axis_weights[jidx+2]) - (t->pts[i][kidx+2] / t->axis_weights[kidx+2]);

			double chamfer = (xj - xk) * (xj - xk) + (yj - yk) * (yj - yk) +
			    wt * (rdiff * rdiff + gdiff * gdiff + bdiff * bdiff);

			if (chamfer < chamfer_min) {
			    xmin = xk;  ymin = yk;
			    chamfer_min = chamfer;
			    rdiff_min = rdiff;
			    gdiff_min = gdiff;
			    bdiff_min = bdiff;
			}
		    }
		}

		// dist += chamfer_min;
		dist += t->axis_weights[jidx+0] * rdiff_min * rdiff_min + 
		    t->axis_weights[jidx+1] * gdiff_min * gdiff_min + 
		    t->axis_weights[jidx+2] * bdiff_min * bdiff_min;

		if (isnan(dist))
		    printf("isnan!\n");
	    }
	}

	if (dist < min_dist) {
	    min_dist = dist;
	    min_pt = t->pts[i];
	}
    }

    ret = vec_new(t->dim);

    for (i = 0; i < q.d; i++)
        ret.p[i] = min_pt[i] / t->axis_weights[i];

    annDeallocPt(query);

    if (distout != NULL)
	*distout = min_dist;

    return ret;
}

/* Print the points in the tree */
void print_ann_tree(ANNkd_tree_t *tree) {
    ANNkd_tree *t = (ANNkd_tree *)tree;

    printf("Axis weights: ");
    for (int i = 0; i < t->dim; i++) {
	printf("%0.3e ", t->axis_weights[i]);
    }
    printf("\n");
    
    printf("%d points:\n", t->n_pts);
    for (int i = 0; i < t->n_pts; i++) {
	printf("%0.3f %0.3f ...\n", t->pts[i][0], t->pts[i][1]);
    }
}

/* Return the approximate nearest neighbor to q out of the points in
 * the tree */
vec_t query_ann_tree_pt(ANNkd_tree_t *tree, vec_t q, double eps, double *distout, int full, int use_sqrt) {
    ANNkd_tree *t = (ANNkd_tree *)tree;
    ANNpoint query = annAllocPt(t->dim);
    ANNidx nn_idx;
    ANNdist dist;
    vec_t ret;
    int i;

    if (q.d != t->dim) {
        printf("[query_ann_tree_pt] Invalid query\n");
        return ret;
    }
    
    for (i = 0; i < q.d; i++)
        query[i] = q.p[i] * t->axis_weights[i];

    t->annkSearch(query, 1, &nn_idx, &dist, eps);

    ret = vec_new(t->dim);

    for (i = 0; i < q.d; i++)
        ret.p[i] = t->pts[nn_idx][i] / t->axis_weights[i];

    annDeallocPt(query);

    if (distout != NULL) {
	double dist2 = 0.0;
	
	if (full) {
	    for (i = 0; i < q.d; i++) {
		double d = q.p[i] * t->axis_weights[i] - t->pts[nn_idx][i];
		dist2 += d * d;
	    }
	} else {
	    for (i = 2; i < q.d; i++) {
		double d = q.p[i] - t->pts[nn_idx][i] / t->axis_weights[i];
		dist2 += d * d;
	    }	    
	}
	
	if (use_sqrt)
	    *distout = sqrt(dist2);
	else
	    *distout = dist2;

#if 0
	if (!use_sqrt && full && dist2 != dist) {
	    printf("error: %0.3e != %0.3e\n", dist, dist2);
	}
#endif
    }
    
    return ret;
}

/* Return the distance to the approximate nearest neighbor to q out of 
 * the points in the tree */
double query_ann_tree_dist(ANNkd_tree_t *tree, vec_t q, double eps, int *idx) {
    ANNkd_tree *t = (ANNkd_tree *)tree;
    ANNpoint query = annAllocPt(t->dim);
    ANNidx nn_idx;
    ANNdist dist;
    int i;

    if (q.d != t->dim) {
        printf("[query_ann_tree_dist] Invalid query\n");
        return 0.0;
    }

    for (i = 0; i < q.d; i++) 
        query[i] = q.p[i] * t->axis_weights[i];

    t->annkSearch(query, 1, &nn_idx, &dist, eps);

    if (idx != NULL)
	*idx = nn_idx;

    annDeallocPt(query);

    return dist;
}

int query_ann_tree_idx(ANNkd_tree_t *tree, vec_t q, double eps, double *dist_out) 
{
    ANNkd_tree *t = (ANNkd_tree *)tree;
    ANNpoint query = annAllocPt(t->dim);
    ANNidx nn_idx;
    ANNdist dist;
    int i;

    if (q.d != t->dim) {
        printf("[query_ann_tree_dist] Invalid query\n");
        return 0;
    }

    for (i = 0; i < q.d; i++) 
        query[i] = q.p[i] * t->axis_weights[i];

    t->annkSearch(query, 1, &nn_idx, &dist, eps);

    annDeallocPt(query);

    if (dist_out != NULL) {
	*dist_out = dist;
    }

    return nn_idx;    
}

int query_pri_ann_tree_idx(ANNkd_tree_t *tree, vec_t q, double eps, double *dist_out) 
{
    ANNkd_tree *t = (ANNkd_tree *)tree;
    ANNpoint query = annAllocPt(t->dim);
    ANNidx nn_idx;
    ANNdist dist;
    int i;

    if (q.d != t->dim) {
        printf("[query_ann_tree_dist] Invalid query\n");
        return 0;
    }

    for (i = 0; i < q.d; i++) 
        query[i] = q.p[i] * t->axis_weights[i];

    t->annkPriSearch(query, 1, &nn_idx, &dist, eps);

    annDeallocPt(query);

    if (dist_out != NULL) {
	*dist_out = dist;
    }

    return nn_idx;    
}

void query_ann_tree_idx2(ANNkd_tree_t *tree, vec_t q, double eps, 
			 int *idx1, int *idx2, 
			 double *dist_out1, double *dist_out2)
{
    ANNkd_tree *t = (ANNkd_tree *)tree;
    ANNpoint query = annAllocPt(t->dim);
    ANNidx nn_idx[2];
    ANNdist dist[2];
    int i;

    if (q.d != t->dim) {
        printf("[query_ann_tree_dist] Invalid query\n");
        return;
    }

    for (i = 0; i < q.d; i++) 
        query[i] = q.p[i] * t->axis_weights[i];

    t->annkPriSearch(query, 2, nn_idx, dist, eps);

    annDeallocPt(query);

    if (idx1 != NULL) {
	*idx1 = nn_idx[0];
	*idx2 = nn_idx[1];
    }

    if (dist_out1 != NULL) {
	*dist_out1 = dist[0];
	*dist_out2 = dist[1];
    }
}

/* Set the maximum number of points ANN will visit */
void set_ann_points_visit(int n) {
    if (n > 0)
	annMaxPtsVisit(n);
}

/* Return the squared distance to the approximate nearest neighbor to q out of 
 * the points in the tree */
double query_ann_tree_distsq(ANNkd_tree_t *tree, vec_t q, double eps) {
    ANNkd_tree *t = (ANNkd_tree *)tree;
    ANNpoint query = annAllocPt(t->dim);
    ANNidx nn_idx;
    ANNdist dist;
    int i;

    if (q.d != t->dim) {
        printf("[query_ann_tree_dist] Invalid query\n");
        return 0.0;
    }

    for (i = 0; i < q.d; i++) 
        query[i] = q.p[i] * t->axis_weights[i];

    t->annkSearch(query, 1, &nn_idx, &dist, eps);

    annDeallocPt(query);

    return dist;
}

/* Return the distance to the approximate nearest neighbor to q out of 
 * the points in the tree, removing any spatial (x,y) distance */
double query_ann_tree_color_dist(ANNkd_tree_t *tree, vec_t q, double eps) {
    ANNkd_tree *t = (ANNkd_tree *)tree;
    ANNpoint query = annAllocPt(t->dim);
    ANNidx nn_idx;
    ANNdist dist;
    int i;
    double dx, dy;

    if (q.d != t->dim) {
        printf("[query_ann_tree_color_dist] Invalid query\n");
        return 0.0;
    }

    for (i = 0; i < q.d; i++) 
        query[i] = q.p[i] * t->axis_weights[i];

    t->annkSearch(query, 1, &nn_idx, &dist, eps);

    /* Remove the spatial error */
    dx = t->pts[nn_idx][0] - query[0];
    dy = t->pts[nn_idx][1] - query[1];
    dist -= dx * dx;
    dist -= dy * dy;

    annDeallocPt(query);

    return dist;
}

#if 0
/* Change the axis weights on the given ANN tree */
void reweight_ann_tree(ANNkd_tree_t *tree, vec_t axis_weights) {
    (ANNkd_tree *t) = (ANNkd_tree *)tree;

    /* Create the new axis weights */
    ANNpoint aw = annAllocPt(d, 1.0);

    for (i = 0; i < d; i++)
        aw[i] = Vn(axis_weights, i);

    annDeallocPt(tree->axis_weights);
    tree->axis_weights = aw;
    
    /* Now change the weights on all points in the tree */
}
#endif

/* Run the deconstructor on the given ANN tree */
void free_ann_tree(ANNkd_tree_t *tree) {
    ANNkd_tree *t = (ANNkd_tree *)tree;

    /* Dealloc the points given to the tree */
    
    annDeallocPts(t->pts);
    delete t;
}
