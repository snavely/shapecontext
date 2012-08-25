#ifndef __lapack_h__
#define __lapack_h__

#ifdef __cplusplus
extern "C" {
#endif

#ifndef F77
#ifdef WIN32
#define F77(x) x ## _
#else
#define F77(func) func ## _
#endif
#endif

/* Routines for inverting matrices */
void F77(dgetrf)(int *m, int *n, double *A, int *lda, int *ipiv, int *info);
void F77(dgetri)(int *n, double *A, int *lda, int *ipiv, double *work, 
	     int *lwork, int *info);

/* Routines for computing the least-squares solution to Ax = b */
void F77(dgelss)(int *m, int *n, int *nrhs, double *A, int *lda, double *b, int *ldb, 
	     double *s, double *rcond, int *rank, double *work, int *lwork, int *info);
void F77(dgelsy)(int *m, int *n, int *nrhs, double *A, int *lda, double *b, int *ldb,
	     int *jpvt, double *rcond, int *rank, double *work, int *lwork, int *info);

void F77(dgesv)(int *n, int *nrhs, double *A, int *lda, int *ipiv, double *b,
	    int *ldb, int *info);

void F77(dgetrs)(char *trans, int *n, int *nrhs, double *A, int *lda, int *ipiv, 
             double *b, int *ldb, int *info);

/* Routine for computing the eigenvalues / eigenvectors of a matrix */
void F77(dgeev)(char *jobvl, char *jobvr, int *n, double *A, int *lda, double *wr, double *wi,
	    double *vl, int *ldvl, double *vr, int *ldvr, double *work, int *lwork,
	    int *info);

/* Routine for singular value decomposition */
void F77(dgesvd)(char *jobu, char *jobvt, int *m, int *n, double *A, int *lda, 
	     double *S, double *U, int *ldu, double *VT, int *ldvt,
	     double *work, int *lwork, int *info);

/* Routine for Cholesky decomposition */
void F77(dpotrf)(char *uplo, int *n, double *A, int *lda, int *info);

/* Routine for QR factorization */
void F77(dgeqrf)(int *m, int *n, double *A, int *lda, double *tau, double *work, 
	     int *lwork, int *info);

/* Routine for RQ factorization */
void F77(dgerqf)(int *m, int *n, double *A, int *lda, double *tau, double *work, 
	     int *lwork, int *info);

#ifdef __cplusplus
}
#endif

#endif /* __lapack_h__ */
