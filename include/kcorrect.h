/*
 * If you do not need or want to link your code to IDL, unset LINKTOIDL
 * in the Makefile
 */
#ifdef LINKTOIDL
#include "export.h"
#else
typedef int IDL_LONG;
#endif

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}

/* fit coefficients, given information about the templates and 
 * the filters (in ematrix and rmatrix) and a set of galaxies with
 * redshifts */
IDL_LONG k_fit_coeffs(double *ematrix,    /* eigentemplates */
											IDL_LONG nt,             /* number of eigentemplates */
											double *zvals,      /* z interpolation */
											IDL_LONG nz,
											double *rmatrix,    /* r matrix */
											IDL_LONG nk,             /* number of bandpasses */
											IDL_LONG nb,             /* number of templates */
											double *amatrix, /* coefficients */
											double *galaxy_maggies, /* galaxy maggies */
											double *galaxy_invvar,
											double *galaxy_z,
											double *constraints_amp,
											double *constraints_mean,
											double *constraints_invvar,
											int nconstraints,
											IDL_LONG ngalaxy);

/* finds best fit eigentemplates, given the galaxy properties and 
 * a previous guess at the best coefficients to the templates (in
 * amatrix) */
IDL_LONG k_fit_templates(double *ematrix,    /* eigentemplates */
												 IDL_LONG nt,           /* number of eigentemplates */
												 double *zvals,      /* z interpolation */
												 IDL_LONG nz,
												 double *rmatrix,    /* r matrix */
												 IDL_LONG nk,             /* number of bandpasses */
												 IDL_LONG nb,             /* number of templates */
												 double *amatrix, /* current coefficients */
												 double *galaxy_maggies, /* galaxy maggies */
												 double *galaxy_invvar,
												 double *galaxy_z,
												 IDL_LONG *galaxy_clip,
												 IDL_LONG ngalaxy,
												 double *dmatrix,    /* d matrix */
												 IDL_LONG initialized_dmatrix);

/* calculate the reconstructed fluxes, given the coeffs of the model */
IDL_LONG k_reconstruct_maggies(double *ematrix,    /* eigentemplates */
															IDL_LONG nt, /* number of eigentemplates */
															double *zvals,      /* z interpolation */
															IDL_LONG nz,
															double *rmatrix,    /* r matrix */
															IDL_LONG nk,      /* number of bandpasses */
															IDL_LONG nb,      /* number of templates */
															double *coeffs, /* coefficients */
															double *galaxy_z,
															 double *band_shift,
															double *rec_flux,
															IDL_LONG ngalaxy);

/* fit redshift and coefficients, given information about the templates and 
 * the filters (in ematrix and rmatrix) and a set of galaxies */
IDL_LONG k_fit_photoz(double *ematrix,    /* eigentemplates */
											IDL_LONG nt,             /* number of eigentemplates */
											double *zvals,      /* z interpolation */
											IDL_LONG nz,
											double *rmatrix,    /* r matrix */
											IDL_LONG nk,             /* number of bandpasses */
											IDL_LONG nb,             /* number of templates */
											double *coeffs, /* coefficients */
											double *galaxy_maggies, /* galaxy */
											double *galaxy_invvar,
											double *galaxy_z,
											IDL_LONG ngalaxy);

/* Create the rmatrix, a lookup table which speeds analysis */
IDL_LONG k_create_r(double *rmatrix,
										IDL_LONG nk,
										IDL_LONG nb,
										double *bmatrix,
										double *lambda,
										IDL_LONG nl,
										double *zvals,
										IDL_LONG nz,
										IDL_LONG *filter_n,
										double *filter_lambda,
										double *filter_pass,
										IDL_LONG maxn);

/* open a file and return file pointer */
FILE *k_fileopen(const char *Filename,
								 const char *Mode);

/* load an ascii table in my standard format */
IDL_LONG k_load_ascii_table(double **table,
														IDL_LONG *ndim,
														IDL_LONG **sizes,
														char filename[]);

/* write an ascii table in my standard format */
IDL_LONG k_write_ascii_table(double *table,
														 IDL_LONG ndim,
														 IDL_LONG *sizes,
														 char filename[]);

/* get a single filter */
IDL_LONG k_load_filter(char filename[],
											 IDL_LONG *filter_n,
											 double **filter_lambda,
											 double **filter_pass);

/* load the filters in their standard format */
IDL_LONG k_load_filters(IDL_LONG **filter_n,
												double **filter_lambda,
												double **filter_pass,
												IDL_LONG *maxn,
												IDL_LONG *nk,
												char filterlist[]);

/* interpolation routines */
double k_interpolate(double currpos,
										 double vals[],
										 double pos[],
										 IDL_LONG n);
double k_interpolate_es(double currpos,
												double vals[],
												double pos[],
												IDL_LONG n);

/* NR routines, renamed so as not to overload */
double *k_vector(long nl, long nh);
void k_free_vector(double *v, long nl, long nh);
double k_midpnt(double (*func)(double), double a, double b, IDL_LONG n);
double k_qromo(double (*func)(double), double a, double b,
							 double (*choose)(double(*)(double), double, double, IDL_LONG));
void k_polint(double xa[], double ya[], IDL_LONG n, double x, double *y, 
							double *dy);
void k_locate(double xx[], unsigned long n, double x, unsigned long *j);
double k_zbrent(double (*func)(double), double x1, double x2, double tol);
double k_brent(double ax, double bx, double cx, double (*f)(double), 
							 double tol, double *xmin);

/* lapack stuff to call from C */
void k_dposv__(char *uplo, IDL_LONG *n, IDL_LONG *nrhs, double *a, 
							 IDL_LONG *lda, double *b, IDL_LONG *ldb, IDL_LONG *info);
void k_dposv_(char *uplo, IDL_LONG *n, IDL_LONG *nrhs, double *a, 
						  IDL_LONG *lda, double *b, IDL_LONG *ldb, IDL_LONG *info);

/* qld */
void ql0001_(IDL_LONG *nconstraints, IDL_LONG *neconstraints,
						 IDL_LONG *nconstraints_max, IDL_LONG *nvar, IDL_LONG *nvar_max,
						 IDL_LONG *mp2n, double *cmatrix, double *dmatrix, double *amatrix,
						 double *bmatrix, double *xl, double *xu, double *x,
						 double *lagrange, IDL_LONG *iout, IDL_LONG *ifail,
						 IDL_LONG *iprint, double *war, IDL_LONG *lwar, IDL_LONG *iwar,
						 IDL_LONG *liwar);

#ifdef IRIX
#define k_dposv__ k_dposv_
#endif
