#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <kcorrect.h>

/*
 * k_create_r.c
 *
 * Makes the rmatrix. This matrix is the impact of each filter curve
 * on each template basis element.
 *
 * Amazingly, this is the only part of the code in which the 
 * filter curves are necessary.
 *
 * AB magnitudes are assumed. The input spectra should be in units of
 * ergs cm^{-2} s^{-1} A^{-1}
 *
 * The sense of the code here is that it produces AB maggies of a
 * source shifted to z and observed through the given bandpasses.
 * To obtain the AB maggies of the rest-frame spectrum in a set
 * of shifted *bandpasses* requires *dividing* this result by (1+z)
 *
 * Mike Blanton
 * 6/2001 */

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}

/* factor of 10^(0.4*2.41) which converts flux in based on 
 * f_lambda=erg/cm^2/s/A into maggies 
 NO LONGER RELEVANT
#define ABSCALE 9.20450
*/

static double cr_glambdaexp=0.;
static double *cr_glambda=NULL;
static double *cr_sppointer=NULL;
static double *cr_bmatrix=NULL;
static double *cr_lambda=NULL;
static IDL_LONG cr_nb;
static IDL_LONG cr_nl;
static double cr_z;
static IDL_LONG cr_b;
static IDL_LONG cr_filter_n;
static double *cr_filter_lambda=NULL;
static double *cr_filter_pass=NULL;

double cr_filter(double lambda) 
{
	double sl,rflambda,filt,spectrum;
	unsigned long i,ip1;

	k_locate(cr_filter_lambda, cr_filter_n, lambda, &i);
	if(i>=cr_filter_n-1 || i<0) return(0.);
	ip1=i+1;
	sl=(lambda-cr_filter_lambda[i])/(cr_filter_lambda[ip1]-cr_filter_lambda[i]);
	filt=cr_filter_pass[i]+sl*(cr_filter_pass[ip1]-cr_filter_pass[i]);

	rflambda=lambda/(1.+cr_z);
	k_locate(cr_lambda, cr_nl, rflambda, &i);
	if(i>=cr_nl-1 || i<0) return(0.);
	ip1=i+1;
	sl=(rflambda-cr_lambda[i])/(cr_lambda[ip1]-cr_lambda[i]);
	spectrum=cr_sppointer[i]+sl*(cr_sppointer[ip1]-cr_sppointer[i]);
	
	filt=filt*lambda*spectrum;
	return(filt);
} /* end filter */

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
										IDL_LONG maxn,
										double *glambda,
										double glambdaexp)
{
	double lammin,lammax,scale;
	IDL_LONG i,l,k,indxoff;

	/* make local copies of bmatrix */
	cr_bmatrix=(double *) malloc(nb*nl*sizeof(double));
	cr_lambda=(double *) malloc((nl+1)*sizeof(double));
	cr_glambda=(double *) malloc(nl*sizeof(double));
	cr_glambdaexp=glambdaexp;
	cr_nb=nb;
	cr_nl=nl;
	for(i=0;i<cr_nl;i++) 
		cr_glambda[i]=glambda[i];
	for(i=0;i<cr_nb*cr_nl;i++) 
		cr_bmatrix[i]=bmatrix[i];
	for(i=0;i<=cr_nl;i++) 
		cr_lambda[i]=lambda[i];
	
	/* make r matrix */
	for(k=0;k<nk;k++) {
		/* allocate memory */
		cr_filter_n=filter_n[k];
		cr_filter_lambda=(double *) malloc(cr_filter_n*sizeof(double));
		cr_filter_pass=(double *) malloc(cr_filter_n*sizeof(double));
		for(l=0;l<cr_filter_n;l++) {
			cr_filter_lambda[l]=filter_lambda[maxn*k+l];
			cr_filter_pass[l]=filter_pass[maxn*k+l];
		} /* end for l */

		/* create scale factor for filter */
		cr_sppointer=cr_glambda; 
		cr_z=0.;
		scale=pow(10.,-glambdaexp)/k_qromo(cr_filter,cr_filter_lambda[0], 
																			 cr_filter_lambda[cr_filter_n-1], 
																			 k_midpnt);
		if(scale!=scale || fabs(scale)>1.e+20) {
			for(i=0;i<cr_filter_n;i++) {
				printf("%e %e\n",cr_filter_lambda[i],cr_filter_pass[i]); 
				fflush(stdout);
			}
			for(i=0;i<cr_nl;i++) {
				printf("%e %e\n",lambda[i],cr_sppointer[i]);
				fflush(stdout);
			}
		}

		/* loop over redshifts */
		for(i=0;i<nz;i++) {
			indxoff=k*cr_nb*nz+i;
			cr_z=zvals[i];
			lammin=cr_filter_lambda[0];
			lammax=cr_filter_lambda[cr_filter_n-1];
			for(cr_b=0;cr_b<cr_nb;cr_b++) {
				cr_sppointer=&(cr_bmatrix[cr_b*cr_nl]);
				rmatrix[indxoff+cr_b*nz]=scale*k_qromo(cr_filter,lammin,lammax, 
																							 k_midpnt)/(1.+cr_z);
			} /* end for */
		} /* end for i */

		/* deallocate memory */
		FREEVEC(cr_filter_lambda);
		FREEVEC(cr_filter_pass);
		cr_filter_n=-1;
	} /* end for k */

	/* deallocate local bmatrix */
	FREEVEC(cr_bmatrix);
	FREEVEC(cr_lambda);
	FREEVEC(cr_glambda);
	cr_sppointer=NULL;
	cr_nb=-1;
	cr_nl=-1;
	cr_z=-1.;
	
	return(1);
} /* end create_r */
