#include <stdio.h>
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
 * Mike Blanton
 * 6/2001 */

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}

static double *cr_bmatrix=NULL;
static double *cr_lambda=NULL;
static IDL_LONG cr_nb;
static IDL_LONG cr_nl;
static double cr_z;
static IDL_LONG cr_b;
static IDL_LONG cr_filter_n;
static double *cr_filter_lambda=NULL;
static double *cr_filter_pass=NULL;

double filter(double lambda) 
{
	double sl,rflambda,filt,spectrum;
	unsigned long i,ip1,indxoff;

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
	indxoff=cr_b*cr_nl;
	spectrum=cr_bmatrix[indxoff+i]
		+sl*(cr_bmatrix[indxoff+ip1]-cr_bmatrix[indxoff+i]);
	
	filt=filt*lambda*spectrum;
	return(filt);
} /* end filter */

double scalefilter(double lambda) 
{
	double sl,filt;
	unsigned long i,ip1;

	k_locate(cr_filter_lambda, cr_filter_n, lambda, &i);
	if(i>=cr_filter_n-1 || i<0) return(0.);
	ip1=i+1;
	sl=lambda-cr_filter_lambda[i];
	filt=cr_filter_pass[i]+sl*(cr_filter_pass[ip1]-cr_filter_pass[i])
		/(cr_filter_lambda[ip1]-cr_filter_lambda[i]);
	filt=filt/lambda;
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
										IDL_LONG maxn)
{
	double lammin,lammax,scale;
	IDL_LONG i,l,k,indxoff;
	char filename[255];
	FILE *fp;

	/* make local copies of bmatrix */
	cr_bmatrix=(double *) malloc(nb*nl*sizeof(double));
	cr_lambda=(double *) malloc((nl+1)*sizeof(double));
	cr_nb=nb;
	cr_nl=nl;
	for(i=0;i<nb*nl;i++)
		cr_bmatrix[i]=bmatrix[i];
	for(i=0;i<=nl;i++)
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
		scale=1./k_qromo(scalefilter,cr_filter_lambda[0],
										 cr_filter_lambda[cr_filter_n-1],k_midpnt);

		/* loop over redshifts */
		for(i=0;i<nz;i++) {
			indxoff=k*nb*nz+i;
			cr_z=zvals[i];
			lammin=cr_filter_lambda[0];
			lammax=cr_filter_lambda[cr_filter_n-1];
			for(cr_b=0;cr_b<nb;cr_b++) 
				rmatrix[indxoff+cr_b*nz]=scale*k_qromo(filter,lammin,lammax,k_midpnt)/
					(1.+cr_z);
		} /* end for i */

		/* deallocate memory */
		FREEVEC(cr_filter_lambda);
		FREEVEC(cr_filter_pass);
	} /* end for k */

	/* deallocate local bmatrix */
	FREEVEC(cr_bmatrix);
	FREEVEC(cr_lambda);
	
	return(1);
} /* end create_r */
