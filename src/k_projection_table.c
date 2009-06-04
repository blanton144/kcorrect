#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <kcorrect.h>

/*
 * k_projection_table.c
 *
 * Makes the rmatrix. This matrix is the impact of each filter curve
 * on a set of spectra.
 *
 * It allows you to shift the bandpasses as well as the spectrum.
 *
 * Amazingly, this is the only part of the code in which the 
 * filter curves are necessary.
 *
 * AB magnitudes are assumed. The input spectra should be in units of
 * ergs cm^{-2} s^{-1} A^{-1}
 *
 * The sense of the code here is that it produces AB maggies of a
 * source shifted to z and observed through the given bandpasses 
 * (where the bandpasses can be shifted by the factor (1+band_shift)
 * blueward if so desired).
 *
 * Mike Blanton
 * 6/2001 */

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}

#define ABFNU 3.631e-20        /* ergs/s/cm^2/Hz */
#define LIGHTSPEED 2.99792e+18 /* Angstrom/s */

static float *cr_vmatrix=NULL;
static float *cr_gmatrix=NULL;
static float *cr_project_matrix=NULL;
static float *cr_lambda=NULL;
static IDL_LONG cr_nv;
static IDL_LONG cr_nl;
static float cr_z;
static IDL_LONG cr_filter_n;
static float *cr_filter_lambda=NULL;
static float *cr_filter_pass=NULL;

/* element of projection --> multplication of (redshifted by cr_z) spectrum
 * and the filter */
float cr_filter(float lambda) 
{
	float sl,rflambda,filt,spectrum;
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
	spectrum=cr_project_matrix[i]
		+sl*(cr_project_matrix[ip1]-cr_project_matrix[i]);
	
	filt=filt*lambda*spectrum/(1.+cr_z);
	return(filt);
} /* end filter */

/* integrate cr_project_matrix*cr_filter_pass */
float cr_integrate() 
{
	float total, filt, sl, mlambda, dl;
	unsigned long ip, il, ilp1;

	/* each range of cr_project_matrix */
	total=0.;
	for(ip=0;ip<cr_nl-1;ip++) {
		/* find the observed wavelength for this bin */
		mlambda= (cr_lambda[ip+1]+cr_lambda[ip])*0.5*(1.+cr_z);
		
		/* find the filter value at that observed wavelength */
		k_locate(cr_filter_lambda, cr_filter_n, mlambda, &il);
		if(il<cr_filter_n-1 && il>=0) {
			ilp1=il+1;
			sl=(mlambda-cr_filter_lambda[il])/ 
				(cr_filter_lambda[ilp1]-cr_filter_lambda[il]);
			filt=cr_filter_pass[il]+sl*(cr_filter_pass[ilp1]-cr_filter_pass[il]);

			/* now add component to integral */
			dl= fabs(cr_lambda[ip+1]-cr_lambda[ip]);
			total+= mlambda*filt*dl*cr_project_matrix[ip];
		}
	}

	return(total);
} /* end filter */

/* Create the rmatrix, a lookup table which speeds analysis */
IDL_LONG k_projection_table(float *rmatrix,
														IDL_LONG nk,
														IDL_LONG nv,
														float *vmatrix,
														float *lambda,
														IDL_LONG nl,
														float *zvals,
														IDL_LONG nz,
														IDL_LONG *filter_n,
														float *filter_lambda,
														float *filter_pass,
														float band_shift,
														IDL_LONG maxn)
{
	float scale,currlam;
	IDL_LONG i,l,k,v,indxoff;

	/* make local copies of vmatrix */
	cr_gmatrix=(float *) malloc(nl*sizeof(float));
	cr_vmatrix=(float *) malloc(nv*nl*sizeof(float));
	cr_lambda=(float *) malloc((nl+1)*sizeof(float));
	cr_nv=nv;
	cr_nl=nl;
	for(i=0;i<nv*nl;i++) 
		cr_vmatrix[i]=vmatrix[i];
	for(i=0;i<=nl;i++) 
		cr_lambda[i]=lambda[i];
	for(i=0;i<nl;i++) {
		currlam=0.5*(cr_lambda[i]+cr_lambda[i+1]);
		cr_gmatrix[i]=ABFNU*LIGHTSPEED/(currlam*currlam);
	}
	
	/* make r matrix */
	for(k=0;k<nk;k++) {
		/* create filter */
		cr_filter_n=filter_n[k];
		cr_filter_lambda=(float *) malloc(cr_filter_n*sizeof(float));
		cr_filter_pass=(float *) malloc(cr_filter_n*sizeof(float));
		for(l=0;l<cr_filter_n;l++) {
			cr_filter_lambda[l]=filter_lambda[maxn*k+l]/(1.+band_shift);
			cr_filter_pass[l]=filter_pass[maxn*k+l];
		} /* end for l */

		/* create scale factor for filter by projection onto the AB source */
		cr_z=0.;
		cr_project_matrix=cr_gmatrix;
		scale=1./cr_integrate();

		/* loop over redshifts */
		for(i=0;i<nz;i++) {
			indxoff=k*nv*nz+i;
			cr_z=zvals[i];
			for(v=0;v<nv;v++) {
				cr_project_matrix=&(cr_vmatrix[v*cr_nl]);
				rmatrix[indxoff+v*nz]=scale*cr_integrate();
			} /* end for */
		} /* end for i */

		/* deallocate memory */
		FREEVEC(cr_filter_lambda);
		FREEVEC(cr_filter_pass);
	} /* end for k */

	/* deallocate local vmatrix */
	FREEVEC(cr_vmatrix);
	FREEVEC(cr_gmatrix);
	FREEVEC(cr_lambda);
	
	return(1);
} /* end create_r */
