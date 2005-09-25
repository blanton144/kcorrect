#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <kcorrect.h>

/*
 * k_binspec.c
 *
 * Bins a spectrum by just integrating over the pixel edges. 
 * Note that this will change a flux density to a flux.
 *
 * Mike Blanton
 * 9/21/2005 */

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}

static float *rb_spectrum=NULL;
static float *rb_lambda=NULL;
static IDL_LONG rb_nl;

/* element of projection --> multplication of (redshifted by cr_z) spectrum
 * and the filter */
float rb_pixel(float lambda) 
{
	float sl,spectrum;
	unsigned long i,ip1;

	k_locate(rb_lambda, rb_nl, lambda, &i);
	if(i>=rb_nl-1 || i<0) return(0.);
	ip1=i+1;
	sl=(lambda-rb_lambda[i])/(rb_lambda[ip1]-rb_lambda[i]);

	spectrum=rb_spectrum[i]+sl*(rb_spectrum[ip1]-rb_spectrum[i]);
	
	return(spectrum);
} /* end filter */

/* Create the rmatrix, a lookup table which speeds analysis */
IDL_LONG k_binspec(float *lambda,
									 float *spectrum,
									 float *newlambda,
									 float *newspectrum,
									 IDL_LONG nl,
									 IDL_LONG nnewl)
{
	float lammin,lammax,scale,currlam;
	IDL_LONG i,l,k,v,indxoff;
	
	/* make local copies of vmatrix */
	rb_spectrum=(float *) malloc(nl*sizeof(float));
	rb_lambda=(float *) malloc(nl*sizeof(float));
	rb_nl=nl;
	for(i=0;i<nl;i++) rb_spectrum[i]=spectrum[i];
	for(i=0;i<nl;i++) rb_lambda[i]=lambda[i];
	
	for(i=0;i<nnewl;i++) 
		newspectrum[i]=k_qromo(rb_pixel,newlambda[i], 
													 newlambda[i+1],k_midpnt);
	
	/* deallocate local vmatrix */
	FREEVEC(rb_spectrum);
	FREEVEC(rb_lambda);
	
	return(1);
} /* end create_r */
