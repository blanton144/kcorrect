#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "kcorrect.h"

/*
 * lf_eep.c
 *
 * Use EEP to calculate the shape of the luminosity function.
 * On output, gives unit-normalized luminosity function as well
 * as phi(z) for each galaxy. Also outputs the covariance matrix.
 *
 * Assumes that selection is pretty much just based on luminosity,
 * rather than luminosity and surface brightness together.
 *
 */

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}

static float tolerance=1.e-5;
static IDL_LONG minn=0;
static IDL_LONG ngals;
static IDL_LONG interp_choice=1;
static float sample_absmmin, sample_absmmax;

static float *galaxy_absmag=NULL;
static float *galaxy_absmmin=NULL;
static float *galaxy_absmmax=NULL;
static float *galaxy_redshift=NULL;
static float *galaxy_lum=NULL;
static float *galaxy_weight=NULL;
static float *factor=NULL;

/* 
 * Interpolation properties
 */
static float *lumk=NULL,*absmk=NULL;

/*
 * Luminosity function
 */
static IDL_LONG nbin;
static IDL_LONG *num=NULL;
static float *phi=NULL;
static float *phi_err=NULL;
static float *covar=NULL;
float lf_eepfit(int interp_choice, float dim, float bright, float lum[],
             float M[], float phi[], int num[], int nparam,
             float galaxy_M[], float galaxy_lum[], float factor[], 
             int ngals,
             float galaxy_Mmin[], float galaxy_Mmax[], float galaxy_weight[],
                float A[],
             float B[], float Asum[], float Bsum[], float tolerance,
             int minn);
void phierrors_lf(float phi[], float errphi[], float covar[], float lum[], 
									float *A, float *B, float *Asum, float *Bsum, 
									int ngals, int nparam);

/*
 * A matrix for the iteteration
 */
static float *A=NULL,*B=NULL;
static float *Asum=NULL,*Bsum=NULL;

float eepfit_fit(float zfit);

IDL_LONG lf_eep(float *in_redshift,
                float *in_absmag, 
                float *in_absmmin, 
                float *in_absmmax, 
                IDL_LONG in_ngals,
                float in_sample_absmmin,
                float in_sample_absmmax,
                float *in_absmk, 
                float *in_phi,
                float *in_phi_err,
                float *in_covar,
                IDL_LONG in_nbin,
                IDL_LONG calc_err,
                float *in_weight)
{
  float like;
  int i,j;

  tolerance=1.e-5;
  minn=0;
  interp_choice=1;
  
  /*
   * copy in data 
   */
  ngals=in_ngals;
  galaxy_redshift=(float *) malloc(ngals*sizeof(float));
  galaxy_absmag=(float *) malloc(ngals*sizeof(float));
  galaxy_absmmin=(float *) malloc(ngals*sizeof(float));
  galaxy_absmmax=(float *) malloc(ngals*sizeof(float));
  galaxy_weight=(float *) malloc(ngals*sizeof(float));
  for(i=0;i<ngals;i++) {
    galaxy_redshift[i]=in_redshift[i];
    galaxy_absmag[i]=in_absmag[i];
    galaxy_absmmin[i]=in_absmmin[i];
    galaxy_absmmax[i]=in_absmmax[i];
    galaxy_weight[i]=in_weight[i];
  } /* end for i */
  nbin=in_nbin;
  absmk=(float *) malloc((nbin+1)*sizeof(float));
  phi=(float *) malloc(nbin*sizeof(float));
  phi_err=(float *) malloc(nbin*sizeof(float));
  galaxy_lum=(float *) malloc(sizeof(float)*ngals);
  factor=(float *) malloc(sizeof(float)*ngals);
  sample_absmmin=in_sample_absmmin;
  sample_absmmax=in_sample_absmmax;
  lumk=(float *) malloc(sizeof(float)*(nbin+1));
  num=(IDL_LONG *) malloc(sizeof(IDL_LONG)*(nbin));
  covar=(float *) malloc(sizeof(float)*nbin*nbin);
  
	/*
	 * Allocate memory
	 */
	A=(float *) malloc(ngals*nbin*sizeof(float));
	B=(float *) malloc(ngals*nbin*sizeof(float));
	Asum=(float *) malloc(ngals*sizeof(float));
	Bsum=(float *) malloc(ngals*sizeof(float));

	/* initialize phi */
	for(i=0;i<nbin;i++)
		phi[i]=0.1;

	/*
	 * Define some parameters
	 */
	like=eepfit_fit(0.);

	/*
	 * Calculate the errors 
	 */
  if(calc_err) {
    printf("Calculate errors ...\n");
    fflush(stdout);
    phierrors_lf(phi,phi_err,covar,lumk,A,B,Asum,Bsum,ngals,nbin);
  }

  for(i=0;i<nbin;i++) 
    for(j=0;j<nbin;j++) 
      in_covar[i*nbin+j]=covar[i*nbin+j];
  for(i=0;i<nbin;i++) {
    in_phi[i]=0.4*log(10.)*exp(0.5*(log(lumk[i])+log(lumk[i+1])))*phi[i];
    in_phi_err[i]=phi_err[i];
  }
  for(i=0;i<=nbin;i++) 
    in_absmk[i]=absmk[i];
  
  FREEVEC(covar);
  FREEVEC(galaxy_lum);
  FREEVEC(factor);
  FREEVEC(num);
  FREEVEC(lumk);
  FREEVEC(A);
  FREEVEC(Asum);
  FREEVEC(B);
  FREEVEC(Bsum);
  FREEVEC(galaxy_absmag);
  FREEVEC(galaxy_absmmin);
  FREEVEC(galaxy_absmmax);
  FREEVEC(galaxy_redshift);
  FREEVEC(phi);
  FREEVEC(phi_err);
  FREEVEC(absmk);
  return(1);
} /* end lf_eep */

float eepfit_fit(float zdep)
{
	int i;
	float like,likeoff;
	for(i=0;i<ngals;i++) {
		factor[i]=1.;
		galaxy_lum[i]=pow(10.,-0.4*(galaxy_absmag[i]+20.));
	} /* end for i */
	printf("Fit using zdep=%f\n",zdep);
	fflush(stdout);
	like=lf_eepfit(interp_choice,sample_absmmax,sample_absmmin,lumk,absmk,
                 phi,num,nbin,galaxy_absmag,galaxy_lum,factor,ngals, 
                 galaxy_absmmin,galaxy_absmmax,galaxy_weight, 
                 A,B,Asum,Bsum,tolerance,minn);
	likeoff=0.;
	printf("%f --> %f\n",zdep,like);
	return(like);
} /* eepfit_fit */
