#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "kcorrect.h"

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}

void set_WH_interp(int interp_choice, float Mdim, float Mbright,
									 float lum[], float M[], int nparam);
void unset_WH_interp(void);
float H(float currlum, int indx, float factor);

static float *galaxy_redshift=NULL;
static float *galaxy_absmag=NULL;
static float *galaxy_absmmin=NULL;
static float *galaxy_absmmax=NULL;
static float *lumk=NULL;
static float *absmk=NULL;
static float *phi=NULL;
static float sample_absmmin, sample_absmmax;
static IDL_LONG interp_choice;
static IDL_LONG ngals;
static IDL_LONG nbin;

/* lf_select_eep.c
 * 
 * Calculates selection function based on an EEP luminosity function
 * 
 */ 

float selfunc_lf(float n1,
                 float phi[],
                 float M[],
                 int nparam,
                 float Mmin,
                 float Mmax)
{
	int i;
	float val;
	float Hfactor;
	float lmin,lmax;

	/* absolute mag limits for this field
		 at this distance */
	if(Mmax<M[nparam]) return(0.);
	if(Mmin>M[0]) return(0.);
	if(Mmin<M[nparam]) Mmin=M[nparam];
	if(Mmax>M[0]) Mmax=M[0];
	lmin=pow(10.,-0.4*(Mmax+20.));
	lmax=pow(10.,-0.4*(Mmin+20.));

	/* do the integral */
	val=0.;
	for(i=0;i<nparam;i++) {
		Hfactor=H(lmin,i,1.)-H(lmax,i,1.);
		val+=n1*phi[i]*Hfactor;
	} /* end for i */

	return(val);
} /* end selfunc_lf */

IDL_LONG lf_select_eep(float *in_redshift, 
                       float *in_absmag, 
                       float *in_absmmin, 
                       float *in_absmmax, 
                       IDL_LONG in_ngals,
                       float in_sample_absmmin,
                       float in_sample_absmmax,
                       float *in_absmk, 
                       float *in_phi,
                       float *sel,
                       IDL_LONG in_nbin)
{
  IDL_LONG i;

  interp_choice=1;
  
  /*
   * copy in data 
   */
  ngals=in_ngals;
  galaxy_redshift=(float *) malloc(ngals*sizeof(float));
  galaxy_absmag=(float *) malloc(ngals*sizeof(float));
  galaxy_absmmin=(float *) malloc(ngals*sizeof(float));
  galaxy_absmmax=(float *) malloc(ngals*sizeof(float));
  for(i=0;i<ngals;i++) {
    galaxy_absmag[i]=in_absmag[i];
    galaxy_redshift[i]=in_redshift[i];
    galaxy_absmmin[i]=in_absmmin[i];
    galaxy_absmmax[i]=in_absmmax[i];
  } /* end for i */
  nbin=in_nbin;
  absmk=(float *) malloc((nbin+1)*sizeof(float));
  phi=(float *) malloc(nbin*sizeof(float));
  for(i=0;i<=nbin;i++) absmk[i]=in_absmk[i];
  lumk=(float *) malloc(sizeof(float)*(nbin+1));
  sample_absmmin=in_sample_absmmin;
  sample_absmmax=in_sample_absmmax;

	set_WH_interp(interp_choice, sample_absmmax, sample_absmmin, lumk, absmk, 
                nbin);
  for(i=0;i<nbin;i++) 
    phi[i]=in_phi[i]/(0.4*log(10.)*exp(0.5*(log(lumk[i])+log(lumk[i+1]))));
	for(i=0;i<ngals;i++) 
    sel[i]=selfunc_lf(1.,phi,absmk,nbin,galaxy_absmmin[i],galaxy_absmmax[i]);
	unset_WH_interp();

  FREEVEC(lumk);
  FREEVEC(galaxy_absmag);
  FREEVEC(galaxy_absmmin);
  FREEVEC(galaxy_absmmax);
  FREEVEC(galaxy_redshift);
  FREEVEC(phi);
  FREEVEC(absmk);
  return(1);
} /* end lf_select_eep */
