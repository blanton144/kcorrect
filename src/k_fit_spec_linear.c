#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <kcorrect.h>

/*
 * k_fit_spec.c
 *
 * given fluxes, inverse variances, and spectrum does a nonnegative fit
 *
 * Mike Blanton
 * 4/2005 */

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}

IDL_LONG k_fit_spec_linear(float *coeffs,
													 float *flux,
													 float *ivar,
													 float *templates,
													 IDL_LONG nt,
													 IDL_LONG nl,
													 float *chi2,
													 IDL_LONG verbose)
{
  int i,j,jp,k;
	float *invcovar=NULL,*bb=NULL,*xx=NULL,model, *pp=NULL;

	invcovar=(float *) malloc(nt*nt*sizeof(float));
	bb=(float *) malloc(nt*sizeof(float));
	xx=(float *) malloc(nt*sizeof(float));
	pp=(float *) malloc(nt*sizeof(float));

  i=0;

  /*    b. calculate A */
  for(j=0;j<nt;j++) 
    for(jp=j;jp<nt;jp++) 
      invcovar[j*nt+jp]=0;
  for(j=0;j<nt;j++) 
    for(jp=j;jp<nt;jp++) 
      for(k=0;k<nl;k++) 
        invcovar[j*nt+jp]+=templates[j*nl+k]*templates[jp*nl+k]*
          ivar[i*nl+k];
  for(j=0;j<nt;j++) 
    for(jp=0;jp<j;jp++) 
      invcovar[j*nt+jp]=invcovar[jp*nt+j];
  
  /*    c. calculate b */
  for(j=0;j<nt;j++) bb[j]=0.;
  for(j=0;j<nt;j++) 
    for(k=0;k<nl;k++) 
      bb[j]+=templates[j*nl+k]*flux[i*nl+k]*ivar[i*nl+k];

  /* 3. solve */
	k_choldc(invcovar,nt,pp);
	k_cholsl(invcovar,nt,pp,bb,xx);
	
	*chi2=0.;
	for(k=0;k<nl;k++) {
		model=0.;
		for(j=0;j<nt;j++)
			model+=templates[j*nl+k]*xx[j];
		(*chi2)+=(model-flux[k])*(model-flux[k])*ivar[k];
	}


	for(j=0;j<nt;j++)
		coeffs[j]=xx[j];
	
	FREEVEC(invcovar);
	FREEVEC(bb);
	FREEVEC(xx);
	FREEVEC(pp);

	return(1);
} /* end k_fit_spec */
