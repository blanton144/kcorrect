#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <kcorrect.h>

/*
 * k_fit_nonneg.c
 *
 * Given the rmatrix, does a nonnegative fit of templates to data. 
 *
 * Mike Blanton
 * 5/2003 */

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}

IDL_LONG k_fit_nonneg(float *coeffs,
											float *rmatrix,
											IDL_LONG nk,
											IDL_LONG nv,
											float *zvals,
											IDL_LONG nz,
											float *maggies,
											float *maggies_ivar,
											float *redshift,
											IDL_LONG ngalaxy,
											float tolerance,
											IDL_LONG maxiter,
											IDL_LONG *niter,
											float *chi2,
											IDL_LONG verbose, 
                      IDL_LONG dontinit)
{
	int i,j,k,jp,ngood,allpos;
	float *local_rmatrix=NULL,*invcovar=NULL,*bb=NULL,*pp=NULL,*xx=NULL,offset;
	float *cholesky=NULL,*bblin=NULL;

	local_rmatrix=(float *) malloc(nk*nv*sizeof(float));
	invcovar=(float *) malloc(nv*nv*sizeof(float));
	cholesky=(float *) malloc(nv*nv*sizeof(float));
	bb=(float *) malloc(nv*sizeof(float));
	bblin=(float *) malloc(nv*sizeof(float));
	xx=(float *) malloc(nv*sizeof(float));
	pp=(float *) malloc(nv*sizeof(float));
	for(i=0;i<ngalaxy;i++) {

		/* how many good bands? */
		for(ngood=0,k=0;k<nk;k++) 
			if(maggies_ivar[i*nk+k]>0.)
				ngood++;

		/* 1. calculate form of problem */
		/*    a. interpolate rmatrix at this redshift */
		for(j=0;j<nv;j++) 
			for(k=0;k<nk;k++) 
				local_rmatrix[k*nv+j]= 
					k_interpolate_es(redshift[i],&(rmatrix[k*nv*nz+j*nz]),zvals,nz);
		
		/*    b. calculate A */
		for(j=0;j<nv;j++) 
			for(jp=j;jp<nv;jp++) 
				invcovar[j*nv+jp]=0;
		for(j=0;j<nv;j++) 
			for(jp=j;jp<nv;jp++) 
				for(k=0;k<nk;k++) 
					invcovar[j*nv+jp]+=local_rmatrix[k*nv+j]*local_rmatrix[k*nv+jp]*
						maggies_ivar[i*nk+k];
		for(j=0;j<nv;j++) 
			for(jp=0;jp<j;jp++) 
				invcovar[j*nv+jp]=invcovar[jp*nv+j];
		for(j=0;j<nv;j++) 
			for(jp=0;jp<nv;jp++)
				cholesky[j*nv+jp]=invcovar[j*nv+jp];

		/*    c. calculate b */
		for(j=0;j<nv;j++) bb[j]=0.;
		for(j=0;j<nv;j++) 
			for(k=0;k<nk;k++) 
				bb[j]-=local_rmatrix[k*nv+j]*maggies[i*nk+k]*maggies_ivar[i*nk+k];
		for(j=0;j<nv;j++) 
      bblin[j]=-bb[j];
		
		/*    d. calculate offset to chi2 */
		offset=0.;
		for(k=0;k<nk;k++) 
			offset+=maggies[i*nk+k]*maggies[i*nk+k]*maggies_ivar[i*nk+k];
		offset*=0.5;

		/* 2. pick a reasonable starting point; solve the linear problem
		 * as best as can be for the first few spectra, take the absolute 
		 * value of all components, and then put in a little bit of all the 
		 * other spectra... */
    if(!dontinit) {
      if(ngood<=nv) {
        ngood=1;
        k_choldc(cholesky,ngood,pp);
        k_cholsl(cholesky,ngood,pp,bb,xx);
        for(j=0;j<ngood;j++) xx[j]=fabs(xx[j]);
        for(j=ngood;j<nv;j++) xx[j]=xx[0]/(float)nv;
        xx[0]/=(float)nv;
        allpos=0;
      } else {
        k_choldc(cholesky,nv,pp);
        k_cholsl(cholesky,nv,pp,bblin,xx);
        allpos=1;
        for(j=0;j<nv;j++) 
          if(xx[j]<0.) allpos=0;
        chi2[i]=k_nonneg_chi2(xx,invcovar,bb,offset,nv);
      }
    } else {
      for(j=0;j<nv;j++) xx[j]=coeffs[i*nv+j];
      allpos=0;
    }
      
		/* 3. run the iteration */
    if(!allpos) {
      for(j=0;j<nv;j++) if(xx[j]<0.) xx[j]=1.e-3*fabs(xx[j]);
      k_nonneg_solve(xx,invcovar,bb,offset,nv,tolerance,maxiter,niter, 
                     &(chi2[i]),verbose);
    }
		
		for(j=0;j<nv;j++) coeffs[i*nv+j]=xx[j];
	} /* end for i */

	FREEVEC(local_rmatrix);
	FREEVEC(invcovar);
	FREEVEC(cholesky);
	FREEVEC(bb);
	FREEVEC(bblin);
	FREEVEC(pp);
	FREEVEC(xx);

	return(1);
} /* end k_fit_nonneg */
