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

/* Create the rmatrix, a lookup table which speeds analysis */
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
											IDL_LONG verbose)
{
	int i,j,k,jp,ngood;
	float *local_rmatrix=NULL,*invcovar=NULL,*bb=NULL,*pp=NULL,*xx=NULL,offset;
	float *cholesky=NULL;


	local_rmatrix=(float *) malloc(nk*nv*sizeof(float));
	invcovar=(float *) malloc(nv*nv*sizeof(float));
	cholesky=(float *) malloc(nk*nk*sizeof(float));
	bb=(float *) malloc(nv*sizeof(float));
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
		for(j=0;j<ngood;j++) 
			for(jp=0;jp<ngood;jp++)
				cholesky[j*ngood+jp]=invcovar[j*nv+jp];

		/*    c. calculate b */
		for(j=0;j<nv;j++) bb[j]=0.;
		for(j=0;j<nv;j++) 
			for(k=0;k<nk;k++) 
				bb[j]-=local_rmatrix[k*nv+j]*maggies[i*nk+k]*maggies_ivar[i*nk+k];
		
		/*    d. calculate offset to chi2 */
		offset=0.;
		for(k=0;k<nk;k++) 
			offset+=maggies[i*nk+k]*maggies[i*nk+k]*maggies_ivar[i*nk+k];
		offset*=0.5;

		/* 2. pick a reasonable starting point; solve the linear problem
		 * as best as can be for the first few spectra, take the absolute 
		 * value of all components, and then put in a little bit of all the 
		 * other spectra... */
#if 0
		for(j=0;j<2;j++) 
			for(jp=0;jp<2;jp++) 
				printf("ch: %d %d %d %e\n",ngood,j,jp,cholesky[j*nv+jp]);
		for(j=0;j<ngood;j++) 
			for(jp=0;jp<ngood;jp++) 
				if(cholesky[j*ngood+jp]<=0.) 
					printf("ch: %d %d %e\n",j,jp,cholesky[j*ngood+jp]);
		fflush(stdout);
#endif
		ngood=1;
		k_choldc(cholesky,ngood,pp);
		k_cholsl(cholesky,ngood,pp,bb,xx);
		for(j=0;j<ngood;j++) xx[j]=fabs(xx[j]);
		for(j=ngood;j<nv;j++) xx[j]=xx[0]/(float)nv;
		xx[0]/=(float)nv;

		/* 3. run the iteration */
		k_nonneg_solve(xx,invcovar,bb,offset,nv,tolerance,maxiter,niter, 
									 &(chi2[i]),verbose);

		for(j=0;j<nv;j++) coeffs[i*nv+j]=xx[j];
	} /* end for i */

	FREEVEC(local_rmatrix);
	FREEVEC(invcovar);
	FREEVEC(cholesky);
	FREEVEC(bb);
	FREEVEC(pp);
	FREEVEC(xx);

	return(1);
} /* end k_fit_nonneg */
