#include <stdio.h>
#include <string.h>
#include <math.h>

float sum_AB_lf(float *AB, float phi[], int k, int nsample, int nparam);
float philike(float Asum[], float Bsum[], int nsample);

/* iterates the EEP equations to maximize the likelihood to determine
   the luminosity function */

float iterate_lf(float phi[],
                 float lum[],
                 float weight[], 
                 float *A,
                 float *B,
                 float *Asum,
                 float *Bsum,
                 int nsample,
									int nparam,
                 float tolerance)
{
	float sum,minphi;
	float bottom,bottomnumer,bottomdenom;
	float top,topnumer,topdenom;
	float newlike,oldlike;
	int i,k,niter,maxiter;

  maxiter=300;

	sum=0.;
	minphi=1.;
	for(i=0;i<nparam;i++)
		if(phi[i]<minphi && phi[i]>0.)
			minphi=phi[i];
	for(i=0;i<nparam;i++)
		if(phi[i]<=0.)
			phi[i]=minphi;
	for(i=0;i<nparam;i++)
		sum+=phi[i]*(lum[i+1]-lum[i]);
	for(i=0;i<nparam;i++)
		phi[i]/=sum;

	oldlike=0.;
	for(k=0;k<nsample;k++) Asum[k]=sum_AB_lf(A,phi,k,nsample,nparam);
	for(k=0;k<nsample;k++) Bsum[k]=sum_AB_lf(B,phi,k,nsample,nparam);
	newlike=philike(Asum,Bsum,nsample);

  niter=0;
	while(fabs(newlike-oldlike)>tolerance && niter<maxiter) {
		printf("like=%f\n",newlike);
		fflush(stdout);

		/* evaluate the next set of phi */
		for(i=0;i<nparam;i++) {
			top=0.;
			for(k=0;k<nsample;k++) {
				topnumer=B[k*nparam+i];
				if(topnumer==0.)
					topdenom=1.;
				else 
					topdenom=Bsum[k]/phi[i];
				top+=weight[k]*topnumer/topdenom;
			} /* end for k */
			if(top==0.) {
				bottom=1.;
			} else {
				bottom=0.;
				for(k=0;k<nsample;k++) {
					bottomnumer=A[k*nparam+i];
					if(bottomnumer==0.)
						bottomdenom=1.;
					else 
						bottomdenom=Asum[k];
					bottom+=weight[k]*bottomnumer/bottomdenom;
				} /* end for k */
			} /* end if..else */
			phi[i]=top/bottom;
		} /* end for i */
		
		/* set constraint */
		sum=0.;
		for(i=0;i<nparam;i++)
			sum+=phi[i]*(lum[i+1]-lum[i]);
		for(i=0;i<nparam;i++)
			phi[i]/=sum;

		/* settings for next iteration */
		for(k=0;k<nsample;k++) Asum[k]=sum_AB_lf(A,phi,k,nsample,nparam);
		for(k=0;k<nsample;k++) Bsum[k]=sum_AB_lf(B,phi,k,nsample,nparam);
		oldlike=newlike;
		newlike=philike(Asum,Bsum,nsample);
		oldlike=oldlike;
    niter++;
	} /* end while */
	printf("like=%f\n",newlike);
	fflush(stdout);
	
	return(newlike);
	
} /* end iterate */
