#include <stdio.h>
#include <string.h>
#include <math.h>
#include "kcorrect.h"

/*
 * eepfit.c
 *
 * Engine for eep
 */

float iterate_lf(float phi[], float lum[], float galaxy_weight[], 
                 float *A, float *B, 
									float *Asum, float *Bsum, int ngals, int nparam, 
									float tolerance);
void set_WH_interp(int interp_choice, float Mdim, float Mbright,
									 float lum[], float M[], int nparam);
void unset_WH_interp(void);
void set_A_lf(float galaxy_mmin[], float galaxy_mmax[],
							float factor[], int ngals, float M[], int nparam, float *A);
void set_B_lf(float galaxy_lum[], float factor[], int ngals, 
							float M[], int nparam, float *B);

float lf_eepfit(int interp_choice,
                float dim,
                float bright,
                float lum[],
                float M[],
                float phi[],
                int num[],
                int nparam,
                float galaxy_M[],
                float galaxy_lum[],
                float factor[],
                int ngals,
                float galaxy_Mmin[],
                float galaxy_Mmax[],
                float galaxy_weight[],
                float A[],
                float B[],
                float Asum[],
                float Bsum[],
                float tolerance,
                int minn)
{
	int i,iM,nOut;
	float like;

	/*
	 * Send interpolation information to W and H routines
	 */
	printf("Send interpolation info to W and H ...\n");
	fflush(stdout);
	set_WH_interp(interp_choice, dim, bright, lum, M, nparam);

	printf("Count objects in bins ...\n");
	fflush(stdout);
	for(i=0;i<nparam;i++) num[i]=0;
	nOut=0;
	for(i=0;i<ngals;i++) {
		iM=(int) floor((float)nparam*(galaxy_M[i]-M[0])/(M[nparam]-M[0]));
		if(iM>=0 && iM<nparam) 
			num[iM]++;
		else 
			nOut++;
	} /* end for i */
	printf("  nOut=%d\n",nOut);
	fflush(stdout);
	if(minn>0) {
		printf("Put any objects in poorly populated bins out of the way ...\n");
		fflush(stdout);
		for(i=0;i<ngals;i++) {
			iM=(int) floor((float)nparam*(galaxy_M[i]-M[0])/(M[nparam]-M[0]));
			if(iM>=0 && iM<nparam) {
				if(num[iM]<minn) {
					galaxy_lum[i]=2.*lum[nparam];
					num[iM]--;
				} /* end if */
			} /* end if */
		} /* end for i */
	} /* end if */

	/* 
	 * Define A and B arrays
	 */
	printf("Define the A array ...\n");
	fflush(stdout);
	set_A_lf(galaxy_Mmin,galaxy_Mmax,factor,ngals,M,nparam,A);
	printf("Define the B array ...\n");
	fflush(stdout);
	set_B_lf(galaxy_lum,factor,ngals,M,nparam,B);

	/*
	 * Do the iteration
	 */
	printf("Iterate to solution ...\n");
	fflush(stdout);
	like=iterate_lf(phi,lum,galaxy_weight,A,B,Asum,Bsum,ngals,nparam,tolerance);
	nOut=0;
	for(i=0;i<nparam;i++)
		if(phi[i]<=0.)
			nOut++;
	printf("%d params\n",nOut);
	fflush(stdout);

	unset_WH_interp();

	return(like);
} /* end main */
