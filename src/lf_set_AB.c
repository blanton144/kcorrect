#include <stdio.h>
#include <string.h>
#include <math.h>

float W(float currlum, int indx, float factor);
float H(float currlum, int indx, float factor);

void set_A_lf(float sample_Mmin[],
							float sample_Mmax[],
							float factor[],
							int nsample,
							float M[],
							int nparam,
							float *A)
{
	int i,k;
	float Mmin,Mmax;
	float lmin,lmax;

	for(k=0;k<nsample;k++) {

		/* set absolute mag limits for this overlap
			 and distance */
		Mmin=sample_Mmin[k];
		Mmax=sample_Mmax[k];
		if(Mmin<M[nparam]) Mmin=M[nparam];
		if(Mmax>M[0]) Mmax=M[0];
		lmin=pow(10.,-0.4*(Mmax+20.));
		lmax=pow(10.,-0.4*(Mmin+20.));

		for(i=0;i<nparam;i++) {
			A[k*nparam+i]=(H(lmin,i,factor[k])-H(lmax,i,factor[k]));
		} /* end for i */
	} /* end for k */
} /* end set_A */

void set_B_lf(float sample_lum[],
							float factor[],
							int nsample,
							float M[],
							int nparam,
							float *B)
{
	int i,k;

	for(k=0;k<nsample;k++) {
		for(i=0;i<nparam;i++) {
			B[k*nparam+i]=W(sample_lum[k],i,factor[k]);
		} /* end for i */
	} /* end for k */
} /* end set_B */

