#include <stdio.h>
#include <string.h>
#include <math.h>

float sum_AB_lf(float *AB,
								float phi[],
								int k,
								int nsample,
								int nparam)
{
	float sum;
	int i;

	sum=0.;
	for(i=0;i<nparam;i++)
		sum+=phi[i]*AB[k*nparam+i];
	return(sum);
} /* end sum_A */
