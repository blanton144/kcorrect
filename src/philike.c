#include <stdio.h>
#include <string.h>
#include <math.h>

float philike(float Asum[],
							float Bsum[],
							int nsample)
{
	float val;
	int k;

	val=0.;
	for(k=0;k<nsample;k++) 
		if(Asum[k]>0. && Bsum[k]>0.) 
			val+=log(Asum[k])-log(Bsum[k]);
	return(2.*val);
} /* end philike */

