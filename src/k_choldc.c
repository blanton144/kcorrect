#include <math.h>
#include <string.h>
#include <stdio.h>

void k_choldc(float *a, int n, float p[])
{
	int i,j,k;
	float sum;

	for (i=0;i<n;i++) {
		for (j=i;j<n;j++) { 
			for (sum=a[i*n+j],k=i-1;k>=0;k--) sum -= a[i*n+k]*a[j*n+k];
			if (i == j) {
				if (sum <= 0.0) return;
				p[i]=sqrt(sum);
			} else a[j*n+i]=sum/p[i];
		}
	}
}
