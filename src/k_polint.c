#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <kcorrect.h>
#define NRANSI

void k_polint(float xa[], float ya[], IDL_LONG n, float x, float *y, 
							float *dy)
{
	IDL_LONG i,m,ns=1;
	float den,dif,dift,ho,hp,w;
	float *c,*d;

	dif=fabs(x-xa[1]);
	c=k_vector(1,n);
	d=k_vector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) {
				printf("Error in routine polint");
				exit(1);
			}
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	k_free_vector(d,1,n);
	k_free_vector(c,1,n);
}
#undef NRANSI
