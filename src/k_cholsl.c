void k_cholsl(float *a, int n, float p[], float b[], float x[])
{
	int i,k;
	float sum;

	for (i=0;i<n;i++) {
		for (sum=b[i],k=i-1;k>=0;k--) sum -= a[i*n+k]*x[k];
		x[i]=sum/p[i];
	}
	for (i=n-1;i>=0;i--) {
		for (sum=x[i],k=i+1;k<n;k++) sum -= a[k*n+i]*x[k];
		x[i]=sum/p[i];
	}
}
