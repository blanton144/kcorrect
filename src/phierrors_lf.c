#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define LN10 2.302585093
#define NR_END 1
#define FREE_ARG char*

void gaussj(float **a, int n, float **b, int m);

/* calculates errors in an EEP luminosity function */

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if (!m) {printf("allocation failure 1 in matrix()"); exit(1);}
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) {printf("allocation failure 2 in matrix()"); exit(1);}
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void phierrors_lf(float phi[],
									float errphi[],
									float covar[],
									float lum[],
									float *A,
									float *B,
									float *Asum,
									float *Bsum,
									int nsample,
									int nparam)
{
	float sum,*dgdphi,**dummy,**inform;
	int ngoodphi,*goodindx,i,j,k,l;

	/* how many good phi values are there? */
	printf("Counting phi values: ");
	fflush(stdout);
	ngoodphi=0;
	for(l=0;l<nparam;l++) 
		if(phi[l]>0.) ngoodphi++;
	printf("%d nonzero phi values out of %d\n",ngoodphi,nparam);
	fflush(stdout);

	/* which are they and what is dgdphi for each? */
	printf("Finding phi values ...\n");
	fflush(stdout);
	goodindx=(int *)malloc(ngoodphi*sizeof(int));
	dgdphi=(float *)malloc(ngoodphi*sizeof(float));
	ngoodphi=0;
	for(l=0;l<nparam;l++) {
		if(phi[l]>0.) {
			goodindx[ngoodphi]=l;
			dgdphi[ngoodphi]=phi[goodindx[ngoodphi]]*(lum[l+1]-lum[l]);
			ngoodphi++;
		} /* end if */
	} /* end for l */
	
	/* make the information matrix */
	printf("Creating information matrix ...\n");
	fflush(stdout);
	inform=matrix(1,ngoodphi+1,1,ngoodphi+1);
	for(i=1;i<=ngoodphi;i++) 
		for(j=1;j<=ngoodphi;j++) {
			sum=0.;
			for(k=0;k<nsample;k++) 
				if(Asum[k]>0. && Bsum[k]>0.) 
					sum+=LN10*LN10*phi[goodindx[i-1]]*phi[goodindx[j-1]]*
						(A[k*nparam+goodindx[i-1]]
						 *A[k*nparam+goodindx[j-1]]/(Asum[k]*Asum[k])
						 -B[k*nparam+goodindx[i-1]]
						 *B[k*nparam+goodindx[j-1]]/(Bsum[k]*Bsum[k]));
			if(i==j) 
				for(k=0;k<nsample;k++) 
					if(Asum[k]>0. && Bsum[k]>0.) 
						sum+=LN10*phi[goodindx[i-1]]*
							(B[k*nparam+goodindx[i-1]]/Bsum[k]
							 -A[k*nparam+goodindx[i-1]]/Asum[k]);
			inform[i][j]=-sum-dgdphi[i-1]*dgdphi[j-1];
		} /* end for ij */
	for(i=1;i<=ngoodphi;i++) {
		inform[i][ngoodphi+1]=-dgdphi[i-1];
		inform[ngoodphi+1][i]=-dgdphi[i-1];
	} /* end for i */
	inform[ngoodphi+1][ngoodphi+1]=0.;

	/* invert the matrix to get covariance matrix */
	printf("Inverting information matrix ...\n");
	fflush(stdout);
	dummy=matrix(1,ngoodphi+1,1,1);
	for(i=1;i<=ngoodphi+1;i++) dummy[i][1]=inform[i][ngoodphi+1];
	gaussj(inform,ngoodphi+1,dummy,1);

	/* extract the errors */
	printf("Extracting errors ...\n");
	fflush(stdout);
	for(i=0;i<nparam;i++) errphi[i]=0.;
	for(i=1;i<=ngoodphi;i++) errphi[goodindx[i-1]]=sqrt(inform[i][i]);
	for(i=0;i<nparam;i++)
		for(j=0;j<nparam;j++)
			covar[i*nparam+j]=0.;
	for(i=0;i<ngoodphi;i++)
		for(j=0;j<ngoodphi;j++)
			covar[goodindx[i]*nparam+goodindx[j]]=inform[i+1][j+1]
				/sqrt(inform[i+1][i+1]*inform[j+1][j+1]);
	
	/* deallocate memory */
  printf("Dealloc ...\n");
  fflush(stdout);
	free_matrix(dummy,1,ngoodphi+1,1,1);
	free_matrix(inform,1,ngoodphi+1,1,ngoodphi+1);
	free((char *)goodindx);
	free((char *)dgdphi);
} /* end phierrors */
							 
