#include <stdio.h>
#include <string.h>
#include <math.h>
#include "utils.h"

/*
 * k_pca_etemplates.c
 *
 * Use PCA on the templates! Assumes they are orthonormal.
 *
 * Mike Blanton
 * 2/2001
 *
 */

static double *covar=NULL;
static double *tmpb=NULL;
static double *tmpa=NULL;
static double *amean=NULL;
static double *ameanold=NULL;
static double *eval=NULL;
static double *fv1=NULL;
static double *f=NULL;

void k_pca_etemplates(double *amatrix,
											IDL_LONG nt,
											IDL_LONG ngalaxy,
											double *ematrix,
											IDL_LONG nb,
											double *bflux,
											IDL_LONG nclip,
											IDL_LONG niter)
{
	double tmp,cliplimit;
	IDL_LONG clip,i,j,jp,b,ierr,unity,order,threent,use,effn;
	char jobz,uplo;

	amean=(double *) malloc((nt-1)*sizeof(double));
	ameanold=(double *) malloc((nt-1)*sizeof(double));
	covar=(double *) malloc((nt-1)*(nt-1)*sizeof(double));
	cliplimit=1.e+10;
	for(j=1;j<nt;j++) 
		amean[j-1]=0.;
	for(clip=0;clip<niter;clip++) {
		
		for(j=1;j<nt;j++) {
			ameanold[j-1]=amean[j-1];
			amean[j-1]=0.;
		} /* end for j */
		
		effn=0;
		for(i=0;i<ngalaxy;i++) {
			use=1;
			for(j=1;j<nt;j++) 
				if(amatrix[i*nt+j]==0. ||
					 fabs(amatrix[i*nt+j]/amatrix[i*nt]-ameanold[j-1])>cliplimit)
					use=0;
			if(use) {
				for(j=1;j<nt;j++) 
					amean[j-1]+=amatrix[i*nt+j]/amatrix[i*nt];
				effn++;
			} /* end if */
		} /* end for i */
		for(j=1;j<nt;j++)
			amean[j-1]/=(double)effn;
		
		for(j=1;j<nt;j++)
			for(jp=1;jp<nt;jp++)
				covar[(j-1)*(nt-1)+(jp-1)]=0.;
		for(i=0;i<ngalaxy;i++) {
				use=1;
				for(j=1;j<nt;j++)
					if(amatrix[i*nt+j]==0. ||
						 fabs(amatrix[i*nt+j]/amatrix[i*nt]-amean[j-1])>cliplimit)
						use=0;
				if(use)
					for(j=1;j<nt;j++)
						for(jp=1;jp<nt;jp++)
							covar[(j-1)*(nt-1)+(jp-1)]+=
								(amatrix[i*nt+j]/amatrix[i*nt]-amean[j-1])
								*(amatrix[i*nt+jp]/amatrix[i*nt]-amean[jp-1]);
			} /* end for i */
			for(j=1;j<nt;j++)
				for(jp=1;jp<nt;jp++)
					covar[(j-1)*(nt-1)+(jp-1)]/=((double)effn-1.);

			cliplimit=0.;
			for(j=1;j<nt;j++)
				cliplimit+=covar[(j-1)*(nt-1)+(j-1)]/((double) nt -1.);
			cliplimit=nclip*sqrt(cliplimit);
			printf("      (clip at %f)\n",cliplimit);
		} /* end for clip */

		/* find its eigenvalues and eigenvectors;
		 * eval lists eigenvalues in ascending order,
		 * evec lists corresponding eigenvectors (evec[ivector][icomponent])*/
		order=nt-1;
		jobz='V';
		uplo='U';
		threent=3*nt;
		eval=(double *) malloc((nt-1)*sizeof(double));
		fv1=(double *) malloc(3*nt*sizeof(double));
		k_dsyev_(&jobz,&uplo,&order,covar,&order,eval,fv1,&threent,&ierr);
		if(ierr!=0) {
			fprintf(stderr,"k_dsyev returned error %d\n",ierr);
		} /* end if */

		/* change ematrix accordingly */
		tmpb=(double *) malloc((nt-1)*nb*sizeof(double));
		for(j=1;j<nt;j++) 
			for(b=0;b<nb;b++)
				tmpb[(j-1)*nb+b]=0.;
		for(j=1;j<nt;j++) 
			for(b=0;b<nb;b++)
				for(jp=1;jp<nt;jp++) 
					tmpb[(j-1)*nb+b]+=ematrix[jp*nb+b]*covar[(j-1)*(nt-1)+(jp-1)];
		for(j=1;j<nt;j++)
			for(b=0;b<nb;b++)
				ematrix[j*nb+b]=tmpb[(j-1)*nb+b];

		/* change amatrix accordingly */
		tmpa=(double *) malloc((nt-1)*sizeof(double));
		for(i=0;i<ngalaxy;i++) {
			for(j=1;j<nt;j++) {
				tmpa[j-1]=0.;
				for(jp=1;jp<nt;jp++) 
					tmpa[j-1]+=amatrix[i*nt+jp]*covar[(jp-1)*(nt-1)+(j-1)];
			} /* end for jp */
			for(j=1;j<nt;j++) 
				amatrix[i*nt+j]=tmpa[j-1];
		} /* end for i */
	} /* end if */
	
	/* clean up */
	FREEVEC(covar);
	FREEVEC(tmpb);
	FREEVEC(tmpa);
	FREEVEC(amean);
	FREEVEC(ameanold);
	FREEVEC(eval);
	FREEVEC(fv1);
	FREEVEC(f);
} /* end pca_templates */
