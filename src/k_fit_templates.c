#include <stdio.h>
#include <string.h>
#include <math.h>
#include <kcorrect.h>

/*
 * fitTemplates.c
 *
 * Updates the eigentemplates given the new coefficients.
 *
 * Mike Blanton
 * 1/2002 */

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}

static double *covar=NULL;
static double *rhs=NULL;

/* finds best fit eigentemplates, given the galaxy properties and 
 * a previous guess at the best coefficients to the templates (in
 * amatrix) */
IDL_LONG k_fit_templates(double *ematrix,    /* eigentemplates */
												 IDL_LONG nt,       /* number of eigentemplates */
												 double *zvals,      /* z interpolation */
												 IDL_LONG nz,
												 double *rmatrix,    /* r matrix */
												 IDL_LONG nk,             /* number of bandpasses */
												 IDL_LONG nb,             /* number of templates */
												 double *amatrix, /* current coefficients */
												 double *galaxy_maggies, /* galaxy maggies */
												 double *galaxy_invvar,
												 double *galaxy_z,
												 IDL_LONG *galaxy_clip,
												 IDL_LONG ngalaxy,
												 double *dmatrix,    /* d matrix */
												 IDL_LONG initialized_dmatrix)
{
	char uplo;
	IDL_LONG i,j,jp,k,b,bp,unity,indx,info,order,use;
	
	/* 
	 * 1. Create system for each galaxy and solve for the appropriate
	 * parameters, and solve it. (Ax=b)
	 */

	/* 1a. create A */
	printf("     (creating A matrix)\n");
	fflush(stdout);
	covar=(double *) malloc(nt*nb*nt*nb*sizeof(double));
	if(!initialized_dmatrix) {
		for(i=0;i<ngalaxy;i++) {
			use=1;
			for(k=0;k<nk;k++)
				if(galaxy_invvar[k+i*nk]==0.)
					use=0;
			if(!galaxy_clip[i] && use) {
				for(k=0;k<nk;k++) 
					for(b=0;b<nb;b++) 
						dmatrix[i*nk*nb+k*nb+b]=
							k_interpolate(galaxy_z[i],&(rmatrix[k*nb*nz+b*nz]),zvals,nz);
			} /* end if */
		} /* end for i */
	} /* end if */
	for(j=0;j<nt;j++) 
		for(b=0;b<nb;b++) 
			for(jp=0;jp<nt;jp++) 
				for(bp=0;bp<nb;bp++) {
					indx=(j*nb+b)*nt*nb+jp*nb+bp;
					covar[indx]=0.;
					/* solving for the best solution */
					for(i=0;i<ngalaxy;i++)  
						if(!galaxy_clip[i]) {
							for(k=0;k<nk;k++) 
								covar[indx]+=
									amatrix[i*nt+j]*amatrix[i*nt+jp]*
									dmatrix[i*nk*nb+k*nb+b]*dmatrix[i*nk*nb+k*nb+bp]
									*galaxy_invvar[k+i*nk];
						} /* end if */
				} /* end for j l jp lp */
	
	/* 1b. create b */
	printf("     (creating b vector)\n");
	fflush(stdout);
	rhs=(double *) malloc(nt*nb*sizeof(double));
	for(j=0;j<nt;j++) 
		for(b=0;b<nb;b++) {
			indx=j*nb+b;
			rhs[indx]=0.;
			for(i=0;i<ngalaxy;i++) {
				use=1;
				for(k=0;k<nk;k++)
					if(galaxy_invvar[k+i*nk]==0.)
						use=0;
				if(!galaxy_clip[i] && use) {
					for(k=0;k<nk;k++) 
						rhs[indx]+=galaxy_maggies[k+i*nk]*amatrix[i*nt+j]*
							k_interpolate(galaxy_z[i],&(rmatrix[k*nb*nz+b*nz]),zvals,nz)
							*galaxy_invvar[k+i*nk];
				} /* end if */
			} /* end for i */
		} /* end for j l */
	
	/* 1c. solve for x */
	printf("     (solving for x vector)\n");
	fflush(stdout);
	uplo='U';
	unity=1;
	order=nb*nt;
	k_dposv__(&uplo,&order,&unity,covar,&order,rhs,&order,&info);
	if(info!=0) {
		printf("info=%d\n",info);
		return(0);
	} /* end if */
	
	/* 1d. store results */
	printf("     (transferring x vector)\n");
	fflush(stdout);
	for(j=0;j<nt;j++) 
		for(b=0;b<nb;b++) 
			ematrix[j*nb+b]=rhs[j*nb+b];
	
	/* 
	 * 2. Clean up
	 */
	FREEVEC(covar);
	FREEVEC(rhs);

	return(1);
} /* end fitTemplates */
	
