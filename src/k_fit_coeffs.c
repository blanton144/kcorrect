#include <stdio.h>
#include <string.h>
#include <math.h>
#include <kcorrect.h>

/*
 * k_fit_coeffs.c
 *
 * Finds the best-fit contribution of each template to each galaxy,
 * given the templates.
 *
 * Mike Blanton
 * 2/2001
 */

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}

static double *cmatrix=NULL;
static double *covar=NULL;
static double *rhs=NULL;

/* fit coefficients, given information about the templates and 
 * the filters (in ematrix and rmatrix) and a set of galaxies with
 * redshifts */
IDL_LONG k_fit_coeffs(double *ematrix,    /* eigentemplates */
											IDL_LONG nt,             /* number of eigentemplates */
											double *zvals,      /* z interpolation */
											IDL_LONG nz,
											double *rmatrix,    /* r matrix */
											IDL_LONG nk,             /* number of bandpasses */
											IDL_LONG nb,             /* number of templates */
											double *amatrix, /* coefficients */
											double *galaxy_flux, /* galaxy fluxes [i][k], 
																							redshifts */
											double *galaxy_invvar,
											double *galaxy_z,
											IDL_LONG ngalaxy)
{
	double currz;
	char uplo;
	IDL_LONG i,j,jp,k,b,indx,unity,info;
	
	/*
	 * 1. Make cmatrix
	 */
#if 0
	printf("    (making cmatrix)\n");
	fflush(stdout);
#endif
	cmatrix=(double *) malloc(ngalaxy*nk*nt*sizeof(double));
	for(i=0;i<ngalaxy;i++) {
		currz=(galaxy_z[i]<0.5*(zvals[0]+zvals[1])) ? 0.5*(zvals[0]+zvals[1]) :
			galaxy_z[i];
		for(k=0;k<nk;k++)
			for(j=0;j<nt;j++) {
				indx=i*nk*nt+k*nt+j;
				cmatrix[indx]=0.;
				for(b=0;b<nb;b++) {
					cmatrix[indx]+=ematrix[j*nb+b]*
						k_interpolate(currz,&(rmatrix[k*nb*nz+b*nz]),zvals,nz);
				}
			} /* end for k j */
	} /* end for i */

	/* 
	 * 2. Create system for each galaxy and solve for the appropriate
	 * parameters (Ax=b)
	 */
	uplo='U';
	unity=1;
#if 0
	printf("    (beginning loop over galaxy)\n");
	fflush(stdout);
#endif
	covar=(double *) malloc(nt*nt*sizeof(double));
	rhs=(double *) malloc(nt*sizeof(double));
	for(i=0;i<ngalaxy;i++) {

		/* 2b. Create A */
		for(j=0;j<nt;j++)
			for(jp=0;jp<nt;jp++) {
				covar[j*nt+jp]=0.;
				for(k=0;k<nk;k++) 
					covar[j*nt+jp]+=cmatrix[i*nk*nt+k*nt+j]
						*cmatrix[i*nk*nt+k*nt+jp]*galaxy_invvar[k+i*nk];
			} /* end for j jp */
																		
		/* 2c. Create b */
		for(j=0;j<nt;j++) {
			rhs[j]=0.;
			for(k=0;k<nk;k++) 
				rhs[j]+=galaxy_flux[k+i*nk]*cmatrix[i*nk*nt+k*nt+j]
						*galaxy_invvar[k+i*nk];
		} /* end for j */

		if(nt==1) {
			amatrix[i]=rhs[0]/covar[0];
		} else {
			/* 2d. Solve for x */
			k_dposv__(&uplo,&nt,&unity,covar,&nt,rhs,&nt,&info);
			if(info!=0) {
				fprintf(stderr,"at galaxy %d, info=%d\n",i,info);
				for(k=0;k<nk;k++) 
					fprintf(stderr,"%e\n",galaxy_flux[k+i*nk]);
				for(k=0;k<nk;k++) 
					fprintf(stderr,"%e\n",galaxy_invvar[k+i*nk]);
				for(j=0;j<nt;j++) 
					for(jp=0;jp<nt;jp++) 
						fprintf(stderr,"%e\n",covar[j*nt+jp]);
				for(j=0;j<nt;j++) 
					fprintf(stderr,"%e\n",rhs[j]);
				for(j=0;j<nt;j++) 
					amatrix[i*nt+j]=0.;
				return(0);
			} /* end if */
			
			/* 2e. Store result */
			for(j=0;j<nt;j++) 
				amatrix[i*nt+j]=rhs[j];
		} /* end if..else */
	} /* end for i */
	
	/* 
	 * 3. Clean up
	 */
	FREEVEC(cmatrix);
	FREEVEC(covar);
	FREEVEC(rhs);

	return(1);
} /* end k_fit_coeffs */
