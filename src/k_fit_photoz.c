#include <stdio.h>
#include <string.h>
#include <math.h>
#include <kcorrect.h>

/*
 * k_fit_photoz.c
 *
 * Finds the best-fit contribution of each template, and 
 * the best redshift estimate, to each galaxy, given the templates.
 *
 * Mike Blanton
 * 1/2002
 */

#define ZRES 0.10
#define TOL 0.005
#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}

static double *pz_ematrix=NULL;
static double *pz_rmatrix=NULL;
static double *pz_zvals=NULL;
static double *pz_coeffs=NULL;
static double *pz_model_flux=NULL;
static double *pz_galaxy_flux=NULL;
static double *pz_galaxy_invvar=NULL;
static double *pz_coeffspnorm=NULL;
static double pz_z,pz_linepos;
static IDL_LONG pz_nt,pz_nz,pz_nk,pz_nb,pz_puse,pz_np,pz_p;

double pz_calc_chi2(double z) 
{
	double *model_flux,chi2;
	IDL_LONG k;
	
	k_model_fluxes(pz_ematrix,pz_nt,pz_zvals,pz_nz,pz_rmatrix,
								 pz_nk,pz_nb,pz_coeffs,&z,pz_model_flux,1);
	chi2=0.;
	for(k=0;k<pz_nk;k++)
		chi2+=(pz_model_flux[k]-pz_galaxy_flux[k])
			*(pz_model_flux[k]-pz_galaxy_flux[k])*pz_galaxy_invvar[k];
	
	return(chi2);
} /* end pz_calc_chi2 */

/* calculates chi2 */
double pz_fit_coeffs(double z) 
{
	double chi2,az,bz,cz;
	IDL_LONG i,k;
	
	/* fit coeffs */
	k_fit_coeffs(pz_ematrix,pz_nt,pz_zvals,pz_nz,pz_rmatrix,
							 pz_nk,pz_nb,pz_coeffs,pz_galaxy_flux,pz_galaxy_invvar,
							 &z,1);
	
	chi2=pz_calc_chi2(z);
	return(chi2);
} /* end pz_fit_coeffs */

/* fit redshift and coefficients, given information about the templates and 
 * the filters (in ematrix and rmatrix) and a set of galaxies */
IDL_LONG k_fit_photoz(double *ematrix,    /* eigentemplates */
											IDL_LONG nt,             /* number of eigentemplates */
											double *zvals,      /* z interpolation */
											IDL_LONG nz,
											double *rmatrix,    /* r matrix */
											IDL_LONG nk,             /* number of bandpasses */
											IDL_LONG nb,             /* number of templates */
											double *coeffs, /* coefficients */
											double *galaxy_flux, /* galaxy fluxes [i][k], 
																							redshifts */
											double *galaxy_invvar,
											double *galaxy_z,
											IDL_LONG ngalaxy)
{
	double *zgrid,*chi2,chi2min,az,bz,cz,sl;
	IDL_LONG i,j,k,b,p,nzsteps,jmin;

	/* allocate memory */
	pz_nk=nk;
	pz_nt=nt;
	pz_nb=nb;
	pz_nz=nz;
	pz_ematrix=(double *) malloc(pz_nb*pz_nt*sizeof(double));
	for(i=0;i<pz_nb*pz_nt;i++)
		pz_ematrix[i]=ematrix[i];
	pz_rmatrix=(double *) malloc(pz_nb*pz_nk*pz_nz*sizeof(double));
	for(i=0;i<pz_nb*pz_nk*pz_nz;i++)
		pz_rmatrix[i]=rmatrix[i];
	pz_zvals=(double *) malloc(pz_nz*sizeof(double));
	for(i=0;i<pz_nz;i++)
		pz_zvals[i]=zvals[i];
	pz_model_flux=(double *) malloc(pz_nk*sizeof(double));
	pz_coeffs=(double *) malloc(pz_nt*sizeof(double));
	pz_galaxy_flux=(double *) malloc(pz_nk*sizeof(double));
	pz_galaxy_invvar=(double *) malloc(pz_nk*sizeof(double));
	
	nzsteps=(IDL_LONG) floor((zvals[nz-1]-zvals[0])/ZRES);
	zgrid=(double *) malloc(nzsteps*sizeof(double));
	chi2=(double *) malloc(nzsteps*sizeof(double));
	for(j=0;j<nzsteps;j++) 
		zgrid[j]=zvals[0]+(zvals[nz-1]-zvals[0])*(double)j/(double)(nzsteps-1);
	for(i=0;i<ngalaxy;i++) {
		
		for(k=0;k<nk;k++) {
			pz_galaxy_flux[k]=galaxy_flux[i*nk+k];
			pz_galaxy_invvar[k]=galaxy_invvar[i*nk+k];
		} /* end for k */

		/* first lay out a grid and get chi^2 for each z value */
		for(j=0;j<nzsteps;j++) 
			chi2[j]=pz_fit_coeffs(zgrid[j]);

		/* find minimum on grid */
		jmin=0;
		chi2min=chi2[0];
		for(jmin=0,chi2min=chi2[0],j=1;j<nzsteps;j++) 
			if(chi2[j]<chi2min) {
				chi2min=chi2[j];
				jmin=j;
			} /* end if */

		/* then search for minimum more intelligently around it */
		if(jmin==0) {
			az=zgrid[jmin];
			bz=0.5*(zgrid[jmin]+zgrid[jmin+1]);
			cz=zgrid[jmin+1];
		} else if(jmin==nzsteps-1) {
			az=zgrid[jmin-1];
			bz=0.5*(zgrid[jmin-1]+zgrid[jmin]);
			cz=zgrid[jmin];
		} else {
			az=zgrid[jmin-1];
			bz=zgrid[jmin];
			cz=zgrid[jmin+1];
		} /* end if..else */
		chi2min=k_brent(az,bz,cz,pz_fit_coeffs,TOL,&(galaxy_z[i]));

		/* evaluate coefficients at final z */
		chi2min=pz_fit_coeffs(galaxy_z[i]);
		for(j=0;j<nt;j++)
			coeffs[i*nt+j]=pz_coeffs[j];
	} /* end for i */


	FREEVEC(pz_rmatrix);
	FREEVEC(pz_ematrix);
	FREEVEC(pz_zvals);
	FREEVEC(pz_galaxy_flux);
	FREEVEC(pz_galaxy_invvar);
	FREEVEC(pz_model_flux);
	FREEVEC(pz_coeffs);
	FREEVEC(zgrid);
	FREEVEC(chi2);
	return(1);
} /* end k_fit_photoz */

