#include <stdio.h>
#include <string.h>
#include <math.h>
#include "kcorrect.h"

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
static double *pz_ematrix1=NULL;
static double *pz_coeffspnorm=NULL;
static double *pz_coeffs1=NULL;
static double pz_z,pz_linepos;
static IDL_LONG pz_nt,pz_nz,pz_nk,pz_nb,pz_puse,pz_np,pz_p;

#if 0
static double pz_coeffsps[12]={
	0.07,0.08,-0.19,
	0.011,0.085,-0.17,
	0.009,0.09,-0.15,
	0.03,0.09,-0.1
};
static double pz_coeffspe[12]={
	0.011,0.085,-0.17,
	0.009,0.09,-0.15,
	0.03,0.09,-0.1,
	0.12,0.08,0.14
};
#endif

static double pz_coeffsps[12]={
	0.04,0.082,-0.185,
	0.04,0.082,-0.175,
	0.04,0.082,-0.15,
	0.04,0.082,-0.1
};
static double pz_coeffspe[12]={
	0.04,0.082,-0.175,
	0.04,0.082,-0.15,
	0.04,0.082,-0.1,
	0.04,0.082,0.0
};

double pz_calc_chi2(double z) 
{
	double *model_flux,chi2;
	IDL_LONG k;
	
	k_model_fluxes(pz_ematrix1,1,pz_zvals,pz_nz,pz_rmatrix,
								 pz_nk,pz_nb,pz_coeffs1,&z,pz_model_flux,1);
	chi2=0.;
	for(k=0;k<pz_nk;k++)
		chi2+=(pz_model_flux[k]-pz_galaxy_flux[k])
			*(pz_model_flux[k]-pz_galaxy_flux[k])*pz_galaxy_invvar[k];
	
	return(chi2);
} /* end pz_calc_chi2 */

double pz_line_fit(double linepos) 
{
	double sl,chi2;
	IDL_LONG b,i,j;

	for(i=0;i<pz_np;i++) {
#if 0
		fprintf(stderr,"%e %e\n",pz_coeffspnorm[i+1]/pz_coeffspnorm[pz_np],linepos);
#endif
		if(pz_coeffspnorm[i]/pz_coeffspnorm[pz_np]<linepos &&
			 pz_coeffspnorm[i+1]/pz_coeffspnorm[pz_np]>linepos)
			break;
	}
	if(i==pz_np) i=pz_np-1;
	sl=(linepos*pz_coeffspnorm[pz_np]-pz_coeffspnorm[i])
		/(pz_coeffspnorm[i+1]-pz_coeffspnorm[i]);
	pz_coeffs[0]=1.;
	for(j=1;j<pz_nt;j++) 
		pz_coeffs[j]=pz_coeffsps[i*(pz_nt-1)+(j-1)]
			+(pz_coeffspe[i*(pz_nt-1)+(j-1)]-pz_coeffsps[i*(pz_nt-1)+(j-1)])*sl;

	for(b=0;b<pz_nb;b++) {
		pz_ematrix1[b]=0.;
		for(j=0;j<pz_nt;j++) 
			pz_ematrix1[b]+=pz_coeffs[j]*pz_ematrix[j*pz_nb+b];
	} /* end for b */
	
	k_fit_coeffs(pz_ematrix1,1,pz_zvals,pz_nz,pz_rmatrix,
							 pz_nk,pz_nb,pz_coeffs1,pz_galaxy_flux,pz_galaxy_invvar,&pz_z,1);
	chi2=pz_calc_chi2(pz_z);

#if 0
	fprintf(stderr,"%e %e %d %e %e\n",pz_z,linepos,i,sl,chi2);
#endif
	
	return(chi2);
} /* end pz_line_fit */

/* calculates chi2 */
double pz_fit_coeffs(double z) 
{
	double chi2min,az,bz,cz;
	IDL_LONG i,k;

	/* find the best along the chain */
	az=0.;
	bz=0.5;
	cz=1.;
	pz_z=z;
	chi2min=k_brent(az,bz,cz,pz_line_fit,TOL,&pz_linepos);
#if 0
	fprintf(stderr,"linepos=%e\n",pz_linepos); fflush(stderr);
#endif

	return(chi2min);

#if 0
	/* try the best of several piece-wise fits */
	chi2min=-1.;
	for(pz_p=0;pz_p<pz_np;pz_p++) {

		/* fit coeffs along this line */
		k_fit_coeffs(&(pz_ematrix1[pz_p*pz_nb*2]),2,pz_zvals,pz_nz,pz_rmatrix,
								 pz_nk,pz_nb,pz_coeffs1,pz_galaxy_flux,pz_galaxy_invvar,&z,1);
								 
		if(pz_coeffs1[1]/pz_coeffs1[0]>0. && pz_coeffs1[1]/pz_coeffs1[0]<1.) {
			/* calculate chi2 ... */
			chi2=pz_calc_chi2(z);
		} else {
			/* unless the best fit is off the edge of the line. in this case
				 pick the best edge */ 
#if 0
			fprintf(stderr,"yuckky\n");
#endif
			pz_coeffs1[1]=0.;
			chi2lo=pz_calc_chi2(z);
			pz_coeffs1[1]=pz_coeffs1[0];
			chi2hi=pz_calc_chi2(z);
			chi2=(chi2lo<chi2hi) ? chi2lo : chi2hi;
		} /* end if..else */

		/* is this the best fit so far? if so, keep it */
		if(chi2min<0. || chi2<chi2min) {
			chi2min=chi2;
			pz_puse=pz_p;
		} /* end if */
#if 0
		fprintf(stderr,"puse=%d p=%d chi2=%e %e %e\n",pz_puse,pz_p,chi2,
						pz_coeffs1[0],pz_coeffs1[1]/pz_coeffs1[0]);
		fflush(stdout);
#endif

	} /* end for p */
#endif
		
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
	IDL_LONG np=1;

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
	pz_np=np;
	pz_coeffspnorm=(double *) malloc((pz_np+1)*sizeof(double));
	pz_ematrix1=(double *) malloc(pz_nb*sizeof(double));
	pz_coeffs1=(double *) malloc(1*sizeof(double));
	
	/* set up piece-wise stuff */
#if 0
	fprintf(stderr,"pz_np=%d\n",pz_np); fflush(stdout);
#endif
	pz_coeffspnorm[0]=0.;
	for(p=0;p<pz_np;p++) {
		pz_coeffspnorm[p+1]=pz_coeffspnorm[p]*pz_coeffspnorm[p];
		for(j=1;j<pz_nt;j++) {
			pz_coeffspnorm[p+1]+=((pz_coeffspe[(pz_nt-1)*p+j-1]
														 -pz_coeffsps[(pz_nt-1)*p+j-1])*
														(pz_coeffspe[(pz_nt-1)*p+j-1]
														 -pz_coeffsps[(pz_nt-1)*p+j-1]));
		} /* end for j */
		pz_coeffspnorm[p+1]=sqrt(pz_coeffspnorm[p+1]);
#if 0
		fprintf(stderr,"pz_coeffspnorm=%e\n",pz_coeffspnorm[p+1]);
#endif
#if 0
		for(b=0;b<pz_nb;b++) {
			pz_ematrix1[pz_nb*2*p+0*pz_nb+b]=pz_ematrix[0*pz_nb+b];
			pz_ematrix1[pz_nb*2*p+1*pz_nb+b]=0.;
			for(j=1;j<pz_nt;j++) {
				pz_ematrix1[pz_nb*2*p+0*pz_nb+b]+=pz_ematrix[pz_nb*j+b]
					*pz_coeffsps[(pz_nt-1)*p+j-1];
				pz_ematrix1[pz_nb*2*p+1*pz_nb+b]+=
					pz_ematrix[pz_nb*j+b]*
					(pz_coeffspe[(pz_nt-1)*p+j-1]-pz_coeffsps[(pz_nt-1)*p+j-1]);
			} /* end for j */
		} /* end for b */
#endif
	} /* end for p */

	nzsteps=(int) floor((zvals[nz-1]-zvals[0])/ZRES);
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

		/* calculate coeffs for this redshift */
		pz_z=galaxy_z[i];
		chi2min=k_brent(0.,0.5,1.,pz_line_fit,TOL,&pz_linepos);

		for(p=0;p<pz_np;p++) 
			if(pz_coeffspnorm[p]/pz_coeffspnorm[pz_np]<pz_linepos &&
				 pz_coeffspnorm[p+1]/pz_coeffspnorm[pz_np]>pz_linepos)
				break;
		if(p==pz_np) p=pz_np-1;
		sl=(pz_linepos*pz_coeffspnorm[pz_np]-pz_coeffspnorm[p])
			/(pz_coeffspnorm[p+1]-pz_coeffspnorm[p]);
		coeffs[0]=pz_coeffs1[0];
		for(j=1;j<pz_nt;j++) 
			coeffs[j]=coeffs[0]*(pz_coeffsps[p*(pz_nt-1)+(j-1)]
													 +(pz_coeffspe[p*(pz_nt-1)+(j-1)]
														 -pz_coeffsps[p*(pz_nt-1)+(j-1)])*sl);

#if 0
		k_fit_coeffs(pz_ematrix,pz_nt,pz_zvals,pz_nz,pz_rmatrix,pz_nk,pz_nb,
								 &(coeffs[i*nt]),pz_galaxy_flux,pz_galaxy_invvar,
								 &(galaxy_z[i]),1);

#endif

#if 0
		pz_fit_coeffs(galaxy_z[i]);
#endif
		
	} /* end for i */


	FREEVEC(pz_rmatrix);
	FREEVEC(pz_ematrix);
	FREEVEC(pz_zvals);
	FREEVEC(pz_galaxy_flux);
	FREEVEC(pz_galaxy_invvar);
	FREEVEC(pz_model_flux);
	FREEVEC(pz_coeffs);
	FREEVEC(pz_coeffs1);
	FREEVEC(pz_coeffspnorm);
	FREEVEC(pz_ematrix1);
	FREEVEC(zgrid);
	FREEVEC(chi2);
	return(1);
} /* end k_fit_photoz */

