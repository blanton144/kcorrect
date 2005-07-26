#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <kcorrect.h>

/*
 * k_fit_photoz.c
 *
 * Given maggies (and maybe a prior) pick the best redshift 
 *
 * Mike Blanton
 * 5/2003 */

#define ZRES 0.1
#define ZTOL 0.02
#define TOL 1.e-4
#define MAXITERFULL 10000
#define MAXITER 2000
#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}

static float *pz_coeffs=NULL;
static float *pz_save_coeffs=NULL;
static float *pz_rmatrix=NULL;
static float *pz_zvals=NULL;
static float *pz_lprior=NULL;
static float *pz_zprior=NULL;
static float *pz_maggies=NULL;
static float *pz_maggies_ivar=NULL;
static IDL_LONG pz_nk=0,pz_nv=0,pz_nz=0,pz_nprior,pz_ncheck;
static IDL_LONG pz_init=1, pz_applyprior=1;

float pz_fit_coeffs(float z) 
{
	IDL_LONG niter,dontinit,maxiter;
	float chi2;
	
  if(pz_init==0) {
    dontinit=0;
		maxiter=MAXITERFULL;
  } else {
    dontinit=1;
		maxiter=MAXITER;
	}
  k_fit_nonneg(pz_coeffs,pz_rmatrix,pz_nk,pz_nv,pz_zvals,pz_nz, 
               pz_maggies,pz_maggies_ivar,&z,1,TOL,maxiter,&niter,
               &chi2,0,dontinit);
	if(pz_applyprior) chi2-=k_interpolate_es(z, pz_lprior, pz_zprior, pz_nprior);
  pz_ncheck++;
	
	return(chi2);
} /* end pz_fit_coeffs */

IDL_LONG k_fit_photoz(float *photoz,
											float *coeffs,
											float *rmatrix,
											IDL_LONG nk,
											IDL_LONG nv,
                      float *zprior, 
											float *lprior,
                      IDL_LONG nprior,
											float *zvals,
											IDL_LONG nz,
											float *maggies,
											float *maggies_ivar,
											IDL_LONG ngalaxy,
											float tolerance,
											IDL_LONG maxiter,
											IDL_LONG *niter,
											float *chi2,
											IDL_LONG verbose)
{
	float *lzgrid, *lchi2grid, *zgrid,*chi2grid,chi2min,az,bz,cz;
  float z_0,z_1,z_2,chi2_0,chi2_1,chi2_2;
	IDL_LONG i,j,k,nzsteps,lnzsteps,jmin;

  pz_applyprior=1;
	pz_nk=nk;
	pz_nv=nv;
	pz_nz=nz;
	pz_nprior=nprior;
	nzsteps=(IDL_LONG) floor((zvals[nz-1]-zvals[0])/ZRES)+1;
	chi2grid=(float *) malloc(nzsteps*sizeof(float));
	pz_maggies=(float *) malloc(nk*sizeof(float));
	pz_maggies_ivar=(float *) malloc(nk*sizeof(float));
	pz_coeffs=(float *) malloc(nv*sizeof(float));
	pz_save_coeffs=(float *) malloc(nv*sizeof(float));
	pz_rmatrix=(float *) malloc(pz_nv*pz_nk*pz_nz*sizeof(float));
	for(i=0;i<pz_nv*pz_nk*pz_nz;i++)
		pz_rmatrix[i]=rmatrix[i];
	pz_zvals=(float *) malloc(pz_nz*sizeof(float));
	for(i=0;i<pz_nz;i++)
		pz_zvals[i]=zvals[i];
	pz_lprior=(float *) malloc(pz_nprior*sizeof(float));
	for(i=0;i<pz_nprior;i++)
		pz_lprior[i]=lprior[i];
	pz_zprior=(float *) malloc(pz_nprior*sizeof(float));
	for(i=0;i<pz_nprior;i++)
		pz_zprior[i]=zprior[i];
	zgrid=(float *) malloc(nzsteps*sizeof(float));
	for(j=0;j<nzsteps;j++) 
		zgrid[j]=zvals[0]+(zvals[nz-1]-zvals[0])*(float)j/(float)(nzsteps-1);
  lchi2grid=(float *) malloc(10000*sizeof(float));
  lzgrid=(float *) malloc(10000*sizeof(float));
	for(i=0;i<ngalaxy;i++) {

    pz_ncheck=0;
		
		for(k=0;k<nk;k++) {
			pz_maggies[k]=maggies[i*nk+k];
			pz_maggies_ivar[k]=maggies_ivar[i*nk+k];
		} /* end for k */

		/* first lay out a grid and get chi^2 for each z value */
		pz_init=1;
		for(j=0;j<nzsteps;j++) 
			chi2grid[j]=pz_fit_coeffs(zgrid[j]);

		/* find minimum on grid */
		for(jmin=0,chi2min=chi2grid[0],j=1;j<nzsteps;j++) 
			if(chi2grid[j]<chi2min) {
				chi2min=chi2grid[j];
				jmin=j;
			} /* end if */

		chi2min=pz_fit_coeffs(zgrid[jmin]);
		for(j=0;j<nv;j++)
			pz_save_coeffs[j]=pz_coeffs[j];
    
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

		/* lay out a grid and get chi^2 for each z value */
    lnzsteps=(IDL_LONG) floor((cz-az)/ZTOL)+1;
    for(j=0;j<lnzsteps;j++) 
      lzgrid[j]=az+(cz-az)*(float)j/(float)(lnzsteps-1);
		pz_init=0;
		for(j=0;j<lnzsteps;j++) {
			for(k=0;k<nv;k++)
				pz_coeffs[k]=pz_save_coeffs[k];
			lchi2grid[j]=pz_fit_coeffs(lzgrid[j]);
		}

		/* find minimum on grid */
		for(jmin=0,chi2min=lchi2grid[0],j=1;j<lnzsteps;j++) 
			if(lchi2grid[j]<chi2min) {
				chi2min=lchi2grid[j];
				jmin=j;
			} /* end if */

    photoz[i]=lzgrid[jmin];
		
    /* near bottom fit a parabola and get min */
    if(photoz[i]>ZTOL) {
      z_0=photoz[i]-ZTOL;
      z_1=photoz[i];
      z_2=photoz[i]+ZTOL;
    } else {
      z_0=1.e-4;
      z_1=ZTOL;
      z_2=2.*ZTOL;
    }
    chi2_0=pz_fit_coeffs(z_0);
    chi2_1=pz_fit_coeffs(z_1);
    chi2_2=pz_fit_coeffs(z_2);
    photoz[i]=z_0+0.5*ZTOL-(chi2_1-chi2_0)*ZTOL/(chi2_2-2.*chi2_1+chi2_0);
    if(photoz[i]<z_0) photoz[i]=z_0;
    if(photoz[i]>z_2) photoz[i]=z_2;

		/* evaluate coefficients at final z (using REAL chi2) */
    pz_applyprior=0;
		chi2[i]=pz_fit_coeffs(photoz[i]);
    pz_applyprior=1;
		for(j=0;j<nv;j++)
			coeffs[i*nv+j]=pz_coeffs[j];

    if(verbose) {
      printf("nredshifts=%d\n", pz_ncheck);
    } 
	} /* end for i */

	FREEVEC(pz_rmatrix);
	FREEVEC(pz_zvals);
	FREEVEC(pz_maggies);
	FREEVEC(pz_maggies_ivar);
	FREEVEC(pz_coeffs);
	FREEVEC(pz_lprior);
	FREEVEC(pz_zprior);
	FREEVEC(zgrid);
	FREEVEC(chi2grid);
	FREEVEC(lzgrid);
	FREEVEC(lchi2grid);
	return(1);
} /* end k_fit_photoz */
