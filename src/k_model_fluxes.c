#include <stdio.h>
#include <string.h>
#include <math.h>
#include "kcorrect.h"

/*
 * k_model_fluxes.c
 *
 * Returns model fluxes, given the model
 *
 * Mike Blanton
 * 1/2002
 */

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}

/* calculate the model fluxes, given the coeffs of the model */
IDL_LONG k_model_fluxes(double *ematrix,    /* eigentemplates */
												IDL_LONG nt,             /* number of eigentemplates */
												double *zvals,      /* z interpolation */
												IDL_LONG nz,
												double *rmatrix,    /* r matrix */
												IDL_LONG nk,             /* number of bandpasses */
												IDL_LONG nb,             /* number of templates */
												double *coeffs, /* coefficients */
												double *galaxy_z,
												double *model_flux,
												IDL_LONG ngalaxy)
{
	IDL_LONG i,j,k,b;
	
	for(i=0;i<ngalaxy;i++) {
		for(k=0;k<nk;k++) {
			model_flux[k]=0.;
			for(b=0;b<nb;b++) {
				for(j=0;j<nt;j++) 
					model_flux[k+i*nk]+=coeffs[i*nt+j]*ematrix[j*nb+b]*
						k_interpolate(galaxy_z[i],&(rmatrix[k*nb*nz+b*nz]),zvals,nz);
			} /* end for b */
		} /* end for k */
	} /* end for i */

	return(1);
} /* end k_fit_photoz */
