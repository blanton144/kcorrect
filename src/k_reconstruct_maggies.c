#include <stdio.h>
#include <string.h>
#include <math.h>
#include <kcorrect.h>

/*
 * k_reconstruct_maggies.c
 *
 * Returns reconstructed fluxes, given the model
 *
 * Mike Blanton
 * 1/2002
 */

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}

/* calculate the reconstructed fluxes, given the coeffs of the model */
IDL_LONG k_reconstruct_maggies(float *zvals,      /* z interpolation */
															 IDL_LONG nz,
															 float *rmatrix,    /* r matrix */
															 IDL_LONG nk,     /* number of bandpasses */
															 IDL_LONG nv,     /* number of templates */
															 float *coeffs, /* coefficients */
															 float *redshift,
															 float *reconstruct_maggies,
															 IDL_LONG ngalaxy)
{
	IDL_LONG i,j,k;
	
	for(i=0;i<ngalaxy;i++) {
		for(k=0;k<nk;k++) {
			reconstruct_maggies[k+i*nk]=0.;
			for(j=0;j<nv;j++) 
				reconstruct_maggies[k+i*nk]+=coeffs[i*nv+j]*
					k_interpolate_es(redshift[i],&(rmatrix[k*nv*nz+j*nz]),zvals,nz);
		} /* end for k */
	} /* end for i */
	
	return(1);
} /* end k_fit_photoz */
