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
IDL_LONG k_reconstruct_maggies(double *ematrix,    /* eigentemplates */
															IDL_LONG nt,     /* number of eigentemplates */
															double *zvals,      /* z interpolation */
															IDL_LONG nz,
															double *rmatrix,    /* r matrix */
															IDL_LONG nk,     /* number of bandpasses */
															IDL_LONG nb,     /* number of templates */
															double *coeffs, /* coefficients */
															double *galaxy_z,
															double *reconstruct_maggies,
															IDL_LONG ngalaxy)
{
	IDL_LONG i,j,k,b;
	
	for(i=0;i<ngalaxy;i++) {
		for(k=0;k<nk;k++) {
			reconstruct_maggies[k]=0.;
			for(b=0;b<nb;b++) {
				for(j=0;j<nt;j++) 
					reconstruct_maggies[k+i*nk]+=coeffs[i*nt+j]*ematrix[j*nb+b]*
						k_interpolate_es(galaxy_z[i],&(rmatrix[k*nb*nz+b*nz]),zvals,nz);
			} /* end for b */
		} /* end for k */
	} /* end for i */

	return(1);
} /* end k_fit_photoz */
