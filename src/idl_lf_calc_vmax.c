#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "kcorrect.h"

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}
static void free_memory()
{
}

/* IDL/C wrapper on code to calculate the maximum volume over which
   you can observe an object, given its absolute magnitude, its
   K-correction coefficients, the K-correction info (zvals, rmatrix),
   the sample limits, the flux limits, and assumed evolution
   parameters */

/********************************************************************/
IDL_LONG idl_lf_calc_vmax(int      argc,
                          void *   argv[])
{
  IDL_LONG nv,nz,nk;
	float absm,*coeffs,*zvals,*rmatrix,sample_zmin,sample_zmax,mmin,mmax;
  float q0,q1,qz0,band_shift,magoffset,omega0,omegal0,*zmin,*zmax;
	
	IDL_LONG i;
	IDL_LONG retval=1;

	/* 0. allocate pointers from IDL */
	i=0;
	absm=*((float *)argv[i]); i++;
	coeffs=((float *)argv[i]); i++;
	nv=*((int *)argv[i]); i++;
	zvals=((float *)argv[i]); i++;
	nz=*((int *)argv[i]); i++;
	rmatrix=((float *)argv[i]); i++;
	nk=*((int *)argv[i]); i++;
	sample_zmin=*((float *)argv[i]); i++;
	sample_zmax=*((float *)argv[i]); i++;
	mmin=*((float *)argv[i]); i++;
	mmax=*((float *)argv[i]); i++;
	q0=*((float *)argv[i]); i++;
	q1=*((float *)argv[i]); i++;
	qz0=*((float *)argv[i]); i++;
	band_shift=*((float *)argv[i]); i++;
	magoffset=*((float *)argv[i]); i++;
	omega0=*((float *)argv[i]); i++;
	omegal0=*((float *)argv[i]); i++;
	zmin=((float *)argv[i]); i++;
	zmax=((float *)argv[i]); i++;
	
	/* 1. run the fitting routine */
	retval=lf_calc_vmax(absm,coeffs,nv,zvals,nz,rmatrix,nk,sample_zmin, 
                      sample_zmax,mmin,mmax,q0,q1,qz0,band_shift,magoffset,
                      omega0,omegal0,zmin,zmax);
	
	/* 2. free memory and leave */
	free_memory();
	return retval;
}

/***************************************************************************/

