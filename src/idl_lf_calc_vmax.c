#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "lf.h"

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}
static void free_memory()
{
}

/********************************************************************/
IDL_LONG idl_lf_calc_vmax(int      argc,
                          void *   argv[])
{
  IDL_LONG nv,nz,nk;
	float absm,*coeffs,*zvals,*rmatrix,sample_zmin,sample_zmax,mmin,mmax;
  float qevolve,qz0,band_shift,magoffset,omega0,omegal0,*zmin,*zmax;
  float absmagdep, ref_absmagdep;
	
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
	qevolve=*((float *)argv[i]); i++;
	qz0=*((float *)argv[i]); i++;
	absmagdep=*((float *)argv[i]); i++;
	ref_absmagdep=*((float *)argv[i]); i++;
	band_shift=*((float *)argv[i]); i++;
	magoffset=*((float *)argv[i]); i++;
	omega0=*((float *)argv[i]); i++;
	omegal0=*((float *)argv[i]); i++;
	zmin=((float *)argv[i]); i++;
	zmax=((float *)argv[i]); i++;
	
	/* 1. run the fitting routine */
	retval=lf_calc_vmax(absm,coeffs,nv,zvals,nz,rmatrix,nk,sample_zmin, 
                      sample_zmax,mmin,mmax,qevolve,qz0,absmagdep, 
                      ref_absmagdep, band_shift,magoffset,
                      omega0,omegal0,zmin,zmax);
	
	/* 2. free memory and leave */
	free_memory();
	return retval;
}

/***************************************************************************/

