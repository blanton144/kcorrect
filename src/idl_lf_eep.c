#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "kcorrect.h"

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}
static void free_memory()
{
}

/* IDL/C wrapper on code to calculate the EEP luminosity function */

/********************************************************************/
IDL_LONG idl_lf_eep(int      argc,
                    void *   argv[])
{
  IDL_LONG ngals, nbin, calc_err;
  float *redshift, *absmag, *absmmin, *absmmax, *absmk, *phi, *phi_err, *covar;
  float *weight;
  float sample_absmmin,sample_absmmax;
	
	IDL_LONG i;
	IDL_LONG retval=1;

	/* 0. allocate pointers from IDL */
	i=0;
	redshift=((float *)argv[i]); i++;
	absmag=((float *)argv[i]); i++;
	absmmin=((float *)argv[i]); i++;
	absmmax=((float *)argv[i]); i++;
	ngals=*((IDL_LONG *)argv[i]); i++;
	sample_absmmin=*((float *)argv[i]); i++;
	sample_absmmax=*((float *)argv[i]); i++;
	absmk=((float *)argv[i]); i++;
	phi=((float *)argv[i]); i++;
	phi_err=((float *)argv[i]); i++;
  covar=((float *)argv[i]); i++;
	nbin=*((IDL_LONG *)argv[i]); i++;
	calc_err=*((IDL_LONG *)argv[i]); i++;
	weight=((float *)argv[i]); i++;
	
	/* 1. run the fitting routine */
	retval=lf_eep(redshift,absmag,absmmin,absmmax,ngals,
                sample_absmmin,sample_absmmax,absmk,
                phi,phi_err,covar,nbin,calc_err,weight);
	
	/* 2. free memory and leave */
	free_memory();
	return retval;
}

/***************************************************************************/

