#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "kcorrect.h"

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}
static void free_memory()
{
}

/* IDL/C wrapper on code to calculate the selection function based on
   an EEP luminosity function */

/********************************************************************/
IDL_LONG idl_lf_select_eep(int      argc,
                           void *   argv[])
{
  IDL_LONG ngals, nbin;
  float *redshift, *absmag, *absmmin, *absmmax, *absmk, *phi, *sel;
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
	sel=((float *)argv[i]); i++;
	nbin=*((IDL_LONG *)argv[i]); i++;
	
	/* 1. run the fitting routine */
	retval=lf_select_eep(redshift,absmag,absmmin,absmmax,ngals,
                       sample_absmmin,sample_absmmax,
                       absmk,phi,sel,nbin);
	
	/* 2. free memory and leave */
	free_memory();
	return retval;
}

/***************************************************************************/

