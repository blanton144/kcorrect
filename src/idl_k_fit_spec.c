#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <kcorrect.h>

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}
static void free_memory()
{
}

/* IDL/C wrapper on code to fit nonnegative coefficients to a set of
   fluxes in a spectrum */

/********************************************************************/
IDL_LONG idl_k_fit_spec
  (int      argc,
   void *   argv[])
{
	IDL_LONG nl,nt,maxiter,*niter,verbose;
	float *coeffs,*templates,*flux,*ivar,*chi2;
	float tolerance;
	
	IDL_LONG i;
	IDL_LONG retval=1;

	/* 0. allocate pointers from IDL */
	i=0;
	coeffs=(float *)argv[i]; i++;
	flux=(float *)argv[i]; i++;
	ivar=(float *)argv[i]; i++;
	templates=(float *)argv[i]; i++;
	nt=*((IDL_LONG *)argv[i]); i++;
	nl=*((IDL_LONG *)argv[i]); i++;
	tolerance=*((float *)argv[i]); i++;
	maxiter=*((IDL_LONG *)argv[i]); i++;
	niter=(IDL_LONG *)argv[i]; i++;
	chi2=(float *)argv[i]; i++;
	verbose=*((IDL_LONG *)argv[i]); i++;
	
	/* 1. run the fitting routine */
	retval=(IDL_LONG) k_fit_spec(coeffs,flux,ivar,templates,nt,nl,tolerance,
                               maxiter,niter,chi2,verbose);
	
	/* 2. free memory and leave */
	free_memory();
	return retval;
}
/***************************************************************************/

