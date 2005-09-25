#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <kcorrect.h>

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}
static void free_memory()
{
}

/* IDL/C wrapper on code to fit nonnegative coefficients to a set of
   maggies */

/********************************************************************/
IDL_LONG idl_k_binspec
  (int      argc,
   void *   argv[])
{
	IDL_LONG nl,nnewl;
	float *spec,*newspec,*lambda,*newlambda;
	
	IDL_LONG i;
	IDL_LONG retval=1;

	/* 0. allocate pointers from IDL */
	i=0;
	lambda=(float *)argv[i]; i++;
	spec=(float *)argv[i]; i++;
	newlambda=(float *)argv[i]; i++;
	newspec=(float *)argv[i]; i++;
	nl=*((IDL_LONG *)argv[i]); i++;
	nnewl=*((IDL_LONG *)argv[i]); i++;
	
	/* 1. run the fitting routine */
	retval=(IDL_LONG) k_binspec(lambda, spec, newlambda, newspec, nl, nnewl);
															
	/* 2. free memory and leave */
	free_memory();
	return retval;
}

/***************************************************************************/

