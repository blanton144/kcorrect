#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <kcorrect.h>

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}
static void free_memory()
{
}

/* IDL/C wrapper on code to create the rmatrix, which expresses how
   each template contributes to each filter at each redshift */

/********************************************************************/
IDL_LONG idl_k_projection_table
  (int      argc,
   void *   argv[])
{
	IDL_LONG nk,nb,nl,nz,*filter_n,maxn;
	float *rmatrix,*bmatrix,*lambda,*zvals,*filter_lambda,*filter_pass;
	float band_shift;
	
	IDL_LONG i;
	IDL_LONG retval=1;

	/* 0. allocate pointers from IDL */
	i=0;
	rmatrix=(float *)argv[i]; i++;
	nk=*((IDL_LONG *)argv[i]); i++;
	nb=*((IDL_LONG *)argv[i]); i++;
	bmatrix=(float *)argv[i]; i++;
	lambda=(float *)argv[i]; i++;
	nl=*((IDL_LONG *)argv[i]); i++;
	zvals=(float *)argv[i]; i++;
	nz=*((IDL_LONG *)argv[i]); i++;
	filter_n=(IDL_LONG *)argv[i]; i++;
	filter_lambda=(float *)argv[i]; i++;
	filter_pass=(float *)argv[i]; i++;
	band_shift=*((float *)argv[i]); i++;
	maxn=*((IDL_LONG *)argv[i]); i++;
	
	/* 1. run the fitting routine */
	retval=(IDL_LONG) k_projection_table(rmatrix,nk,nb,bmatrix,lambda,nl,zvals,
																			 nz,filter_n,filter_lambda,filter_pass, 
																			 band_shift,maxn);
	
	/* 2. free memory and leave */
	free_memory();
	return retval;
}

/***************************************************************************/

