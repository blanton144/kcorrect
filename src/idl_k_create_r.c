#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "kcorrect.h"

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}
static void free_memory()
{
}

/********************************************************************/
IDL_LONG idl_k_create_r
  (int      argc,
   void *   argv[])
{
	IDL_LONG nk,nb,nl,nz,*filter_n,maxn;
	double *rmatrix,*bmatrix,*lambda,*zvals,*filter_lambda,*filter_pass;
	
	IDL_LONG i,j;
	IDL_LONG retval=1;

	/* 0. allocate pointers from IDL */
	i=0;
	rmatrix=(double *)argv[i]; i++;
	nk=*((IDL_LONG *)argv[i]); i++;
	nb=*((IDL_LONG *)argv[i]); i++;
	bmatrix=(double *)argv[i]; i++;
	lambda=(double *)argv[i]; i++;
	nl=*((IDL_LONG *)argv[i]); i++;
	zvals=(double *)argv[i]; i++;
	nz=*((IDL_LONG *)argv[i]); i++;
	filter_n=(IDL_LONG *)argv[i]; i++;
	filter_lambda=(double *)argv[i]; i++;
	filter_pass=(double *)argv[i]; i++;
	maxn=*((IDL_LONG *)argv[i]); i++;
	
	/* 1. run the fitting routine */
	retval=(IDL_LONG) k_create_r(rmatrix,nk,nb,bmatrix,lambda,nl,zvals,nz,
															 filter_n,filter_lambda,filter_pass,maxn);
	
	/* 2. free memory and leave */
	free_memory();
	return retval;
}

/***************************************************************************/

