#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <kcorrect.h>

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}
static void free_memory()
{
}

/* IDL/C wrapper on code to reconstruct maggies from a set of
   templates and coefficients */

/********************************************************************/
IDL_LONG idl_k_reconstruct_maggies
  (int      argc,
   void *   argv[])
{
   IDL_LONG nz,nk,nv,ngalaxy;
   float *zvals,*rmatrix,*coeffs,*reconstruct_maggies,*galaxy_z;

	 IDL_LONG i;
	 IDL_LONG retval=1;

   /* 0. allocate pointers from IDL */
	 i=0;
	 zvals=(float *)argv[i]; i++;
   nz=*((IDL_LONG *)argv[i]); i++;
	 rmatrix=(float *)argv[i]; i++;
   nk=*((IDL_LONG *)argv[i]); i++;
   nv=*((IDL_LONG *)argv[i]); i++;
	 coeffs=(float *)argv[i]; i++;
	 galaxy_z=(float *)argv[i]; i++;
	 reconstruct_maggies=(float *)argv[i]; i++;
   ngalaxy=*((IDL_LONG *)argv[i]); i++;

	 /* 1. run the fitting routine */
	 retval=(IDL_LONG) k_reconstruct_maggies(zvals,nz,rmatrix,nk,nv,coeffs,
																					 galaxy_z,reconstruct_maggies,
																					 ngalaxy);

	 /* 2. free memory and leave */
	 free_memory();
   return retval;
}

/***************************************************************************/

