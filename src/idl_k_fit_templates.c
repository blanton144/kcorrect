#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <kcorrect.h>

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}
static void free_memory()
{
}

/********************************************************************/
IDL_LONG idl_k_fit_templates
  (int      argc,
   void *   argv[])
{
   IDL_LONG nt,nz,nk,nb,ngalaxy,*galaxy_clip,initialized_dmatrix;
   double *ematrix,*zvals,*rmatrix,*amatrix,*dmatrix;
	 double *galaxy_maggies,*galaxy_invvar,*galaxy_z;

	 IDL_LONG i;
	 IDL_LONG retval=1;

   /* 0. allocate pointers from IDL */
	 i=0;
	 ematrix=(double *)argv[i]; i++;
   nt=*((IDL_LONG *)argv[i]); i++;
	 zvals=(double *)argv[i]; i++;
   nz=*((IDL_LONG *)argv[i]); i++;
	 rmatrix=(double *)argv[i]; i++;
   nk=*((IDL_LONG *)argv[i]); i++;
   nb=*((IDL_LONG *)argv[i]); i++;
	 amatrix=(double *)argv[i]; i++;
	 galaxy_maggies=(double *)argv[i]; i++;
	 galaxy_invvar=(double *)argv[i]; i++;
	 galaxy_z=(double *)argv[i]; i++;
	 galaxy_clip=(IDL_LONG *)argv[i]; i++;
   ngalaxy=*((IDL_LONG *)argv[i]); i++;
	 dmatrix=(double *)argv[i]; i++;
   initialized_dmatrix=*((IDL_LONG *)argv[i]); i++;

	 /* 1. run the fitting routine */
	 retval=(IDL_LONG) k_fit_templates(ematrix,nt,zvals,nz,rmatrix,nk,nb,
																		 amatrix,galaxy_maggies,galaxy_invvar,
																		 galaxy_z,galaxy_clip,ngalaxy,dmatrix,
																		 initialized_dmatrix);

	 /* 2. free memory and leave */
	 free_memory();
   return retval;
}

/***************************************************************************/

