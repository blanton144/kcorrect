#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <kcorrect.h>

static double *war;
static IDL_LONG *iwar;

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}
static void free_memory()
{
	FREEVEC(war);
	FREEVEC(iwar);
}

/********************************************************************/
IDL_LONG idl_k_qld
  (int      argc,
   void *   argv[])
{
	IDL_LONG nconstraints,neconstraints,nconstraints_max;
	IDL_LONG nvar,nvar_max,iout,*ifail,*iprint;
	double *cmatrix, *dmatrix, *amatrix, *bmatrix;
	double *xl, *xu, *x, *lagrange;

	IDL_LONG i,mp2n,lwar,liwar;
	IDL_LONG retval=1;
	
	/* 0. allocate pointers from IDL */
	i=0;
	nconstraints=*((IDL_LONG *)argv[i]); i++;
	neconstraints=*((IDL_LONG *)argv[i]); i++;
	nconstraints_max=*((IDL_LONG *)argv[i]); i++;
	nvar=*((IDL_LONG *)argv[i]); i++;
	nvar_max=*((IDL_LONG *)argv[i]); i++;
	cmatrix=(double *)argv[i]; i++;
	dmatrix=(double *)argv[i]; i++;
	amatrix=(double *)argv[i]; i++;
	bmatrix=(double *)argv[i]; i++;
	xl=(double *)argv[i]; i++;
	xu=(double *)argv[i]; i++;
	x=(double *)argv[i]; i++;
	lagrange=(double *)argv[i]; i++;
	ifail=((IDL_LONG *)argv[i]); i++;
	iprint=((IDL_LONG *)argv[i]); i++;

	mp2n=nconstraints+nvar+nvar;
	lwar=3*nvar_max*nvar_max/2+10*nvar_max+2*nconstraints_max+10;
	war=(double *) malloc(sizeof(double)*lwar);
	liwar=2*nvar_max;
	iwar=(IDL_LONG *) malloc(sizeof(IDL_LONG)*liwar);
	iwar[0]=1;
	iout=6;

	/* 1. run the qld routine */
	ql0001_(&nconstraints,&neconstraints,&nconstraints_max, 
					&nvar,&nvar_max,&mp2n,cmatrix,dmatrix,amatrix, 
					bmatrix,xl,xu,x,lagrange,&iout,ifail,iprint, 
					war,&lwar,iwar,&liwar);
	
	/* 2. free memory and leave */
	free_memory();
	return retval;
}

/***************************************************************************/

