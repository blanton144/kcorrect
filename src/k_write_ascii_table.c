#include <stdio.h>
#include <string.h>
#include <math.h>
#include "kcorrect.h"

/*
 * k_write_ascii_table.c
 *
 * Write an ascii table, of the form:
 *
 * <ndim> <size_{0}> ... <size_{ndim-1}>
 * <entry_0>
 * <entry_1>
 * ...
 * <entry_n>
 * ...
 * <entry_{size_0*size_1*..*size_{ndim-1}-1>
 *
 * where the table element [i][j][k] (for ndim==3) would be the 
 * entry element n, where n=i*size_2*size1+j*size2+k
 *
 * Mike Blanton
 * 1/2002
 */

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}

/* write an ascii table in my standard format */
IDL_LONG k_write_ascii_table(double *table,
														 IDL_LONG ndim,
														 IDL_LONG *sizes,
														 char filename[])
{
	IDL_LONG i,nelem;
	FILE *fp;

	fp=k_fileopen(filename,"w");
	fprintf(fp,"%d ",ndim);
	for(i=0;i<ndim;i++) 
		fprintf(fp,"%d ",sizes[i]);
	fprintf(fp,"\n");
	nelem=1;
	for(i=0;i<ndim;i++) 
		nelem*=sizes[i];
	for(i=0;i<nelem;i++) 
		fprintf(fp,"%e\n",table[i]);
	fclose(fp);

} /* end k_write_ascii_table */
