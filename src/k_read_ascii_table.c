#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <kcorrect.h>

/*
 * k_read_ascii_table.c
 *
 * Load an ascii table, of the form:
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

/* load an ascii table in my standard format */
IDL_LONG k_read_ascii_table(float **table,
														IDL_LONG *ndim,
														IDL_LONG **sizes,
														char filename[])
{
	IDL_LONG i,nelem;
	FILE *fp;

	FREEVEC((*table));
	FREEVEC((*sizes));

	fp=k_fileopen(filename,"r");
	fscanf(fp,"%d",ndim);
	if((*ndim)<=0) {
		fprintf(stderr,"k_read_ascii_table found 0-dimensional table in %s\n",
						filename);
		exit(1);
	} /* end if */
	nelem=1;
	(*sizes)=(IDL_LONG *) malloc((*ndim)*sizeof(IDL_LONG));
	for(i=0;i<(*ndim);i++) {
		fscanf(fp,"%d",&((*sizes)[i]));
		nelem*=(*sizes)[i];
	} /* end for i */
	if(nelem<=0) {
		fprintf(stderr,"k_read_ascii_table found 0-size table in %s\n",
						filename);
		exit(1);
	} /* end if */
	(*table)=(float *) malloc(nelem*sizeof(float));
	for(i=0;i<nelem;i++) 
		fscanf(fp,"%f",&((*table)[i]));
	fclose(fp);

	return(1);
} /* end k_read_ascii_table */
