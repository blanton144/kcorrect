#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <kcorrect.h>

/*
 * k_load_filters.c
 *
 * Mike Blanton
 * 1/2002
 */

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}

/* load the filters in their standard format */
IDL_LONG k_load_filters(IDL_LONG **filter_n,
												double **filter_lambda,
												double **filter_pass,
												IDL_LONG *maxn,
												IDL_LONG *nk,
												char filterlist[])
{
	IDL_LONG i,j;
	FILE *listfp,*fp;
	char filename[2000],fullfilename[2000],filterdir[1000];

	FREEVEC((*filter_n));
	FREEVEC((*filter_lambda));
	FREEVEC((*filter_pass));
	
	if(getenv("KCORRECT_DIR")!=NULL) {
		sprintf(filterdir,"%s/data/filters",getenv("KCORRECT_DIR"));
	} else {
		strcpy(filterdir,".");
	} /* end if */

	listfp=k_fileopen(filterlist,"r");
	fscanf(listfp,"%d",nk);
	(*filter_n)=(IDL_LONG *) malloc((*nk)*sizeof(IDL_LONG));
	(*maxn)=0;
	for(i=0;i<(*nk);i++) {
		fscanf(listfp,"%s",filename);
		sprintf(fullfilename,"%s/%s",filterdir,filename);
		fp=k_fileopen(fullfilename,"r");
		fscanf(fp,"%d",&((*filter_n)[i]));
		fclose(fp);
		if((*filter_n)[i]>(*maxn)) (*maxn)=(*filter_n)[i];
	} /* end for i */
	fclose(listfp);
		
	(*filter_lambda)=(double *) malloc((*nk)*(*maxn)*sizeof(double));
	(*filter_pass)=(double *) malloc((*nk)*(*maxn)*sizeof(double));
	listfp=k_fileopen(filterlist,"r");
	fscanf(listfp,"%d",nk);
	for(i=0;i<(*nk);i++) {
		fscanf(listfp,"%s",filename);
		sprintf(fullfilename,"%s/%s",filterdir,filename);
		fp=k_fileopen(fullfilename,"r");
		fscanf(fp,"%d",&((*filter_n)[i]));
		for(j=0;j<(*filter_n)[i];j++)
			fscanf(fp,"%lf %lf",&((*filter_lambda)[i*(*maxn)+j]),
						 &((*filter_pass)[i*(*maxn)+j]));
		fclose(fp);
	} /* end for i */
	fclose(listfp);
	
	return(1);
} /* end k_load_ascii_table */
