#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <kcorrect.h>

/*
 * k_load_filters.c
 *
 * Load filters into a filter structure. Used by the stand-alone C
 * code.
 *
 * Mike Blanton
 * 1/2002
 */

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}

/* load the filters in their standard format */
IDL_LONG k_load_filters(IDL_LONG **filter_n,
												float **filter_lambda,
												float **filter_pass,
												IDL_LONG *maxn,
												IDL_LONG *nk,
												char filterfile[])
{
	void **filter_struct=NULL;
	IDL_LONG i,j,tmp_filter_n;
	float *tmp_filter_lambda=NULL,*tmp_filter_pass=NULL;
	FILE *listfp;
	char filterdirenv[1000],filterdir[1000],filter[1000],fullfilter[1000];

	FREEVEC((*filter_n));
	FREEVEC((*filter_lambda));
	FREEVEC((*filter_pass));
	
	/* read in once to determine number */
	listfp=k_fileopen(filterfile,"r");
	fscanf(listfp,"%s",filterdirenv);
	strcpy(filterdir,getenv(filterdirenv));
	fscanf(listfp,"%s",filter);
	(*nk)=0;
	(*maxn)=0;
	while(!feof(listfp)) {
		sprintf(fullfilter,"%s/%s",filterdir,filter);
		k_yanny_readone(fullfilter,&filter_struct,&tmp_filter_n,
										(void *) k_add_row_filter_struct);
		k_free_filter_struct(&filter_struct,tmp_filter_n);
		if(tmp_filter_n>(*maxn)) (*maxn)=tmp_filter_n;
		(*nk)++;
		fscanf(listfp,"%s",filter);
	}
	fclose(listfp);
	
	(*filter_n)=(IDL_LONG *) malloc((*nk)*sizeof(IDL_LONG));
	(*filter_lambda)=(float *) malloc((*nk)*(*maxn)*sizeof(float));
	(*filter_pass)=(float *) malloc((*nk)*(*maxn)*sizeof(float));
	listfp=k_fileopen(filterfile,"r");
	fscanf(listfp,"%s",filterdirenv);
	strcpy(filterdir,getenv(filterdirenv));
	fscanf(listfp,"%s",filter);
	for(i=0;i<(*nk);i++) {
		sprintf(fullfilter,"%s/%s",filterdir,filter);
		k_yanny_readone(fullfilter,&filter_struct,&tmp_filter_n,
										(void *) k_add_row_filter_struct);
		(*filter_n)[i]=tmp_filter_n;
		for(j=0;j<(*maxn);j++) {
			(*filter_lambda)[i*(*maxn)+j]=0.;
			(*filter_pass)[i*(*maxn)+j]=0.;
		} /* end for j */
		k_copy_filter_struct(filter_struct,(*filter_n)[i], 
												 &((*filter_lambda)[i*(*maxn)]), 
												 &((*filter_pass)[i*(*maxn)]));
		FREEVEC(tmp_filter_lambda);
		FREEVEC(tmp_filter_pass);
		k_free_filter_struct(&filter_struct,tmp_filter_n);
    fscanf(listfp,"%s",filter);
	} /* end for i */
	fclose(listfp);
	
	return(1);
} /* end k_load_ascii_table */
