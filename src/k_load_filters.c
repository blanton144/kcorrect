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

IDL_LONG k_load_filter(char filename[],
											 IDL_LONG *filter_n,
											 double **filter_lambda,
											 double **filter_pass)
{
	char line[100000],format[100000];
	int started,i,j,nlines,ncolumns,defaultcol;
	FILE *fp;
	
	started=0;
	fp=k_fileopen(filename,"r");
	fgets(line,100000,fp);
	j=0;
	while(!feof(fp)) {
		if(line[0]!='#') {
			if(!started) {
				sscanf(line,"%d %d %d",&nlines,&ncolumns,&defaultcol);
				(*filter_n)=nlines;
				(*filter_lambda)=(double *) malloc(nlines*sizeof(double));
				(*filter_pass)=(double *) malloc(nlines*sizeof(double));
				started=1;
				j=0;
			} else if(j<(*filter_n)) {
				sprintf(format,"%%lf");
				for(i=1;i<defaultcol-2;i++) 
					sprintf(format,"%s%%*lf ",format);
				sprintf(format,"%s%%lf",format);
				sscanf(line,format,&((*filter_lambda)[j]),&((*filter_pass)[j]));
				j++;
			} /* end if..else */
		} /* end if */
		fgets(line,100000,fp);
	} /* end while */
	fclose(fp);

	return(1);
}


/* load the filters in their standard format */
IDL_LONG k_load_filters(IDL_LONG **filter_n,
												double **filter_lambda,
												double **filter_pass,
												IDL_LONG *maxn,
												IDL_LONG *nk,
												char filterlist[])
{
	IDL_LONG i,j,tmp_filter_n;
	double *tmp_filter_lambda,*tmp_filter_pass;
	FILE *listfp;
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
		k_load_filter(fullfilename,&tmp_filter_n,&tmp_filter_lambda, 
									&tmp_filter_pass);
		if(tmp_filter_n>(*maxn)) (*maxn)=tmp_filter_n;
		FREEVEC(tmp_filter_lambda);
		FREEVEC(tmp_filter_pass);
	} /* end for i */
	fclose(listfp);
		
	(*filter_lambda)=(double *) malloc((*nk)*(*maxn)*sizeof(double));
	(*filter_pass)=(double *) malloc((*nk)*(*maxn)*sizeof(double));
	listfp=k_fileopen(filterlist,"r");
	fscanf(listfp,"%d",nk);
	for(i=0;i<(*nk);i++) {
		fscanf(listfp,"%s",filename);
		sprintf(fullfilename,"%s/%s",filterdir,filename);
		k_load_filter(fullfilename,&tmp_filter_n,&tmp_filter_lambda, 
									&tmp_filter_pass);
		(*filter_n)[i]=tmp_filter_n;
		for(j=0;j<(*maxn);j++) {
			(*filter_lambda)[i*(*maxn)+j]=0.;
			(*filter_pass)[i*(*maxn)+j]=0.;
		} /* end for j */
		for(j=0;j<(*filter_n)[i];j++) {
			(*filter_lambda)[i*(*maxn)+j]=tmp_filter_lambda[j];
			(*filter_pass)[i*(*maxn)+j]=tmp_filter_pass[j];
		} /* end for j */
		FREEVEC(tmp_filter_lambda);
		FREEVEC(tmp_filter_pass);
	} /* end for i */
	fclose(listfp);
	
	return(1);
} /* end k_load_ascii_table */
