#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <kcorrect.h>

/*
 * k_filter_struct
 *
 * Deal with the filter structure. 
 *
 * Mike Blanton
 * 6/2003 
 */

typedef struct {
	double lambda;
	double pass;
} FILTER_STRUCT;

void k_print_filter_struct(void **input_struct, 
													 int nrows) 
{
	int i;
	FILTER_STRUCT **ptr_struct;

	ptr_struct=(FILTER_STRUCT **) input_struct;
	for(i=0;i<nrows;i++) 
		printf("%e %e\n",ptr_struct[i]->lambda,ptr_struct[i]->pass);
}

void k_add_row_filter_struct(void ***input_struct,
														 int *nrows, 
														 char *columns[], 
														 int ncolumns,
														 char **wrd,
														 int nwrds)
{
	int i;
	FILTER_STRUCT **ptr_struct;

	ptr_struct=(FILTER_STRUCT **) (*input_struct);
	(*nrows)++;
	ptr_struct=(FILTER_STRUCT **) realloc((void **) ptr_struct, 
																				(*nrows)*sizeof(FILTER_STRUCT *));
	ptr_struct[(*nrows)-1]=(FILTER_STRUCT *) malloc(sizeof(FILTER_STRUCT));
	for(i=0;i<nwrds;i++) {
		if(!strcmp("lambda",columns[i])) 
			ptr_struct[(*nrows)-1]->lambda=atof(wrd[i]);
		if(!strcmp("pass",columns[i])) 
			ptr_struct[(*nrows)-1]->pass=atof(wrd[i]);
	} /* end for */
	(*input_struct)=(void **) ptr_struct;
}

void k_copy_filter_struct(void **filter_struct,
													int filter_n,
													float *filter_lambda,
													float *filter_pass)
{
	int i;
	FILTER_STRUCT **ptr_struct;

	ptr_struct=(FILTER_STRUCT **) filter_struct;
	for(i=0;i<filter_n;i++) {
		filter_lambda[i]=ptr_struct[i]->lambda;
		filter_pass[i]=ptr_struct[i]->pass;
	}
}

void k_free_filter_struct(void ***input_struct,
													int nrows) 
{
	int i;
	if(input_struct!=NULL) {
		if(*input_struct!=NULL) {
			for(i=0;i<nrows;i++) {
				if((*input_struct)[i]!=NULL)
					free((*input_struct)[i]);
				(*input_struct)[i]=NULL;
			}
			free((*input_struct));
		}
		(*input_struct)=NULL;
	}
} 
