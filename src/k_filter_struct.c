#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <kcorrect.h>

/*
 * k_filter_struct
 *
 * Deal with the filter structure
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
