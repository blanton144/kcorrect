#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <kcorrect.h>

/*
 * k_test_filter
 *
 * test filter reading
 *
 * Mike Blanton
 * 6/2003
 */

int main(int argc,
				 char **argv)
{
	int res,nrows;
	char filename[255];
	void **filter_struct=NULL;
	
	sprintf(filename,"/home/blanton/kcorrect/data/filters/twomass_Ks.par");
	res=k_yanny_readone(filename,&filter_struct, &nrows, 
											(void *) k_add_row_filter_struct);
	res=k_yanny_readone(filename,&filter_struct, &nrows, 
											(void *) k_add_row_filter_struct);
	res=k_yanny_readone(filename,&filter_struct, &nrows, 
											(void *) k_add_row_filter_struct);
	res=k_yanny_readone(filename,&filter_struct, &nrows, 
											(void *) k_add_row_filter_struct);
	res=k_yanny_readone(filename,&filter_struct, &nrows, 
											(void *) k_add_row_filter_struct);
	printf("%d\n",nrows);
	k_print_filter_struct(filter_struct, nrows); 
	return(0);
} /* end main */
