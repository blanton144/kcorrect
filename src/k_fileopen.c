#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * k_fileopen.c : contains routine to open a file 
 *              Created 9/28/95 MB
 */

/* open a file and return file pointer */
FILE *k_fileopen(const char *Filename,
								 const char *Mode)
{
	FILE *fp;
	
	if (!(fp=fopen(Filename,Mode))) {
		fprintf(stderr,"Could not open file for %s: %s",Mode,Filename);
		exit(1);
	} /* end if */

	return(fp);
} /* end fileopen() */
