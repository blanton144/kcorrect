#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <kcorrect.h>

/*
 * k_yanny_readone
 *
 * Code to read the Yanny parameter files and store them in a
 * particular structure. This code is not completely general and
 * depends on some coordination by the user. The structure passed in
 * must have all the desired columns defined and the desired columns
 * must be specified by an array of strings. Other columns in the file
 * are ignored. It will only work on a the first structure defined in
 * the file and it will ignore header information.
 *
 * Basically, it is what I need to read in the filters and nothing
 * more.
 *
 * Mike Blanton
 * 6/2003 
 */

IDL_LONG k_yanny_readone(char filename[255],
												 void ***input_struct,
												 int *nrows,
												 void (*add_row_struct)(void ***input_struct, 
																								int *nrows, 
																								char *columns[],
																								int ncolumns,
																								char **wrd,
																								int nwrds))
{
	int nwrds, ncolumns, typedef_done, typedef_doing;
	char str[255], struct_name[255];
	char **wrd, **columns=NULL;
	FILE *fp;

	fp=k_fileopen(filename,"r");
	typedef_done=0;
	typedef_doing=0;
	ncolumns=0;
	(*nrows)=0;
	while(!feof(fp)) {
		fgets(str,255,fp);

		/* ignore lines starting with # */ 
		if(str[0] != '#' && !feof(fp)) {
			str[strlen(str)-1]='\0';
			k_strparse(str,"\t ",&nwrds,&wrd);

			if(nwrds>0) {
				if(!typedef_done && !typedef_doing) {
					/* if no typedef yet, search for it, and pass if not done yet */
					if(!strcmp(wrd[0],"typedef")) {
						typedef_doing=1;
					} 
				} else if(typedef_doing) {
					/* if we are in the middle of typedefing, deal with that */
					if(!strcmp(wrd[0],"}")) {
						/* if we have reached the end, record the name */
						/* (get rid of semicolon first) */
						if(wrd[1][strlen(wrd[1])-1] == ';') wrd[1][strlen(wrd[1])-1]='\0';
						strcpy(struct_name,wrd[1]);
						typedef_done=1;
						typedef_doing=0;
					} else {
						/* if not, add a column */
						ncolumns++;
						if(columns==NULL) {
							columns=(char **) malloc(ncolumns*sizeof(char *));
						} else {
							columns=(char **) realloc((void *) columns, 
																				ncolumns*sizeof(char *));
						}
						/* get rid of semicolon */
						if(wrd[1][strlen(wrd[1])-1] == ';') wrd[1][strlen(wrd[1])-1]='\0';
						columns[ncolumns-1]=(char *) malloc(255*sizeof(char));
						strcpy(columns[ncolumns-1],wrd[1]);
					}
				} else if(typedef_done) { 
					/* if we have finished check structure name */
					if(!strcmp(wrd[0],struct_name)) 
						add_row_struct(input_struct, nrows, columns, ncolumns, &(wrd[1]), 
													 nwrds-1);
				} else {
					fprintf(stderr, 
									"k_yanny_readone.c: You really really shouldn't get here.\n");
					fflush(stderr);
					exit(1);
				}
			}
		} 
	}
	fclose(fp);
	
	k_strfree(columns,ncolumns);
	k_strfree(wrd,nwrds);
} /* end k_yanny_readone */
