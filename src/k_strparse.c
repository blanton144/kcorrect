#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <kcorrect.h>

#define STRMAX 255

/* allocate a char matrix with subscript range m[nrl..nrh][ncl..nch] */
char **k_cmatrix(long a, long b)
{
	char **list;
	long i;
 
	list=(char **) malloc((size_t) a*sizeof(char *));
	for(i=0;i<a;i++) 
		list[i]=(char *) malloc((size_t) b*sizeof(char));
	return(list);
}

/* free a char matrix allocated by cmatrix() */
void k_free_cmatrix(char **m, long a, long b)
{
	int i;
	if(m!=NULL) {
		for(i=0;i<a;i++) 
			if(m[i]!=NULL) free(m[i]);
		free(m);
	}
}

/*
 * strparse()
 *
 * Parse a string into separate words.  Do allocate 
 * the wrd array. Use the characters in the string
 * signal as the delimiting characters. 
 *
 * CRK 6/2/94
 *
 */
void k_strparse(char *str,
								char *signal,
								int *nwrd,
								char ***wrd)
{
	int i,istr,iword,len;

	/* count words and allocate memory */

	(*nwrd)=0;
	istr=0;
	len=strlen(str);
	while(istr<len) {
		while((strchr(signal,str[istr])) && istr<len) istr++;
		if(istr>=len) break;
		iword=0;
		while((!strchr(signal,str[istr])) && istr<len) {
			istr++;
			iword++;
		} /* end while */
		(*nwrd)++;
	} /* end while */

	if((*nwrd)==0) return;

	/* allocate words */
	(*wrd)=k_cmatrix((*nwrd),STRMAX+1);

	i=0;
	istr=0;
	iword=0;
	len=strlen(str);
	while(istr<len) {
		while((strchr(signal,str[istr])) && istr<len) istr++;
		if(istr>=len) break;
		iword=0;
		while((!strchr(signal,str[istr])) && istr<len) {
			(*wrd)[i][iword]=str[istr];
			iword++;
			istr++;
		} /* end while */
		(*wrd)[i][iword]='\0';
		i++;
	} /* end wHile */
} /* end strparse */

/* free a parsed string */
void k_strfree(char **str, int nwrd)
{
	k_free_cmatrix(str,nwrd,STRMAX+1);
} /* end strfree */

