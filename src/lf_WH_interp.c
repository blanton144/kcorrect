#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*
 * WH_interp.c: do W and H function for interpolation
 *
 * interp_choice == 1 means step-wise interpolation
 * interp_choice == 2 means linear interpolation (NOT IMPLEMENTED)
 */

static int WH_interp_choice;
static int WH_nparam;
static float *WH_lum=NULL;
static float WH_lumbright, WH_lumdim;
static float WH_loglumbright, WH_loglumdim;

void set_WH_interp(int interp_choice,
									 float Mdim,
									 float Mbright,
									 float lum[],
									 float M[],
									 int nparam)
{
	int i;

	if(interp_choice==1)
		WH_interp_choice=1;
	else if (WH_interp_choice==2) {
		printf("Linear interpolation not implemented yet!\n"); 
		exit(1);
	} else {
		printf("No such interpolation choice (%d)!\n",interp_choice); 
		exit(1);
	} /* end if */

	WH_nparam=nparam;
	WH_loglumbright=-0.4*(Mbright+20.);
	WH_lumbright=pow(10.,WH_loglumbright);
	WH_loglumdim=-0.4*(Mdim+20.);
	WH_lumdim=pow(10.,WH_loglumdim);

	WH_lum=(float *) malloc((WH_nparam+1)*sizeof(float));
	for(i=0;i<=WH_nparam;i++) {
		WH_lum[i]=WH_loglumdim+(WH_loglumbright-WH_loglumdim)*(float)i
			/(float)WH_nparam;
		lum[i]=WH_lum[i]=pow(10.,WH_lum[i]);
		M[i]=-20.-2.5*log10(lum[i]);
	} /* end for */
} /* end set_WH_interp */

void unset_WH_interp(void)
{
	free((char *)WH_lum);
  WH_lum=NULL;
} /* end unset_WH_interp */

float W(float currlum,
				 int indx,
				 float factor)
		 
{
	float val,lowlum,highlum;

	if (WH_interp_choice==1) {
		lowlum=WH_lum[indx]*factor;
		if(lowlum<WH_lum[0]) lowlum=WH_lum[0];
		if(lowlum>WH_lum[WH_nparam]) lowlum=WH_lum[WH_nparam];
		highlum=WH_lum[indx+1]*factor;
		if(highlum<WH_lum[0]) highlum=WH_lum[0];
		if(highlum>WH_lum[WH_nparam]) highlum=WH_lum[WH_nparam];
		if(currlum>lowlum && currlum<=highlum)
			val=1.;
		else 
			val=0.;
		return(val);
	} else if (WH_interp_choice==2) {
		printf("Linear interpolation not implemented yet!\n"); 
		exit(1);
	} else {
		printf("No such interpolation choice (%d)!\n",WH_interp_choice); 
		exit(1);
	} /* end if */

	return(1.e+30);
} /* end W */

float H(float currlum,
				 int indx,
				 float factor)
{
	float val,lowlum,highlum;

	if (WH_interp_choice==1) {
		lowlum=WH_lum[indx]*factor;
		if(lowlum<WH_lum[0]) lowlum=WH_lum[0];
		if(lowlum>WH_lum[WH_nparam]) lowlum=WH_lum[WH_nparam];
		highlum=WH_lum[indx+1]*factor;
		if(highlum<WH_lum[0]) highlum=WH_lum[0];
		if(highlum>WH_lum[WH_nparam]) highlum=WH_lum[WH_nparam];

		if(currlum>highlum)
			val=0.;
		else if (currlum>lowlum && currlum<=highlum)
			val=highlum-currlum;
		else 
			val=highlum-lowlum;
		return(val);
	} else if (WH_interp_choice==2) {
		printf("Linear interpolation not implemented yet!\n"); 
		exit(1);
	} else {
		printf("No such interpolation choice (%d)!\n",WH_interp_choice); 
		exit(1);
	} /* end if */
	
	return(1.e+30);
} /* end H */
