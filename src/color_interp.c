#include <stdio.h>
#include <string.h>
#include <math.h>
#include "utils.h"
#include "color_interp.h"

/*
 * color_interp.c
 *
 * This code is an implementation of the dumbest possible 
 * interpolation scheme between colors.
 *
 * Mike Blanton
 * 6/2001
 */

#ifndef ZZERO
#define ZZERO 0.
#endif

static int band,currband,igal;
static int band0,band1;
static double mag0,mag1;
static double ci_A,ci_alpha,ci_z,ci_break;
static double filt[5][300];
static double lfilt[5][300];
static int nfilt[5];

/* Galaxy properties */
static double galaxy_mag[5][NALL];
static double galaxy_rfmag[5][NALL];
static double galaxy_magerr[5][NALL];
static double galaxy_z[NALL];
static int galaxy_use[NALL];

static double filtfunc(double lambda)
{
	unsigned long i,ip1;
	double s,val;

	locate(lfilt[currband], nfilt[currband], lambda, &i);
	if(i>=nfilt[currband]-1 || i<0)
		return(0.);
		
	ip1=i+1;
	s=lambda-lfilt[currband][i];
	val=filt[currband][i]+s*(filt[currband][ip1]-filt[currband][i])/
		(lfilt[currband][ip1]-lfilt[currband][i]);
	return(val);
} /* end filtfunc */

static double scalefiltfunc(double lambda)
{
	double val;
	val=filtfunc(lambda)/lambda;
	return(val);
} /* end tempfiltfunc */

static double tfunc(double lambda)
{
	unsigned long i,ip1;

	return(ci_A*pow(lambda/(1.+ci_z),ci_alpha)/(1.+ci_z));
} /* end tfunc */

static double tempfiltfunc(double lambda)
{
	double val;
	val=filtfunc(lambda)*tfunc(lambda)*lambda;
	return(val);
} /* end tempfiltfunc */

static double bfunc(double lambda)
{
	unsigned long i,ip1;

	if(lambda/(1.+ci_z)<4100.)
		return(ci_break*pow(lambda/(1.+ci_z),2.)/(1.+ci_z));
	else 
		return(ci_A*pow(lambda/(1.+ci_z),ci_alpha)/(1.+ci_z));
} /* end tfunc */

static double btempfiltfunc(double lambda)
{
	double val,bval;
	bval=bfunc(lambda);
	val=filtfunc(lambda)*bval*lambda;
	return(val);
} /* end tempfiltfunc */

double alpha_fit(double alpha) 
{
	double scaleval0,val0;
	double scaleval1,val1;
	double diff;

	ci_A=1.;
	ci_alpha=alpha;

	currband=band0;
	scaleval0=qromo(scalefiltfunc,lfilt[currband][0],
									lfilt[currband][nfilt[currband]-1],midpnt);
	val0=qromo(tempfiltfunc,lfilt[currband][0],
						 lfilt[currband][nfilt[currband]-1],midpnt);
	val0/=scaleval0;

	currband=band1;
	scaleval1=qromo(scalefiltfunc,lfilt[currband][0],
									lfilt[currband][nfilt[currband]-1],midpnt);
	val1=qromo(tempfiltfunc,lfilt[currband][0],
						 lfilt[currband][nfilt[currband]-1],midpnt);
	val1/=scaleval1;

	diff=(mag0-mag1)-2.5*log10(val1/val0);

	return(diff);
} /* end alpha_fit */

double balpha_fit(double alpha) 
{
	double scaleval0,val0;
	double scaleval1,val1;
	double diff;

	ci_A=1.;
	ci_alpha=alpha;

	currband=band1;
	scaleval1=qromo(scalefiltfunc,lfilt[currband][0],
									lfilt[currband][nfilt[currband]-1],midpnt);
	val1=qromo(btempfiltfunc,lfilt[currband][0],
						 lfilt[currband][nfilt[currband]-1],midpnt);
	val1/=scaleval1;
	ci_A=pow(10.,-0.4*(galaxy_mag[currband][igal]-17.))/val1;
	scaleval1=qromo(scalefiltfunc,lfilt[currband][0],
									lfilt[currband][nfilt[currband]-1],midpnt);
	val1=qromo(btempfiltfunc,lfilt[currband][0],
						 lfilt[currband][nfilt[currband]-1],midpnt);
	val1/=scaleval1;

	currband=band0;
	scaleval0=qromo(scalefiltfunc,lfilt[currband][0],
									lfilt[currband][nfilt[currband]-1],midpnt);
	val0=qromo(btempfiltfunc,lfilt[currband][0],
						 lfilt[currband][nfilt[currband]-1],midpnt);
	val0/=scaleval0;

	diff=(mag0-mag1)-2.5*log10(val1/val0);

	return(diff);
} /* end balpha_fit */

int main(int argc, char **argv)
{
	double currz,cz,redden,alpha,A,val,scaleval;
	int i,ii,j,k,l;
	char filename[255];
	FILE *fp;

	/*
	 * 0. Memory allocation 
	 */

	/*
	 * 1. Read in observed fluxes of each galaxy, observed redshifts,
	 *    and whether to use it in the sample; then make the sample
	 */

	/* 1a. Read in magnitudes, and reddening-adjust if desired */
	printf("Reading magnitudes ...\n");
	fflush(stdout);
	fp=fileopen(MAGFILE,"r");
	for(i=0;i<NALL;i++) {
		for(k=0;k<5;k++)
			fscanf(fp,"%lf",&(galaxy_mag[k][i]));
		for(k=0;k<5;k++) 
			fscanf(fp,"%lf",&(galaxy_magerr[k][i]));
	} /* end for i */
	fclose(fp);
#ifdef REDDENFILE
	printf("Reading reddening ...\n");
	fflush(stdout);
	fp=fileopen(REDDENFILE,"r");
	for(i=0;i<NALL;i++) 
		for(k=0;k<5;k++) {
			fscanf(fp,"%lf",&redden);
			galaxy_mag[k][i]-=redden;
		} /* end for k */
	fclose(fp);
#endif

	/* 1b. Read in redshifts */
	printf("Reading redshifts ...\n");
	fflush(stdout);
	fp=fileopen(POSITIONFILE,"r");
	for(i=0;i<NALL;i++) {
		fscanf(fp,"%*lf %*lf %lf %*d %*d",&cz);
		galaxy_z[i]=cz/LIGHTSPEED;
	} /* end for i */
	fclose(fp);
	printf("Reading use ...\n");
	fflush(stdout);

	/* 
	 * 2. Read in filters 
	 */
	for(k=0;k<5;k++) {
		/* read in the filter */
		sprintf(filename,"%s/%s%d.dat",FILTERPATH,FILTERBASE,k);
		printf("Processing filter %s ...\n",filename);
		fflush(stdout);
		fp=fileopen(filename,"r");
		fscanf(fp,"%d",&(nfilt[k]));
		printf("     (%d elements)\n",nfilt[k]);
		fflush(stdout);
		for(i=0;i<nfilt[k];i++)
			fscanf(fp,"%lf %lf",&(lfilt[k][i]),&(filt[k][i]));
		fclose(fp);
	} /* end for k */

	/* 
	 * 3. For each galaxy, interpolate to each restframe color 
	 */
	for(i=0;i<NALL;i++) {
		igal=i;
		for(band=0;band<5;band++) {
#ifdef INCLUDEBREAK
			if(band==0) {
				ci_A=1.;
				ci_break=1.;
				ci_alpha=0.;
				ci_z=galaxy_z[i];
				
				currband=band;
				ci_z=galaxy_z[i];
				scaleval=qromo(scalefiltfunc,lfilt[band][0],
											 lfilt[band][nfilt[band]-1],midpnt);
				val=qromo(btempfiltfunc,lfilt[band][0],lfilt[band][nfilt[band]-1],
									midpnt);
				val/=scaleval;
				ci_break=pow(10.,-0.4*(galaxy_mag[band][i]-17.))/val;
				
				ci_z=ZZERO;
				currband=band;
				scaleval=qromo(scalefiltfunc,lfilt[band][0],
											 lfilt[band][nfilt[band]-1],midpnt);
				val=qromo(btempfiltfunc,lfilt[band][0],lfilt[band][nfilt[band]-1],
									midpnt);
				val/=scaleval;
				galaxy_rfmag[band][i]=17.-2.5*log10(val);
				galaxy_rfmag[band][i]=val;
			} else if(band==1) {
				band0=band;
				band1=band+1;
				mag0=galaxy_mag[band0][i];
				mag1=galaxy_mag[band1][i];
				ci_A=1.;
				ci_z=galaxy_z[i];
				alpha=zbrent(balpha_fit,-18.,18.,1.e-6);
				ci_A=1.;
				ci_alpha=alpha;
				currband=band1;
				scaleval=qromo(scalefiltfunc,lfilt[band1][0],
											 lfilt[band1][nfilt[band1]-1], midpnt);
				val=qromo(btempfiltfunc,lfilt[band1][0],lfilt[band1][nfilt[band1]-1],
									midpnt);
				val/=scaleval;
				A=pow(10.,-0.4*(galaxy_mag[band1][i]-17.))/val;
				
				ci_z=ZZERO;
				ci_A=A;
				currband=band;
				scaleval=qromo(scalefiltfunc,lfilt[band][0],
											 lfilt[band][nfilt[band]-1],midpnt);
				val=qromo(btempfiltfunc,lfilt[band][0],lfilt[band][nfilt[band]-1],
									midpnt);
				val/=scaleval;
				galaxy_rfmag[band][i]=17.-2.5*log10(val);
				galaxy_rfmag[band][i]=val;
			} else {
				if(band<4) {
					band0=band;
					band1=band+1;
				} else {
					band0=band-1;
					band1=band;
				} /* end if..else */
				mag0=galaxy_mag[band0][i];
				mag1=galaxy_mag[band1][i];
				ci_A=1.;
				ci_z=galaxy_z[i];
				alpha=zbrent(alpha_fit,-18.,18.,1.e-6);
				ci_alpha=alpha;
				currband=band;
				scaleval=qromo(scalefiltfunc,lfilt[band][0],
											 lfilt[band][nfilt[band]-1], midpnt);
				val=qromo(tempfiltfunc,lfilt[band][0],lfilt[band][nfilt[band]-1],
									midpnt);
				val/=scaleval;
				A=pow(10.,-0.4*(galaxy_mag[band][i]-17.))/val;
				
				ci_z=ZZERO;
				ci_A=A;
				currband=band;
				scaleval=qromo(scalefiltfunc,lfilt[band][0],
											 lfilt[band][nfilt[band]-1],midpnt);
				val=qromo(tempfiltfunc,lfilt[band][0],lfilt[band][nfilt[band]-1],
									midpnt);
				val/=scaleval;
				galaxy_rfmag[band][i]=17.-2.5*log10(val);
				galaxy_rfmag[band][i]=val;
			} /* end if */
#else 
			if(band<4) {
				band0=band;
				band1=band+1;
			} else {
				band0=band-1;
				band1=band;
			} /* end if..else */
			mag0=galaxy_mag[band0][i];
			mag1=galaxy_mag[band1][i];
			ci_A=1.;
			ci_z=galaxy_z[i];
			alpha=zbrent(alpha_fit,-18.,18.,1.e-6);
			ci_alpha=alpha;
			currband=band;
			scaleval=qromo(scalefiltfunc,lfilt[band][0],lfilt[band][nfilt[band]-1],
										 midpnt);
			val=qromo(tempfiltfunc,lfilt[band][0],lfilt[band][nfilt[band]-1],
								midpnt);
			val/=scaleval;
			A=pow(10.,-0.4*(galaxy_mag[band][i]-17.))/val;

			ci_z=ZZERO;
			ci_A=A;
			currband=band;
			scaleval=qromo(scalefiltfunc,lfilt[band][0],lfilt[band][nfilt[band]-1],
										 midpnt);
			val=qromo(tempfiltfunc,lfilt[band][0],lfilt[band][nfilt[band]-1],
								midpnt);
			val/=scaleval;
			galaxy_rfmag[band][i]=17.-2.5*log10(val);
			galaxy_rfmag[band][i]=val;
#endif
		} /* end for band */
	} /* end for i */
	
	fp=fileopen(RFFILE,"w");
	for(i=0;i<NALL;i++) 
		fprintf(fp,"%le %le %le %le %le\n", galaxy_rfmag[0][i], 
						galaxy_rfmag[1][i],galaxy_rfmag[2][i],galaxy_rfmag[3][i], 
						galaxy_rfmag[4][i]);
	fclose(fp);
	
	printf("Done.\n");
	fflush(stdout);
	
} /* end main */

