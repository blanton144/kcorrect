#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <kcorrect.h>

/*
 * reconstruct_maggies.c
 *
 * Calculates the reconstructed flux, given the coefficients and the template
 * information; one can specify two parameters which characterize the 
 * magnitude:
 *
 * band_z : the shift applied to the band pass
 * at_z   : the redshift at which to observe the galaxy
 *
 * To use reconstruct_maggies as part of a system to produce absolute
 * magnitudes, use at_z=0. and set band_z to the shifted bandpass you
 * want to K-correct to. Eg. to K-correct to ^{0.2}b (in the notation
 * of the paper) call it like:
 *
 *  cat coeffs.dat | reconstruct_maggies 0.2 0. > maggies.dat
 *
 * where we specify at_z = 0. because this guarantees that we 
 * observe the galaxy at z=0. (appropriate for K-correction defn).
 *
 * To reconstruct the observed magnitudes at the observed redshift:
 *
 *  cat coeffs.dat | reconstruct_maggies > maggies.dat
 *
 * (Assuming coeffs.dat is actually output of fit_coeffs)
 *
 * Mike Blanton
 * 1/2002 */

static IDL_LONG nz=1000;
static double zmin=1.e-4, zmax=1.e-0;

static double *reconstruct_maggies=NULL;
static double *bmatrix=NULL;
static double *lambda=NULL;
static double *rmatrix=NULL;
static double *ematrix=NULL;  
static double *zvals=NULL;
static IDL_LONG nk,nb,nl,nt;
static double *filter_lambda=NULL;
static double *filter_pass=NULL;
static IDL_LONG *filter_n=NULL;
static IDL_LONG maxn;
static double *galaxy_z=NULL;
static double *band_shift=NULL;
static double *coeffs=NULL;

int main(int argc,
				 char **argv)
{
	double at_z=-1.,band_z=0.;
	IDL_LONG i,j,k,*sizes=NULL,ndim,nchunk,ncurrchunk;
	char ematrixfile[2000],bmatrixfile[2000],filterlist[2000],lambdafile[2000];
	char version[1000],versionpath[2000];

	/* set defaults */
	strcpy(version,"default");
	strcpy(versionpath,".");
	if(getenv("KCORRECT_DIR")!=NULL) {
		sprintf(versionpath,"%s/data/etemplates",getenv("KCORRECT_DIR"));
	} /* end if */

	/* read arguments */
	if(argc<0) {
		fprintf(stderr,"Usage: cat <coeff file> | reconstruct_maggies [band_z] [at_z] [version [version path]]\n");
		exit(1);
	} /* end if */
	i=1;
	if(argc>=2) 
		band_z=atof(argv[i]); i++; 
	if(argc>=3) 
		at_z=atof(argv[i]); i++; 
	if(argc>=4) 
		strcpy(version,argv[i]); i++;
	if(argc>=5) 
		strcpy(versionpath,argv[i]); i++;

	sprintf(ematrixfile,"%s/ematrix.%s.dat",versionpath,version);
	sprintf(bmatrixfile,"%s/bmatrix.%s.dat",versionpath,version);
	sprintf(lambdafile,"%s/lambda.%s.dat",versionpath,version);
	sprintf(filterlist,"%s/filterlist.%s.dat",versionpath,version);

	/* load the ematrix; this gives the eigentemplates in terms of the
	 * basis set */
	k_load_ascii_table(&ematrix,&ndim,&sizes,ematrixfile);
	nt=sizes[0];
	nb=sizes[1];
	FREEVEC(sizes);

	/* input the bmatrix; this is the spectrum for each basis element */
	k_load_ascii_table(&bmatrix,&ndim,&sizes,bmatrixfile);
	if(sizes[0]!=nb) {
		fprintf(stderr,"ematrix and bmatrix files incompatible.\n");
		exit(1);
	} /* end if */
	nl=sizes[1];
	FREEVEC(sizes);

	/* input the lambda; this is the wavelength scale for the bmatrix */
	k_load_ascii_table(&lambda,&ndim,&sizes,lambdafile);
	if(sizes[0]!=nl+1) {
		fprintf(stderr,"bmatrix and lambda files incompatible.\n");
		exit(1);
	} /* end if */
	FREEVEC(sizes);

	/* load the filters; these are the response curves, they can have
     arbitrary normalization */
	k_load_filters(&filter_n,&filter_lambda,&filter_pass,&maxn,&nk,filterlist);

	/* create the rmatrix; this is a big matrix which tabulates the
     projection of each basis element b onto each filter k, as a
     function of redshift; you only have to project onto the filters
     here, since every other spectrum you will project onto the
     filters will be a linear combination of the basis elements b; you
     interpolate the rmatrix to get a particular redshift */
	rmatrix=(double *) malloc(nz*nb*nk*sizeof(double));
	zvals=(double *) malloc(nz*sizeof(double));
	for(i=0;i<nz;i++)
		zvals[i]=zmin+(zmax-zmin)*((double)i+0.5)/(double)nz;
	k_create_r(rmatrix,nk,nb,bmatrix,lambda,nl,zvals,nz,filter_n,
						 filter_lambda,filter_pass,maxn);

	/* get the reconstructed fluxes given the coefficients; if you wanted, you
     could put the call to k_reconstruct_maggies directly after k_fit_coeffs
     in fit_coeffs.c; this is what you would do if you wanted to
     calculate K-corrections within a C code of your own */
	nchunk=2;
	reconstruct_maggies=(double *) malloc(nk*(nchunk+1)*sizeof(double));
	galaxy_z=(double *) malloc((nchunk+1)*sizeof(double));
	band_shift=(double *) malloc((nchunk+1)*sizeof(double));
	coeffs=(double *) malloc((nchunk+1)*nt*sizeof(double));
	fscanf(stdin,"%lf",&(coeffs[0]));
	while(!feof(stdin)) {
		for(i=0;i<nchunk && !feof(stdin);i++) {
			for(j=1;j<nt;j++)
				fscanf(stdin,"%lf",&(coeffs[i*nt+j]));
			fscanf(stdin,"%lf",&(galaxy_z[i]));
			if(at_z!=-1.) galaxy_z[i]=at_z;
			band_shift[i]=band_z;
			fscanf(stdin,"%lf",&(coeffs[(i+1)*nt+0]));
		} /* end for i */
		ncurrchunk=i;
		k_reconstruct_maggies(ematrix,nt,zvals,nz,rmatrix,nk,nb,coeffs,galaxy_z,
													band_shift,reconstruct_maggies,ncurrchunk);
		for(i=0;i<ncurrchunk;i++) {
			for(k=0;k<nk;k++)
				fprintf(stdout,"%e ",reconstruct_maggies[i*nk+k]);
			fprintf(stdout,"\n");
		} /* end for i */
		coeffs[0]=coeffs[ncurrchunk*nt+0];
	}
	
	FREEVEC(zvals);
	FREEVEC(rmatrix);
	FREEVEC(bmatrix);
	FREEVEC(ematrix);
	FREEVEC(lambda);
	FREEVEC(filter_n);
	FREEVEC(filter_lambda);
	FREEVEC(filter_pass);
	FREEVEC(galaxy_z);
	FREEVEC(reconstruct_maggies);
	FREEVEC(coeffs);

	return(0);
} /* end main */
