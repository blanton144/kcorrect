#include <stdio.h>
#include <string.h>
#include <math.h>
#include <kcorrect.h>

/*
 * reconstruct_maggies.c
 *
 * Calculates the reconstructed flux, given the coefficients and the template
 * information
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
static double *coeffs=NULL;

int main(int argc,
				 char **argv)
{
	double to_z=-1.;
	IDL_LONG i,j,k,*sizes=NULL,ndim,nchunk;
	char ematrixfile[2000],bmatrixfile[2000],filterlist[2000],lambdafile[2000];

	if(argc<5) {
		fprintf(stderr,"Usage: cat <coeff file> | reconstruct_maggies <ematrix file> <bmatrix file> <lambdafile> <filterlist> [to_z] [zmin] [zmax] [nz]\n");
		exit(1);
	} /* end if */
	i=1;
	strcpy(ematrixfile,argv[i]); i++;
	strcpy(bmatrixfile,argv[i]); i++;
	strcpy(lambdafile,argv[i]); i++;
	strcpy(filterlist,argv[i]); i++;
	if(argc>=6) { to_z=atof(argv[i]); i++; }

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
	nchunk=1;
	reconstruct_maggies=(double *) malloc(nk*nchunk*sizeof(double));
	galaxy_z=(double *) malloc(nchunk*sizeof(double));
	coeffs=(double *) malloc(nchunk*nt*sizeof(double));
	fscanf(stdin,"%lf",&(coeffs[0]));
	while(!feof(stdin)) {
		for(j=1;j<nt;j++)
			fscanf(stdin,"%lf",&(coeffs[j]));
		fscanf(stdin,"%lf",&(galaxy_z[0]));
		if(to_z!=-1.) galaxy_z[0]=to_z;
		k_reconstruct_maggies(ematrix,nt,zvals,nz,rmatrix,nk,nb,coeffs,galaxy_z,
												 reconstruct_maggies,nchunk);
		for(k=0;k<nk;k++)
			fprintf(stdout,"%e ",reconstruct_maggies[k]);
		fprintf(stdout,"\n");
		fscanf(stdin,"%lf",&(coeffs[0]));
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
} /* end main */
