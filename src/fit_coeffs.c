#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <kcorrect.h>

/*
 * fit_coeffs.c
 *
 * Given the locations of the eigentemplate matrix, the basis set matrix,
 * and galaxy maggies, errors, and redshift, output the best fit coefficients.
 *
 * Mike Blanton
 * 2/2001
 */

static IDL_LONG nz=1000;
static double zmin=1.e-4, zmax=1.e-0;

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
static double *galaxy_maggies=NULL;
static double *galaxy_invvar=NULL;
static double *coeffs=NULL;

int main(int argc,
				 char **argv)
{
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
		fprintf(stderr,"Usage: cat <galaxy file> | fit_coeffs [version [version path]]\n");
		exit(1);
	} /* end if */
	i=1;
	if(argc>=2) 
		strcpy(version,argv[i]); i++;
	if(argc>=3) 
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

	/* get the coefficients; sequentially calls k_fit_coeffs to get the
     coefficients */
	nchunk=20;
	galaxy_z=(double *) malloc((nchunk+1)*sizeof(double));
	galaxy_maggies=(double *) malloc((nchunk+1)*nk*sizeof(double));
	galaxy_invvar=(double *) malloc((nchunk+1)*nk*sizeof(double));
	coeffs=(double *) malloc((nchunk+1)*nt*sizeof(double));
	fscanf(stdin,"%lf",&(galaxy_z[0]));
	while(!feof(stdin)) {
		for(i=0;i<nchunk && !feof(stdin);i++) {
			for(k=0;k<nk;k++)
				fscanf(stdin,"%lf",&(galaxy_maggies[i*nk+k]));
			for(k=0;k<nk;k++)
				fscanf(stdin,"%lf",&(galaxy_invvar[i*nk+k]));
			fscanf(stdin,"%lf",&(galaxy_z[i+1]));
		} /* end for i */
		ncurrchunk=i;
		/* no direct constraints on the coeffs are included in this fit */
		k_fit_coeffs(ematrix,nt,zvals,nz,rmatrix,nk,nb,coeffs,galaxy_maggies,
								 galaxy_invvar,galaxy_z,NULL,NULL,NULL,0,ncurrchunk);
		for(i=0;i<ncurrchunk;i++) {
			for(j=0;j<nt;j++)
				fprintf(stdout,"%e ",coeffs[i*nt+j]);
			fprintf(stdout,"%e\n",galaxy_z[i]);
		} /* end for i */
		galaxy_z[0]=galaxy_z[ncurrchunk];
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
	FREEVEC(galaxy_maggies);
	FREEVEC(galaxy_invvar);
	FREEVEC(coeffs);
} /* end main */
