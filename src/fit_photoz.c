#include <stdio.h>
#include <string.h>
#include <math.h>
#include <kcorrect.h>

/*
 * fit_photoz.c
 *
 * Given the locations of the eigentemplate matrix, the basis set matrix,
 * and galaxy fluxes and errors, output the best-fit redshift.
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
	IDL_LONG i,j,k,*sizes=NULL,ndim,nchunk;
	char ematrixfile[2000],bmatrixfile[2000],filterlist[2000],lambdafile[2000];

	if(argc<5) {
		fprintf(stderr,"Usage: cat <galaxy file> | fit_photoz <ematrix file> <bmatrix file> <lambdafile> <filterlist> [zmin] [zmax] [nz]\n");
		exit(1);
	} /* end if */
	i=1;
	strcpy(ematrixfile,argv[i]); i++;
	strcpy(bmatrixfile,argv[i]); i++;
	strcpy(lambdafile,argv[i]); i++;
	strcpy(filterlist,argv[i]); i++;

	/* load the ematrix */
	k_load_ascii_table(&ematrix,&ndim,&sizes,ematrixfile);
	nt=sizes[0];
	nb=sizes[1];
	FREEVEC(sizes);

	/* input the bmatrix */
	k_load_ascii_table(&bmatrix,&ndim,&sizes,bmatrixfile);
	if(sizes[0]!=nb) {
		fprintf(stderr,"ematrix and bmatrix files incompatible.\n");
		exit(1);
	} /* end if */
	nl=sizes[1];
	FREEVEC(sizes);

	/* input the lambda */
	k_load_ascii_table(&lambda,&ndim,&sizes,lambdafile);
	if(sizes[0]!=nl+1) {
		fprintf(stderr,"bmatrix and lambda files incompatible.\n");
		exit(1);
	} /* end if */
	FREEVEC(sizes);

	/* load the filters */
	k_load_filters(&filter_n,&filter_lambda,&filter_pass,&maxn,&nk,filterlist);

	/* create the rmatrix */
	rmatrix=(double *) malloc(nz*nb*nk*sizeof(double));
	zvals=(double *) malloc(nz*sizeof(double));
	for(i=0;i<nz;i++)
		zvals[i]=zmin+(zmax-zmin)*((double)i+0.5)/(double)nz;
	k_create_r(rmatrix,nk,nb,bmatrix,lambda,nl,zvals,nz,filter_n,
						 filter_lambda,filter_pass,maxn);

	/* get the redshift estimates */
	nchunk=1;
	galaxy_z=(double *) malloc(nchunk*sizeof(double));
	galaxy_maggies=(double *) malloc(nchunk*nk*sizeof(double));
	galaxy_invvar=(double *) malloc(nchunk*nk*sizeof(double));
	coeffs=(double *) malloc(nchunk*nt*sizeof(double));
	fscanf(stdin,"%lf",&(galaxy_maggies[0]));
	while(!feof(stdin)) {
		for(k=1;k<nk;k++)
			fscanf(stdin,"%lf",&(galaxy_maggies[k]));
		for(k=0;k<nk;k++)
			fscanf(stdin,"%lf",&(galaxy_invvar[k]));
		k_fit_photoz(ematrix,nt,zvals,nz,rmatrix,nk,nb,coeffs,galaxy_maggies,
								 galaxy_invvar,galaxy_z,nchunk);
		for(j=0;j<nt;j++)
			fprintf(stdout,"%e ",coeffs[j]);
		fprintf(stdout,"%e\n",galaxy_z[0]);
		fscanf(stdin,"%lf",&(galaxy_maggies[0]));
	} /* end for while */
	
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
