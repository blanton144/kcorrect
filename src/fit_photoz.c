#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <kcorrect.h>

/*
 * fit_photoz.c
 *
 * Given the locations of the vmatrix and lambda, the filter list,
 * and the galaxy maggies and invvar, return the best guess at the 
 * redshift.
 *
 * Mike Blanton
 * 7/2003
 */

static IDL_LONG nz=1000;
static IDL_LONG nprior=2;
static IDL_LONG maxiter=10000;
static float tolerance=1.e-6;
static float zmin=1.e-4, zmax=1.e-0;
static float band_shift=0.;

static float *lprior=NULL;
static float *zprior=NULL;
static float *lambda=NULL;
static float *vmatrix=NULL;
static float *rmatrix=NULL;
static float *zvals=NULL;
static IDL_LONG nk,nv,nl;
static float *filter_lambda=NULL;
static float *filter_pass=NULL;
static IDL_LONG *filter_n=NULL;
static IDL_LONG maxn;
static float *redshift=NULL;
static float *maggies=NULL;
static float *maggies_ivar=NULL;
static float *coeffs=NULL;
static float *chi2=NULL;

#define USAGE \
		{ fprintf(stderr,"Usage: cat <maggies file> | fit_photoz [--vfile <vfile> --lfile <lfile>\n"); \
		  fprintf(stderr,"            --ffile <ffile> ]\n"); }

int main(int argc,
				 char **argv)
{
	IDL_LONG i,j,k,c,ndim,niter,nchunk,ncurrchunk,*sizes=NULL;
	char vfile[2000],lfile[2000],ffile[2000],path[2000];
	char vmatrixfile[2000],lambdafile[2000],filterfile[2000];
	char priorfile[2000];
  FILE *fp;

	/* read arguments */
	strcpy(priorfile,"NULL");
	strcpy(vfile,"vmatrix.default.dat");
	strcpy(lfile,"lambda.default.dat");
	strcpy(ffile,"sdss_filters.dat");
	sprintf(path,"%s/data/templates",getenv("KCORRECT_DIR"));
	i=0;
	while(1) {
		int option_index = 0;
 		static struct option long_options[] =
			{
				{"vfile", 1, 0, 0}, 
				{"lfile", 1, 0, 0},
				{"path", 1, 0, 0},
				{"ffile", 1, 0, 0},
				{"prior", 1, 0, 0},
				{"help", 0, 0, 0}
			};
		static const char short_options[]="v:l:p:f:c:h";
		static const char short_options_c[]="vlpfch";

		c=getopt_long(argc, argv, short_options, long_options, &option_index);
		if(c==-1) break;
		if(c==0) c=short_options_c[option_index];
		switch(c) {
		case 'v':
			strcpy(vfile,optarg);
			break; 
		case 'l':
			strcpy(lfile,optarg);
			break; 
		case 'p':
			strcpy(path,optarg);
			break; 
		case 'f':
			strcpy(ffile,optarg);
			break; 
		case 'c':
			strcpy(priorfile,optarg);
			break; 
		case 'h':
      USAGE;
      exit(1);
			break; 
		case '?':
			break;
		default: 
			printf("fit_photoz: getopt returned character code 0%o ??\n", c);
		}
		i++;
	}
	if(argc<0) {
    USAGE;
		exit(1);
	} /* end if */

	/* get file names */
	sprintf(vmatrixfile,"%s/%s",path,vfile);
	sprintf(lambdafile,"%s/%s",path,lfile);
	sprintf(filterfile,"%s/%s",path,ffile);

  /* set up prior */
  if(strcmp(priorfile,"NULL")) {
    fp=k_fileopen(priorfile,"r");
    nprior=0;
    while(!feof(fp)) {
      fscanf(fp,"%*e %*e");
      nprior++;
    } /* end while */
    fclose(fp);
    lprior=(float *) malloc(nprior*sizeof(float));
    zprior=(float *) malloc(nprior*sizeof(float));
    fp=k_fileopen(priorfile,"r");
    nprior=0;
    while(!feof(fp)) {
      fscanf(fp,"%e %e", &(zprior[i]), &(lprior[i]));
      nprior++;
    } /* end while */
    fclose(fp);
  } else {
    nprior=2;
    lprior=(float *) malloc(nprior*sizeof(float));
    zprior=(float *) malloc(nprior*sizeof(float));
    zprior[0]=0.;
    zprior[1]=1000.;
  } /* end if..else */

	/* read in templates */
	k_read_ascii_table(&vmatrix,&ndim,&sizes,vmatrixfile);
	nl=sizes[1];
	nv=sizes[0];
	FREEVEC(sizes);
	k_read_ascii_table(&lambda,&ndim,&sizes,lambdafile);
	if(sizes[0]!=nl+1) {
		fprintf(stderr,"vmatrix and lambda files incompatible (nl==%d vs sizes[0]=%d).\n",nl,sizes[0]);
		exit(1);
	} /* end if */
	FREEVEC(sizes);

	/* load in the filters */
	k_load_filters(&filter_n,&filter_lambda,&filter_pass,&maxn,&nk,filterfile);

	/* create the rmatrix; this is a big matrix which tabulates the
     projection of each basis element b onto each filter k, as a
     function of redshift; you only have to project onto the filters
     here, since every other spectrum you will project onto the
     filters will be a linear combination of the basis elements b; you
     interpolate the rmatrix to get a particular redshift */
	rmatrix=(float *) malloc(nz*nv*nk*sizeof(float));
	zvals=(float *) malloc(nz*sizeof(float));
	for(i=0;i<nz;i++)
		zvals[i]=zmin+(zmax-zmin)*((float)i+0.5)/(float)nz;
	k_projection_table(rmatrix,nk,nv,vmatrix,lambda,nl,zvals,nz,filter_n,
										 filter_lambda,filter_pass,band_shift,maxn);

	/* get the coefficients; sequentially calls k_fit_nonneg to get the
     coefficients */
	nchunk=20;
	redshift=(float *) malloc((nchunk+1)*sizeof(float));
	maggies=(float *) malloc((nchunk+1)*nk*sizeof(float));
	maggies_ivar=(float *) malloc((nchunk+1)*nk*sizeof(float));
	coeffs=(float *) malloc((nchunk+1)*nv*sizeof(float));
	chi2=(float *) malloc((nchunk+1)*sizeof(float));
	fscanf(stdin,"%f",&(maggies[0]));
	while(!feof(stdin)) {
		for(i=0;i<nchunk && !feof(stdin);i++) {
			for(k=1;k<nk;k++)
				fscanf(stdin,"%f",&(maggies[i*nk+k]));
			for(k=0;k<nk;k++)
				fscanf(stdin,"%f",&(maggies_ivar[i*nk+k]));
			fscanf(stdin,"%f",&(maggies[(i+1)*nk+0]));
		} /* end for i */
		ncurrchunk=i;
		/* no direct constraints on the coeffs are included in this fit */
		k_fit_photoz(redshift,coeffs,rmatrix,nk,nv,zprior,lprior,nprior,zvals, 
                 nz,maggies,maggies_ivar,ncurrchunk,tolerance,maxiter,&niter,
                 chi2,0);
		for(i=0;i<ncurrchunk;i++) {
			fprintf(stdout,"%e ",redshift[i]);
			for(j=0;j<nv;j++)
				fprintf(stdout,"%e ",coeffs[i*nv+j]);
			fprintf(stdout,"\n");
		} /* end for i */
		maggies[0]=maggies[ncurrchunk*nk+0];
	}

	return(0);
} /* end main */
