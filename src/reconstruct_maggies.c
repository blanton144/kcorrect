#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <kcorrect.h>

/*
 * reconstruct_maggies.c
 *
 * given coefficients, a redshift, a filterlist, and the templates,
 * reconstruct the maggies at particular redshift
 *
 * Mike Blanton
 * 7/2003
 */

static IDL_LONG nz=1000;
static float zmin=1.e-4, zmax=1.e-0;
static float band_shift=0.;
static float all_redshift=-1.;

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
static float *coeffs=NULL;

#define USAGE \
   { fprintf(stderr,"Usage: cat <coeffs file> | reconstruct_maggies [--vfile <vfile> --lfile <lfile>\n"); \
     fprintf(stderr,"            --ffile <ffile> --band-shift <band_shift> --redshift <redshift>]\n"); }

int main(int argc,
				 char **argv)
{
	IDL_LONG i,j,k,c,ndim,nchunk,ncurrchunk,*sizes=NULL;
	char vfile[2000],lfile[2000],ffile[2000],path[2000];
	char vmatrixfile[2000],lambdafile[2000],filterfile[2000];

	/* read arguments */
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
				{"band-shift", 1, 0, 0},
				{"redshift", 1, 0, 0},
				{"help", 0, 0, 0}
			};
		static const char short_options[]="v:l:p:f:b:r:h";
		static const char short_options_c[]="vlpfbrh";

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
		case 'b':
			band_shift=atof(optarg);
			break; 
		case 'r':
			all_redshift=atof(optarg);
			break; 
		case 'h':
      USAGE;
      exit(1);
			break; 
		case '?':
			break;
		default: 
			printf("fit_coeffs: getopt returned character code 0%o ??\n", c);
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

	/* reconstruct the magggies */
	nchunk=20;
	redshift=(float *) malloc((nchunk+1)*sizeof(float));
	maggies=(float *) malloc((nchunk+1)*nk*sizeof(float));
	coeffs=(float *) malloc((nchunk+1)*nv*sizeof(float));
	fscanf(stdin,"%f",&(redshift[0]));
	while(!feof(stdin)) {
		for(i=0;i<nchunk && !feof(stdin);i++) {
      if(all_redshift!=-1.) redshift[i]=all_redshift;
			for(j=0;j<nv;j++)
				fscanf(stdin,"%f",&(coeffs[i*nv+j]));
			fscanf(stdin,"%f",&(redshift[i+1]));
		} /* end for i */
		ncurrchunk=i;
		/* no direct constraints on the coeffs are included in this fit */
		k_reconstruct_maggies(zvals,nz,rmatrix,nk,nv,coeffs,redshift,maggies, 
                          ncurrchunk);
		for(i=0;i<ncurrchunk;i++) {
			fprintf(stdout,"%e ",redshift[i]);
			for(k=0;k<nk;k++)
				fprintf(stdout,"%e ",maggies[i*nk+k]);
			fprintf(stdout,"\n");
		} /* end for i */
		redshift[0]=redshift[ncurrchunk];
	}

	return(0);
} /* end main */
