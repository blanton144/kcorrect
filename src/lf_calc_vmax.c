#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <kcorrect.h>
#include "lf.h"

/*
 * lf_calc_vmax.c
 *
 * Calculate vmax for a particular galaxy.
 * Uses K-correct product for K-corrections.
 *
 * Mike Blanton
 * 9/2003 */

#define PI 3.1415926535897932384626433832975

static float absm; 
static float *coeffs=NULL;
static int nv; 
static float *rmatrix=NULL; 
static float *zvals=NULL; 
static int nk;
static int nz; 
static float sample_zmin; 
static float sample_zmax; 
static float mlimit; 
static float qevolve; 
static float qz0; 
static float absmagdep; 
static float ref_absmagdep; 
static float magoffset;
static float omega0; 
static float omegal0; 
static float band_shift; 

float lf_calc_vmax_func(float z)
{
  int ngals;
  float dm,evol,recm,recm0,kcorrect,mag,curr_zdep;
  ngals=1;
  dm=z2dm(z,omega0,omegal0);
  curr_zdep=qevolve;
#if 0
  if(absm<ref_absmagdep)
    curr_zdep=absmagdep;
  else 
    curr_zdep=qevolve;
  curr_zdep=qevolve*(1.+absmagdep*(absm-ref_absmagdep));
#endif
  evol=curr_zdep*(z-qz0);
  k_reconstruct_maggies(zvals,nz,rmatrix,nk,nv,coeffs,&band_shift,&recm,ngals);
  recm/=(1.+band_shift);
  k_reconstruct_maggies(zvals,nz,rmatrix,nk,nv,coeffs,&z,&recm0,ngals);
  kcorrect=2.5*log10(recm/recm0);
  mag=absm+dm+kcorrect-evol-magoffset;
#if 0
  printf("%d %d %d %e %e %e %e %e %e %e %e %e %e %e\n",nz,nk,nv,z,absm,dm,kcorrect,evol,magoffset,recm,recm0,coeffs[0],coeffs[1],coeffs[2]);
#endif
  return(mag-mlimit);
}

int lf_calc_vmax(float in_absm, 
                 float *in_coeffs,
                 int in_nv, 
                 float *in_zvals, 
                 int in_nz, 
                 float *in_rmatrix, 
                 int in_nk,
                 float in_sample_zmin, 
                 float in_sample_zmax, 
                 float mmin, 
                 float mmax,
                 float in_qevolve, 
                 float in_qz0, 
                 float in_absmagdep, 
                 float in_ref_absmagdep, 
                 float in_band_shift,
                 float in_magoffset,
                 float in_omega0, 
                 float in_omegal0, 
                 float *zmin, 
                 float *zmax)
{
  absm=in_absm;
  coeffs=in_coeffs;
  rmatrix=in_rmatrix;
  zvals=in_zvals;
  nv=in_nv;
  nz=in_nz;
  nk=in_nk;
  sample_zmin=in_sample_zmin; 
  sample_zmax=in_sample_zmax; 
  qevolve=in_qevolve; 
  qz0=in_qz0; 
  absmagdep=in_absmagdep; 
  ref_absmagdep=in_ref_absmagdep; 
  band_shift=in_band_shift; 
  magoffset=in_magoffset;
  omega0=in_omega0; 
  omegal0=in_omegal0; 

  mlimit=mmin;
  if(lf_calc_vmax_func(sample_zmin)<0.) {
    if(lf_calc_vmax_func(sample_zmax)>=0.) {
      *zmin=lf_zbrent(lf_calc_vmax_func,sample_zmin,sample_zmax,1.e-5);
    } else {
    *zmin=sample_zmax;
    }
  } else {
    *zmin=sample_zmin;
  }

  mlimit=mmax;
  if(lf_calc_vmax_func(sample_zmin)<=0.) {
    if(lf_calc_vmax_func(sample_zmax)>0.) {
      *zmax=lf_zbrent(lf_calc_vmax_func,sample_zmin,sample_zmax,1.e-5);
    } else {
    *zmax=sample_zmax;
    }
  } else {
    *zmax=sample_zmin;
  }

  return(1);
  
} /* end lf_calc_vmax */
