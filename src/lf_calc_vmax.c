#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <kcorrect.h>

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
static float q0; 
static float q1; 
static float qz0; 
static float magoffset;
static float omega0; 
static float omegal0; 
static float band_shift; 

float lf_calc_vmax_func(float z)
{
  int ngals;
  float dm,evol,recm,recm0,kcorrect,mag;
  ngals=1;
  dm=z2dm(z,omega0,omegal0);
  
  evol=k_evolve(0., z, q0, q1, qz0);
  k_reconstruct_maggies(zvals,nz,rmatrix,nk,nv,coeffs,&band_shift,&recm,ngals);
  recm/=(1.+band_shift);
  k_reconstruct_maggies(zvals,nz,rmatrix,nk,nv,coeffs,&z,&recm0,ngals);
  kcorrect=2.5*log10(recm/recm0);
  mag=absm+dm+kcorrect-evol-magoffset;
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
                 float in_q0, 
                 float in_q1, 
                 float in_qz0, 
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
  q0=in_q0; 
  q1=in_q1; 
  qz0=in_qz0; 
  band_shift=in_band_shift; 
  magoffset=in_magoffset;
  omega0=in_omega0; 
  omegal0=in_omegal0; 

  mlimit=mmin;
  if(lf_calc_vmax_func(sample_zmin)<0.) {
    if(lf_calc_vmax_func(sample_zmax)>=0.) {
      *zmin=k_zbrent(lf_calc_vmax_func,sample_zmin,sample_zmax,1.e-5);
    } else {
    *zmin=sample_zmax;
    }
  } else {
    *zmin=sample_zmin;
  }

  mlimit=mmax;
  if(lf_calc_vmax_func(sample_zmin)<=0.) {
    if(lf_calc_vmax_func(sample_zmax)>0.) {
      *zmax=k_zbrent(lf_calc_vmax_func,sample_zmin,sample_zmax,1.e-5);
    } else {
    *zmax=sample_zmax;
    }
  } else {
    *zmax=sample_zmin;
  }

  return(1);
  
} /* end lf_calc_vmax */
