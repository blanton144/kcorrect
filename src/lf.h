#include "ztransform.h"
#ifndef IDL_LONG_DEFINED
typedef int IDL_LONG;
#define IDL_LONG_DEFINED
#endif
float *lf_vector(long nl, long nh);
void lf_free_vector(float *v, long nl, long nh);
float lf_zbrent(float (*func)(float), float x1, float x2, float tol);
void lf_polint(float xa[], float ya[], IDL_LONG n, float x, float *y, 
               float *dy);
int lf_calc_vmax(float in_absm, float *in_coeffs, int in_nv, float *in_zvals, 
                 int in_nz, float *in_rmatrix, int in_nk,
                 float in_sample_zmin, float in_sample_zmax, float mmin, 
                 float mmax, float in_qevolve, float in_qz0, 
                 float in_absmagdep, float in_ref_absmagdep, 
                 float in_band_shift, float in_magoffset, float in_omega0, 
                 float in_omegal0, float *zmin, float *zmax);
IDL_LONG lf_eep(float *in_redshift, float *in_absmag, 
                float *in_absmmin, float *in_absmmax, IDL_LONG in_ngals,
                float in_qevolve, float in_qz0, float in_absmagdep,
                float in_absmagdep0, float in_sample_absmmin,
                float in_sample_absmmax, float *in_absmk, float *in_phi,
                float *in_phi_err, float *in_covar, IDL_LONG in_nbin, 
                IDL_LONG calc_err, float *in_weight);
IDL_LONG lf_select_eep(float *in_redshift, float *in_absmag, 
                       float *in_absmmin, float *in_absmmax, 
                       IDL_LONG in_ngals, float in_qevolve, float in_qz0,
                       float in_absmagdep, float in_absmagdep0,
                       float sample_absmmin, float sample_absmmax, 
                       float *in_absmk, float *in_phi, float *sel,
                       IDL_LONG in_nbin);
