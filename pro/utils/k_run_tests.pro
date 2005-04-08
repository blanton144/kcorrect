;+
; NAME:
;   k_run_tests
; PURPOSE:
;   runs tests on test data
; CALLING SEQUENCE:
;   k_run_tests
; COMMENTS:
;   Reads in spspec_test.fits and runs:
;      sdss_kcorrect to get corrections, evaluates chi2
;      does own corrections, runs kcorrect, evaluates chi2
;      compares these two
;   Reads in spall_test.fits and runs:
;      sdss_kcorrect to get corrections, evaluates chi2
;      does own corrections, runs kcorrect, evaluates chi2
;      compares these two
;   Compares all four of the above results.
;   Reads in galex_test.fits and runs:
;      galex_kcorrect, evaluates chi2
;      galex_kcorrect, without SDSS bands, evaluates chi2
;      compares these two
;   Reads in twomass_test.fits and runs:
;      twomass_kcorrect, evaluates chi2
;      twomass_kcorrect, without SDSS bands, evaluates chi2
;      compares these two
;   Reads in deeptest.fits and runs:
;      deep_kcorrect, evaluates chi2
;   Requires that:
;     all the chi2 be good
;     all the comparisons be close
; BUGS:
;   Not finished yet.
; REVISION HISTORY:
;   2005-04-07 MRB, NYU
;-
;------------------------------------------------------------------------------
pro k_run_tests

k_print, filename='k_run_tests.ps'

spspec=mrdfits(getenv('KCORRECT_DIR')+'/data/test/spspec_test.fits',1)
kc_spspec1=sdss_kcorrect(spspec.z, tsobj=spspec)
petrocounts=spspec.petrocounts-spspec.reddening
kcorrect, spspec.z, petrocounts, spspec.petrocountserr, kc_spspec2, $
  /sdssfix
for i=0L, 4L do $
  djs_plot, spspec.z, kc_spspec1[i,*]-kc_spspec2[i,*], psym=4, $
  xra='z', yra='kc1-kc2 (spspec, band='+strtrim(string(i),2)+')'

spall=mrdfits(getenv('KCORRECT_DIR')+'/data/test/spall_test.fits',1)
kc_spall1=sdss_kcorrect(spall.z, tsobj=spall)
petroflux=spall.petroflux*10.^(0.4*spspec.extinction)
petroflux_ivar=spall.petroflux_ivar*10.^(-0.8*spspec.extinction)
kcorrect, spall.z, petroflux, petroflux_ivar, kc_spall2, $
  /abfix, minerrors=[0.05,0.02,0.02,0.02,0.03]
for i=0L, 4L do $
  djs_plot, spspec.z, kc_spall1[i,*]-kc_spall2[i,*], psym=4, $
  xra='z', yra='kc1-kc2 (spall, band='+strtrim(string(i),2)+')'

for i=0L, 4L do $
  djs_plot, spspec.z, kc_spall1[i,*]-kc_spspec1[i,*], psym=4, $
  xra='z', yra='kcall-kcspec (band='+strtrim(string(i),2)+')'

k_end_print


end
