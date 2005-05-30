;+
; NAME:
;   k_sdss_tests_generic
; PURPOSE:
;   runs tests on SDSS test data
; CALLING SEQUENCE:
;   k_sdss_tests
; COMMENTS:
;   Reads in spobj_test.fits and runs:
;      sdss_kcorrect to get corrections, evaluates chi2
;      does own corrections, runs kcorrect, evaluates chi2
;      compares these two
;   Reads in obj_test.fits and runs:
;      sdss_kcorrect to get corrections, evaluates chi2
;      does own corrections, runs kcorrect, evaluates chi2
;      compares these two
;   Compares all four of the above results.
; DATA DEPENDENCIES:
;   $KCORRECT_DIR/data/test/obj_test.fits    (calibObj version)
;   $KCORRECT_DIR/data/test/spobj_test.fits  (tsObj version)
; REVISION HISTORY:
;   2005-04-07 MRB, NYU
;-
;------------------------------------------------------------------------------
pro k_sdss_tests_generic, vname=vname

k_print, filename='k_sdss_tests_pvf.ps', pold=pold, xold=xold, yold=yold

obj=mrdfits(getenv('KCORRECT_DIR')+'/data/test/obj_test.fits',1)
spspec=mrdfits(getenv('KCORRECT_DIR')+'/data/test/spobj_test.fits',1)

kc_spall1=sdss_kcorrect(spspec.zfinal, calibobj=obj, vname=vname)
petroflux=obj.petroflux*10.^(0.4*obj.extinction)
petroflux_ivar=obj.petroflux_ivar*10.^(-0.8*obj.extinction)
kcorrect, petroflux, petroflux_ivar, spspec.zfinal, kc_spall2, $
  /abfix, minerrors=[0.05,0.02,0.02,0.02,0.03], vname=vname
for i=0L, 4L do $
  djs_plot, spspec.zfinal, kc_spall1[i,*]-kc_spall2[i,*], psym=4, $
  xti='z', yti='kc1-kc2 (spall, band='+strtrim(string(i),2)+')'

kc_spspec1=sdss_kcorrect(spspec.zfinal, tsobj=spspec, vname=vname)
petrocounts=spspec.petrocounts-spspec.reddening
kcorrect, petrocounts, spspec.petrocountserr, spspec.zfinal, kc_spspec2, $
  /sdssfix, vname=vname
for i=0L, 4L do $
  djs_plot, spspec.zfinal, kc_spspec1[i,*]-kc_spspec2[i,*], psym=4, $
  xti='z', yti='kc1-kc2 (spspec, band='+strtrim(string(i),2)+')'

!P.MULTI=[0,1,2]
for i=0L, 4L do begin
    djs_plot, spspec.zfinal, kc_spall1[i,*], psym=4, $
      xti='z', yti='kc (band='+strtrim(string(i),2)+')', xra=[0.01,0.4]
    djs_oplot, spspec.zfinal, kc_spspec1[i,*], psym=4, $
      xti='z', yti='kc (band='+strtrim(string(i),2)+')', color='red', $
      xra=[0.01, 0.4]
    djs_plot, spspec.zfinal, kc_spall1[i,*]-kc_spspec1[i,*], psym=4, $
      xti='z', yti='kcall-kcspec (band='+strtrim(string(i),2)+')', $
      xra=[0.01, 0.4]
endfor

k_end_print, xold=xold, yold=yold, pold=pold

end
