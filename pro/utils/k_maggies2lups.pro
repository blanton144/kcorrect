;+
; NAME:
;   k_maggies2lups
; PURPOSE:
;   Convert maggies to SDSS luptitudes ("asinh" magnitudes) 
; CALLING_SEQUENCE:
;   lups= k_maggies2lups(maggies [, maggies_err, lups_err=, bvalues=])
; INPUTS:
;   maggies - [Nb,N] maggies (10.^(-0.4*magnitude)). Unless you set
;             bvalues yourself, Nb=5
; OUTPUTS:
;   lups - [Nb, N] luptitudes (the "asinh" magnitudes used by the
;          SDSS databases). 
; OPTIONAL INPUTS:
;   bvalues - "b" value for luptitude definition for each band
;             (default [1.4D-10, 0.9D-10, 1.2D-10, 1.8D-10, 7.4D-10])
;   maggies_err - [Nb, N] 1 sigma errors in maggies
; OPTIONAL OUTPUTS:
;   lups_err - [Nb, N] 1 sigma errors in luptitudes
;              (fails if maggies_err not set)
; COMMENTS:
;   Conversion from maggies to lups is:
;      lups = - 2.5 ln(b) / ln(10) - 2.5 asinh( 0.5 maggies/b) / ln(10)
;   This is linear when maggies ~ b, identical to maggies 
;   for maggies >> b
;
;   Conversion of the errors is:
;      lups_err = 2.5 maggies_err / ( 2 b ln(10) sqrt(1 + (m/2b)^2)) 
; 
;   If you set bvalues yourself, the code assumes Nb=n_elements(bvalues)
;   for the purposes of intepreting the "lups" input.
; REVISION HISTORY:
;   28-Mar-2002  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function k_maggies2lups,maggies,maggies_err, lups_err=lups_err, bvalues=bvalues

if(NOT keyword_set(bvalues)) then $
  bvalues=[1.4D-10, 0.9D-10, 1.2D-10, 1.8D-10, 7.4D-10]

nb=n_elements(bvalues)
nmaggies=long(n_elements(maggies)/nb)
maggies=reform([maggies],nb,nmaggies)
lups=dblarr(nb,nmaggies)
for b=0L, nb-1L do $
  lups[b,*]=2.5D*alog10(1.0D/bvalues[b])- $
      asinh2(0.5D*maggies[b,*]/bvalues[b])/(0.4D*alog(10.D));
if(arg_present(maggies_err)) then begin
    if(n_elements(maggies_err) ne n_elements(maggies)) then $
      message, 'maggies_err must be same size as maggies'
    lups_err=reform([maggies_err],nb,nmaggies)
    lups_err=dblarr(nb,nmaggies)
    for b=0L, nb-1L do $
      lups_err[b,*]=2.5*maggies_err[b,*]/(2.*bvalues[b]*alog(10.)* $
                                          sqrt(1.+(0.5*maggies[b,*]/ $
	bvalues[b])^2))
endif
return,lups

end
