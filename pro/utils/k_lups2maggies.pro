;+
; NAME:
;   k_lups2maggies
; PURPOSE:
;   Convert SDSS luptitudes ("asinh" magnitudes) to maggies
; CALLING_SEQUENCE:
;   maggies= k_lups2maggies(lups [, lups_err, maggies_err=, bvalues=])
; INPUTS:
;   lups - [Nb, N] luptitudes (the "asinh" magnitudes from the
;          SDSS databases). Unless you set bvalues yourself, Nb=5
; OUTPUTS:
;   maggies - [Nb,N] maggies (10.^(-0.4*magnitude))
; OPTIONAL INPUTS:
;   lups_err - [Nb, N] 1 sigma errors in luptitudes
;   bvalues - "b" value for luptitude definition for each band
;             (default [1.4D-10, 0.9D-10, 1.2D-10, 1.8D-10, 7.4D-10])
; OPTIONAL OUTPUTS:
;   maggies_err - [Nb, N] 1 sigma errors in maggies
;                 (fails if lups_err not set)
; COMMENTS:
;   Conversion from luptitudes to maggies is:
;      maggies = 2 b sinh( - ln(b) - 0.4 ln(10) lups)
;   This is linear when maggies ~ b, identical to maggies 
;   for maggies >> b
;
;   Conversion of the errors is:
;      maggies_err = 2 b cosh( -ln(b) - 0.4 ln(10) lups) (0.4 ln 10) lups_err
; 
;   If you set bvalues yourself, the code assumes Nb=n_elements(bvalues)
;   for the purposes of intepreting the "lups" input.
; REVISION HISTORY:
;   28-Mar-2002  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function k_lups2maggies,lups, lups_err, maggies_err=maggies_err, $
                        bvalues=bvalues

if(NOT keyword_set(bvalues)) then $
  bvalues=[1.4D-10, 0.9D-10, 1.2D-10, 1.8D-10, 7.4D-10]

nb=n_elements(bvalues)
nlups=long(n_elements(lups)/nb)
lups=reform([lups],nb,nlups)
maggies=dblarr(nb,nlups)
for b=0L, nb-1L do $
  maggies[b,*]=2.D*bvalues[b]*sinh(-alog(bvalues[b])-0.4D*alog(10.D)*lups[b,*])
if(arg_present(maggies_err)) then begin
    if(n_elements(lups_err) ne n_elements(lups)) then $
      message, 'lups_err must be same size as lups'
    lups_err=reform([lups_err],nb,nlups)
    maggies_err=dblarr(nb,nlups)
    for b=0L, nb-1L do $
      maggies_err[b,*]=2.D*bvalues[b]*cosh(-alog(bvalues[b])- $
                                           0.4D*alog(10.D)*lups[b,*])* $
      0.4*alog(10.)*lups_err[b,*]
endif
return,maggies

end
