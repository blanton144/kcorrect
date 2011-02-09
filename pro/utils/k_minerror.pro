;+
; NAME:
;   k_minerror
; PURPOSE:
;   apply a minimum magnitude to error to a set of inverse variances
; CALLING SEQUENCE:
;   k_minerror, maggies, maggies_ivar [, minerrors]
; INPUTS:
;   maggies - [nk, n] maggies
; INPUTS/OUTPUTS:
;   maggies_ivar - [nk, n] inverse variances (changed on output)
; OPTIONAL INPUTS:
;   minerrors - [nk] minimum errors to apply (default [0.05, 0.02,
;               0.02, 0.02, 0.03] to represent calibration
;               uncertainties in SDSS ugriz bands)
; COMMENTS:
;   Adds minerror in quadrature to each band.
; BUGS:
;   check "factor"
; REVISION HISTORY:
;   07-Feb-2002  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_minerror, maggies, maggies_ivar, minerrors

if(n_elements(minerrors) eq 0) then $
  minerrors=[0.05, 0.02, 0.02, 0.02, 0.03]
nk=n_elements(minerrors)

for k=0L, nk-1L do begin
    igood=where(maggies[k,*] gt 0 and maggies_ivar[k,*] gt 0, ngood)
    if(ngood gt 0) then begin
        factor=(2.5/alog(10.))
        err=factor/sqrt(maggies_ivar[k,igood])/maggies[k,igood]
        err2=err^2+minerrors[k]^2
        maggies_ivar[k,igood]=factor^2/(maggies[k,igood]^2*err2)
    endif
endfor

end
;------------------------------------------------------------------------------
