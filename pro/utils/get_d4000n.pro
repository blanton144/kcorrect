;+
; NAME:
;   get_d4000n
; PURPOSE:
;   get "narrow" D4000 measure from a spectrum
; CALLING SEQUENCE:
;   d4000n= get_d4000n(loglam, flux, flux_ivar, ivar= )
; INPUTS:
;   loglam - log base 10 rest-frame wavelength (angstroms)
;   flux - flux in flambda units
;   flux_ivar - inverse variance of flux
; OUTPUTS:
;   d4000n - ratio of red to blue continuum
;   ivar - inverse variance of d4000n
; REVISION HISTORY:
;   05-May-2005  Written by M. Blanton, NYU
;-  
;------------------------------------------------------------------------------
function get_d4000n, loglam, flux, flux_ivar, ivar=ivar

if(n_elements(flux_ivar) eq 0) then $
  flux_ivar=fltarr(n_elements(flux))+1.

bluelim=alog10([3850., 3950.])
redlim=alog10([4000., 4100.])

iblue=where(loglam gt bluelim[0] and loglam lt bluelim[1], nblue)
if(nblue eq 0) then begin
    ivar=0.
    return, 1.
endif
sweight=total(flux_ivar[iblue])
if(sweight le 0.) then begin
    ivar=0.
    return, 1.
endif
sflux=total(flux[iblue]*flux_ivar[iblue])
bluec=sflux/sweight
if(bluec le 0.) then begin
    ivar=0. 
    return, 1.e+5
endif
sig2bluec=1./sweight

ired=where(loglam gt redlim[0] and loglam lt redlim[1], nred)
if(nred eq 0) then begin
    ivar=0.
    return, 1.
endif
sweight=total(flux_ivar[ired])
if(sweight le 0.) then begin
    ivar=0.
    return, 1.
endif
sflux=total(flux[ired]*flux_ivar[ired])
redc=sflux/sweight
sig2redc=1./sweight

d4000n=redc/bluec
ivar=(1./d4000n^2)/(sig2redc/redc^2+sig2bluec/bluec^2)
return, d4000n

end
