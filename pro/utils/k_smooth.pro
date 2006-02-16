;+
; NAME:
;   k_smooth
; PURPOSE:
;   gaussian smooth a spectrum
; CALLING SEQUENCE:
;   outflux= k_smooth(loglam, flux, vdisp)
; INPUTS:
;   loglam - log_{10} of the wavelength in Angstroms
;   flux - input flux
;   vdisp - gaussian velocity width in km/s (ignores if < 10 km/s)
; OUTPUTS:
;   outflux - smoothed flux
; WARNINGS:
;   Does NOTHING if vdisp < 1 km/s
; REVISION HISTORY:
;   05-May-2005  Michael Blanton (NYU)
;-
;------------------------------------------------------------------------------
function k_smooth, loglam,flux,vdisp

if vdisp GT 1.0 then begin
    nlambda= n_elements(loglam)
    pixsize= $
      abs(alog(10.)*2.99792e+5*(loglam[nlambda-1]-loglam[0])/double(nlambda))
    smoothing= vdisp/pixsize    ; pixels
    npix= long(4.0*ceil(smoothing))*2L+3
    klam= findgen(npix)-float(npix-1.)/2.
    kernel= exp(-0.5*(klam/smoothing)^2)/sqrt(2.*!DPI)/smoothing
    kernel= kernel/total(kernel)
    smoothed_spec= convol(flux,kernel,/edge_truncate)
endif else begin
    smoothed_spec= flux
endelse

return, smoothed_spec

end
