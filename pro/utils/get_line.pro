;+
; NAME:
;   get_line
; PURPOSE:
;   measures a line flux and equivalent width
; CALLING SEQUENCE:
;   lstruct= get_line(loglam, flux, flux_ivar [, blim=, rlim=, llim=, lname=]
; INPUTS:
;   loglam - log base 10 rest-frame wavelength (angstroms)
;   flux - flux in flambda units
;   flux_ivar - inverse variance of flux
; OPTIONAL INPUTS
;   lname - name of line (default to 'LINE')
;   blim - [2] limits in Angstroms of blue continuum defn
;   rlim - [2] limits in Angstroms of red continuum defn
;   llim - [2] limits to integrate line over
; OUTPUTS:
;   lstruct - structure with the elements:
;      .[LNAME]_FLUX - flux in line (units are flambda*Angstroms)
;      .[LNAME]_FLUX_IVAR - flux in line (units are flambda*Angstroms)
;      .[LNAME]_CONTINUUM - flux density in continuum
;      .[LNAME]_EQW - equivalent widths
; COMMENTS:
;   If "lname" is known, uses a standard set of rlim, blim, and
;   llmin. Otherwise those must be set.
;
;   Known lnames are:
;       HDELTA
; 
;   Only *really* makes sense to measure absorption lines with such a
;   crude measure of the continnum.
; REVISION HISTORY:
;   05-May-2005  Written by M. Blanton, NYU
;-  
;------------------------------------------------------------------------------
function get_line, loglam, flux, flux_ivar, blim=blim, rlim=rlim, llim=llim, $
                   lname=lname

if(n_elements(flux_ivar) eq 0) then $
  flux_ivar=fltarr(n_elements(flux))+1.

if(keyword_set(lname)) then begin
    if(lname eq 'HDELTA') then begin
        blim=alog10([4041., 4080.])
        rlim=alog10([4128., 4161.])
        llim=alog10([4083., 4122.])
    endif
endif else begin
    lname='LINE'
endelse

iblue=where(loglam gt blim[0] and loglam lt blim[1], nblue)
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

ired=where(loglam gt rlim[0] and loglam lt rlim[1], nred)
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

continuum=0.5*(redc+bluec)
sig2continuum=0.5*(sig2redc+sig2bluec)

iline=where(loglam gt llim[0] and loglam lt llim[1], nline)
if(nline eq 0) then begin
    ivar=0.
    return, 1.
endif
sweight=total(flux_ivar[iline] gt 0.)
if(sweight le 0.) then begin
    ivar=0.
    return, 1.
endif
dlam=0.5*(10.^(loglam[iline+1L])-10.^(loglam[iline-1L]))
lflux=total((flux[iline]-continuum[0])*float(flux_ivar[iline] gt 0.)*dlam)
sig2lflux=total(dlam^2*float(flux_ivar[iline] gt 0.)/ $
                (flux_ivar[iline]+1.*float(flux_ivar[iline] le 0.)))
lflux_ivar=1./sig2lflux
continuum_ivar=1./sig2continuum
eqw=lflux/continuum
eqw_ivar=(1./eqw^2)/(1./(lflux^2*lflux_ivar)+1./(continuum^2*continuum_ivar))

lstruct=create_struct(lname+'_FLUX', lflux, $
                      lname+'_FLUX_IVAR', lflux_ivar, $
                      lname+'_CONTINUUM', continuum, $
                      lname+'_CONTINUUM_IVAR', continuum_ivar, $
                      lname+'_EQW', eqw, $
                      lname+'_EQW_IVAR', eqw_ivar)
return, lstruct

end
