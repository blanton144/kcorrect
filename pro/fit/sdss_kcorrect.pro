;+
; NAME:
;   sdss_kcorrect
; PURPOSE:
;   calculate K-corrections for standard SDSS input
; CALLING SEQUENCE:
;   kcorrect= sdss_kcorrect(redshift [, nmgy=, ivar=, mag=, err=, $
;                           calibobj=, tsobj=, flux=, band_shift=,
;                           chi2=])
; INPUTS:
;   calibobj - photoop-style structure, containing:
;                  .PETROFLUX[5]
;                  .PETROFLUX_IVAR[5]
;                  .MODELFLUX[5]
;                  .MODELFLUX_IVAR[5]
;                  .PSFFLUX[5]
;                  .PSFFLUX_IVAR[5]
;                  .EXTINCTION[5]
;   tsobj - opdb-style structure, containing:
;                  .PETROCOUNTS[5]
;                  .PETROCOUNTSERR[5]
;                  .COUNTS_MODEL[5]
;                  .COUNTS_MODELERR[5]
;                  .PSFCOUNTS[5]
;                  .PSFCOUNTSERR[5]
;                  .REDDENINg[5]
;   nmgy, ivar - nanomaggies, Galactic-reddening corrected, and inverse
;                variance of same
;   mag, err - asinh magnitudes, Galactic-reddening corrected and
;              errors of same
; OPTIONAL INPUTS:
;   flux - use this version of the fluxes ('PETRO', 'MODEL', or 'PSF')
;          [defaults to 'PETRO']
;   band_shift    - blueshift of bandpasses to apply (to get ^{z}b
;                   type bands) [default 0.]
; OUTPUTS:
;   kcorrect   - [5, ngals] K-corrections satisfying
;                   m = M + DM(z) + K(z)
;                based on the best fit sum of templates
; OPTIONAL OUTPUTS:
;   chi2       - chi^2 of fit
; COMMENTS:
;   This is a simple wrapper on kcorrect.pro which is almost always
;   just what you want. It keeps a version of rmatrix and zvals in
;   memory to save time, recalculating them each time you change
;   band_shift.
;
;   You must specify nmgy,ivar OR mag,err OR calibobj OR tsobj.  This
;   code ALWAYS assumes you are inputting pipeline quantities and
;   converts them to AB quantities! It also ALWAYS applies a minimum
;   error of [0.05, 0.02, 0.02, 0.02, 0.03] in ugriz respectively.
;   calibobj will be interpreted as a photoop-style SDSS structure.
;   tsobj will be interpreted as a opdb-style SDSS structure.  In either
;   the case of calibobj or tsobj, the Petrosian flux will be used by
;   default --- or you can specify flux="model" or flux="psf".  In
;   addition, in the case of calibobj or tsobj, this code uses the
;   Galactic extinction values in those structures to correct the fluxes.
; REVISION HISTORY:
;   07-Apr-2005  Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function sdss_kcorrect, redshift, nmgy=nmgy, ivar=ivar, mag=mag, err=err, $
                        calibobj=calibobj, tsobj=tsobj, flux=flux, $
                        band_shift=in_band_shift, chi2=chi2

common com_sdss_kcorrect, rmatrix, zvals, band_shift

if(n_params() lt 1 OR $
   (((keyword_set(nmgy) eq 0 OR keyword_set(ivar) eq 0)) AND $
    ((keyword_set(mag) eq 0 OR keyword_set(err) eq 0)) AND $
    (n_tags(calibobj) eq 0) AND $
    (n_tags(tsobj) eq 0))) $
  then begin
    print, 'Usage: kcorrect= sdss_kcorrect(redshift [, nmgy=, ivar=, mag=, err=, $'
    print, '                   calibobj=, tsobj=, flux=, band_shift=, chi2= ])'
    print, ''
    print, 'You must specify nmgy,ivar OR mag,err OR calibobj OR tsobj.'
    print, 'This code ALWAYS assumes you are inputting pipeline quantities and'
    print, 'converts them to AB quantities! It also ALWAYS applies a minimum'
    print, 'error of [0.05, 0.02, 0.02, 0.02, 0.03] in ugriz respectively.' 
    print, 'calibobj will be interpreted as a photoop-style SDSS structure.'
    print, 'tsobj will be interpreted as a opdb-style SDSS structure.'
    print, 'In either the case of calibobj or tsobj, the Petrosian flux will be'
    print, 'used by default --- or you can specify flux="model" or flux="psf".'
    print, 'In addition, in the case of calibobj or tsobj, this code uses the'
    print, 'Galactic extinction values in those structures to correct the fluxes.'
    return, -1
endif 

if(NOT keyword_set(flux)) then flux='petro'

;; if we have input a tsobj structure, interpret it 
if(n_tags(tsobj) gt 0) then begin
    if(flux ne 'model') then begin
        magname=flux+'counts'
        errname=fluxname+'err'
    endif else begin
        magname='counts_model'
        errname='counts_modelerr'
    endelse 
    imag=tag_indx(tsobj[0], magname)
    if(imag eq -1) then begin
        splog, 'ERROR: tsobj MUST have .'+magname
        return,-1
    endif
    ierr=tag_indx(tsobj[0], errname)
    if(ierr eq -1) then begin
        splog, 'ERROR: tsobj MUST have .'+errname
        return,-1
    endif
    if(n_elements(tsobj[0].(imag)) ne 5) then begin
        splog, 'ERROR: tsobj.'+magname+' MUST have five elements!'
        return,-1
    endif
    if(n_elements(tsobj[0].(ierr)) ne 5) then begin
        splog, 'ERROR: tsobj.'+errname+' MUST have five elements!'
        return,-1
    endif
    mag=tsobj.(imag)
    err=tsobj.(ierr)
    ireddening=tag_indx(tsobj[0], 'reddening')
    if(ireddening eq -1) then begin
        splog, 'ERROR: tsobj MUST have .reddening'
        return,-1
    endif
    if(n_elements(tsobj[0].(ireddening)) ne 5) then begin
        splog, 'ERROR: tsobj.reddening MUST have five elements!'
        return,-1
    endif
    mag=mag-tsobj.reddening
endif

;; if we have input a calibobj structure, interpret it 
if(n_tags(calibobj) gt 0) then begin
    fluxname=flux+'flux'
    ivarname=fluxname+'_ivar'
    iflux=tag_indx(calibobj[0], fluxname)
    if(iflux eq -1) then begin
        splog, 'ERROR: calibobj MUST have .'+fluxname
        return,-1
    endif
    iivar=tag_indx(calibobj[0], ivarname)
    if(iivar eq -1) then begin
        splog, 'ERROR: calibobj MUST have .'+ivarname
        return,-1
    endif
    if(n_elements(calibobj[0].(iflux)) ne 5) then begin
        splog, 'ERROR: calibobj.'+fluxname+' MUST have five elements!'
        return,-1
    endif
    if(n_elements(calibobj[0].(iivar)) ne 5) then begin
        splog, 'ERROR: calibobj.'+ivarname+' MUST have five elements!'
        return,-1
    endif
    nmgy=calibobj.(iflux)
    ivar=calibobj.(iivar)
    iextinction=tag_indx(calibobj[0], 'extinction')
    if(iextinction eq -1) then begin
        splog, 'ERROR: calibobj MUST have .extinction'
        return,-1
    endif
    if(n_elements(calibobj[0].(iextinction)) ne 5) then begin
        splog, 'ERROR: calibobj.extinction MUST have five elements!'
        return,-1
    endif
    nmgy=nmgy*10.^(0.4*calibobj.extinction)
    ivar=ivar*10.^(-0.8*calibobj.extinction)
endif

;; need to reset rmatrix if band_shift changes
if(n_elements(band_shift) ne 0) then begin
    if(band_shift ne in_band_shift) then begin
       rmatrix=0
       zvals=0
    endif
endif

;; interpret band_shift
if(keyword_set(in_band_shift)) then $
  band_shift=in_band_shift $
else $
  band_shift=0.

;; if you are using SDSS Archive values, fix em
if(n_elements(mag) gt 0) then begin
    k_sdssfix, mag, err, mgy, mgy_ivar 
endif else begin
    mgy=nmgy*1.e-9
    mgy_ivar=ivar*1.e+18
    k_minerror, mgy, mgy_ivar
    k_abfix, mgy, mgy_ivar
endelse

;; call kcorrect
kcorrect, mgy, mgy_ivar, redshift, kcorrect, band_shift=band_shift, $
  rmatrix=rmatrix, zvals=zvals

return, kcorrect

end
