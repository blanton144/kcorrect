;+
; NAME:
;   sdss_to_maggies
; PURPOSE:
;   convert SDSS data to AB, Galactic extinction corrected maggies
; CALLING SEQUENCE:
;   sdss_to_maggies, maggies, ivar [, tsobj=, calibobj=, flux= ]
; INPUTS:
;   calibobj - [N] photoop-style SDSS structure, containing:
;                  .PETROFLUX[5]
;                  .PETROFLUX_IVAR[5]
;                  .MODELFLUX[5]
;                  .MODELFLUX_IVAR[5]
;                  .PSFFLUX[5]
;                  .PSFFLUX_IVAR[5]
;                  .EXTINCTION[5]
;   tsobj - [N] opdb-style SDSS structure, containing:
;                  .PETROCOUNTS[5]
;                  .PETROCOUNTSERR[5]
;                  .COUNTS_MODEL[5]
;                  .COUNTS_MODELERR[5]
;                  .PSFCOUNTS[5]
;                  .PSFCOUNTSERR[5]
;                  .REDDENINg[5]
;   cas - [N] CAS-style structure, with:
;                  .MODELMAG_[X]  (for X = U,G,R,I, and Z)
;                  .MODELMAGERR[X]  (for X = U,G,R,I, and Z)
;                  .PSFMAG_[X]  (for X = U,G,R,I, and Z)
;                  .PSFMAGERR[X]  (for X = U,G,R,I, and Z)
;                  .PETROMAG_[X]  (for X = U,G,R,I, and Z)
;                  .PETROMAGERR[X]  (for X = U,G,R,I, and Z)
;                  .EXTINCTION_[X]  (for X = U,G,R,I, and Z)
; OPTIONAL INPUTS:
;   flux - use this version of the SDSS fluxes ('PETRO', 'MODEL', or 'PSF')
;          [defaults to 'PETRO'] if tsobj or calibobj keywords are
;          used 
; OUTPUTS:
;   maggies - [5, N] output in AB maggies in ugriz
;   ivar - [5, N] inverse variance of maggies
; COMMENTS:
;   You must specify calibobj OR tsobj.
;
;   This code ALWAYS assumes you are inputting pipeline SDSS quantities and
;   converts them to AB quantities!
;
;   It ALWAYS applies a minimum error of [0.05, 0.02, 0.02, 0.02, 0.03] 
;   in ugriz respectively.
;
;   calibobj will be interpreted as a photoop-style SDSS structure.
;
;   tsobj will be interpreted as a opdb-style SDSS structure.
;
;   In either the case of calibobj or tsobj, the Petrosian flux will
;   be used by default --- or you can specify flux="model" or
;   flux="psf".
;
;   In addition, in the case of calibobj or tsobj, this code uses the
;   Galactic extinction values in those structures to correct the
;   fluxes. 
; REVISION HISTORY:
;   07-Apr-2005  Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro sdss_to_maggies, maggies, ivar, tsobj=tsobj, calibobj=calibobj, flux=flux, $
                     cas=cas

if(NOT keyword_set(flux)) then flux='petro'

if(long(n_tags(tsobj) gt 0)+long(n_tags(calibobj) gt 0)+long(n_tags(cas) gt 0) ne 1) then begin
    splog, 'One and only one of tsobj OR calibobj must be set.'
    maggies=-1.
    ivar=-1.
    return
endif

;; if we have input a tsobj structure, interpret it 
if(n_tags(tsobj) gt 0) then begin
    if(flux ne 'model') then begin
        magname=flux+'counts'
        errname=magname+'err'
    endif else begin
        magname='counts_model'
        errname='counts_modelerr'
    endelse 
    imag=tag_indx(tsobj[0], magname)
    if(imag eq -1) then begin
        splog, 'ERROR: tsobj MUST have .'+magname
        return
    endif
    ierr=tag_indx(tsobj[0], errname)
    if(ierr eq -1) then begin
        splog, 'ERROR: tsobj MUST have .'+errname
        return
    endif
    if(n_elements(tsobj[0].(imag)) ne 5) then begin
        splog, 'ERROR: tsobj.'+magname+' MUST have five elements!'
        return
    endif
    if(n_elements(tsobj[0].(ierr)) ne 5) then begin
        splog, 'ERROR: tsobj.'+errname+' MUST have five elements!'
        return
    endif
    sdss_mag=tsobj.(imag)
    sdss_err=tsobj.(ierr)
    ireddening=tag_indx(tsobj[0], 'reddening')
    if(ireddening eq -1) then begin
        splog, 'ERROR: tsobj MUST have .reddening'
        return
    endif
    if(n_elements(tsobj[0].(ireddening)) ne 5) then begin
        splog, 'ERROR: tsobj.reddening MUST have five elements!'
        return
    endif
    sdss_mag=sdss_mag-tsobj.reddening
    k_sdssfix, sdss_mag, sdss_err, maggies, ivar 
    return
endif

;; if it is CAS output
if(n_tags(cas) gt 0) then begin
    filters=['u', 'g', 'r', 'i', 'z']
    sdss_mag=fltarr(n_elements(filters), n_elements(cas))
    sdss_err=fltarr(n_elements(filters), n_elements(cas))
    for i=0L, n_elements(filters)-1L do begin
        magname=flux+'mag_'+filters[i]
        errname=flux+'magerr_'+filters[i]
        imag=tag_indx(cas[0], magname)
        if(imag eq -1) then begin
            splog, 'ERROR: cas MUST have .'+magname
            return
        endif
        ierr=tag_indx(cas[0], errname)
        if(ierr eq -1) then begin
            splog, 'ERROR: cas MUST have .'+errname
            return
        endif
        sdss_mag[i,*]=cas.(imag)
        sdss_err[i,*]=cas.(ierr)
        iext=tag_indx(cas[0], 'extinction_'+filters[i])
        if(iext eq -1) then begin
            splog, 'ERROR: cas MUST have .extinction_'+filters[i]
            return
        endif
        sdss_mag[i,*]=sdss_mag[i,*]-cas.(iext)
        k_sdssfix, sdss_mag, sdss_err, maggies, ivar 
        return
    endfor
endif

;; if we have input a calibobj structure, interpret it 
if(n_tags(calibobj) gt 0) then begin
    fluxname=flux+'flux'
    ivarname=fluxname+'_ivar'
    iflux=tag_indx(calibobj[0], fluxname)
    if(iflux eq -1) then begin
        splog, 'ERROR: calibobj MUST have .'+fluxname
        return
    endif
    iivar=tag_indx(calibobj[0], ivarname)
    if(iivar eq -1) then begin
        splog, 'ERROR: calibobj MUST have .'+ivarname
        return
    endif
    if(n_elements(calibobj[0].(iflux)) ne 5) then begin
        splog, 'ERROR: calibobj.'+fluxname+' MUST have five elements!'
        return
    endif
    if(n_elements(calibobj[0].(iivar)) ne 5) then begin
        splog, 'ERROR: calibobj.'+ivarname+' MUST have five elements!'
        return
    endif
    sdss_nmgy=calibobj.(iflux)
    sdss_ivar=calibobj.(iivar)
    iextinction=tag_indx(calibobj[0], 'extinction')
    if(iextinction eq -1) then begin
        red_fac = [5.155, 3.793, 2.751, 2.086, 1.479 ]
        euler,calibobj.ra,calibobj.dec,ll,bb,1
        extinction= red_fac # dust_getval(ll, bb, /interp, /noloop)
    endif else begin
        extinction=calibobj.extinction
    endelse
    if(n_elements(extinction[*,0]) ne 5) then begin
        splog, 'ERROR: calibobj.extinction MUST have five elements!'
        return
    endif
    sdss_nmgy=sdss_nmgy*10.^(0.4*extinction)
    sdss_ivar=sdss_ivar*10.^(-0.8*extinction)
    maggies=sdss_nmgy*1.e-9
    ivar=sdss_ivar*1.e+18
    k_abfix, maggies, ivar
    k_minerror, maggies, ivar
    return
endif

end
