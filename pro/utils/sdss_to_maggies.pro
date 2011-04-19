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
;   be used by default --- or you can specify flux="model",
;   flux="psf", or flux="cmodel".
;
;   In addition, in the case of calibobj or tsobj, this code uses the
;   Galactic extinction values in those structures to correct the
;   fluxes. 
; REVISION HISTORY:
;   07-Apr-2005  Mike Blanton, NYU
;   09-Feb-2011  J. Moustakas, UCSD - deal with cmodel magnitudes
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
     endfor
    k_sdssfix, sdss_mag, sdss_err, maggies, ivar 
    return
endif

;; if we have input a calibobj structure, interpret it 
if(n_tags(calibobj) gt 0) then begin
; cmodel maggies are a special case
    if strtrim(strlowcase(flux),2) eq 'cmodel' then begin
       ifracpsf = tag_indx(calibobj[0],'fracpsf')
       idevflux = tag_indx(calibobj[0],'devflux')
       iexpflux = tag_indx(calibobj[0],'expflux')
       idevflux_ivar = tag_indx(calibobj[0],'devflux_ivar')
       iexpflux_ivar = tag_indx(calibobj[0],'expflux_ivar')
; check for everything we need
       if (ifracpsf eq -1) then begin
          splog, 'ERROR: calibobj MUST have .FRACPSF'
          return
       endif
       if (idevflux eq -1) then begin
          splog, 'ERROR: calibobj MUST have .DEVFLUX'
          return
       endif
       if (iexpflux eq -1) then begin
          splog, 'ERROR: calibobj MUST have .EXPFLUX'
          return
       endif
       if (idevflux_ivar eq -1) then begin
          splog, 'ERROR: calibobj MUST have .DEVFLUX_IVAR'
          return
       endif
       if (iexpflux_ivar eq -1) then begin
          splog, 'ERROR: calibobj MUST have .EXPFLUX_IVAR'
          return
       endif
; check for the right number of elements
       if (n_elements(calibobj[0].(ifracpsf)) ne 5) then begin
          splog, 'ERROR: calibobj.FRACPSF MUST have five elements!'
          return
       endif
       if (n_elements(calibobj[0].(idevflux)) ne 5) then begin
          splog, 'ERROR: calibobj.DEVFLUX MUST have five elements!'
          return
       endif
       if (n_elements(calibobj[0].(iexpflux)) ne 5) then begin
          splog, 'ERROR: calibobj.EXPFLUX MUST have five elements!'
          return
       endif
       if (n_elements(calibobj[0].(idevflux_ivar)) ne 5) then begin
          splog, 'ERROR: calibobj.DEVFLUX_IVAR MUST have five elements!'
          return
       endif
       if (n_elements(calibobj[0].(iexpflux_ivar)) ne 5) then begin
          splog, 'ERROR: calibobj.EXPFLUX_IVAR MUST have five elements!'
          return
       endif
; construct the composite flux according to
; http://www.sdss3.org/dr8/algorithms/magnitudes.php#cmodel 
       sdss_nmgy = fltarr(5,n_elements(calibobj))
       sdss_ivar = sdss_nmgy
       for ii = 0, 4 do begin
          good = where((calibobj.(ifracpsf))[ii,*] ge 0.0,ngood)
          if (ngood ne 0L) then begin
             frac = (calibobj[good].(ifracpsf))[ii,*]
             sdss_nmgy[ii,good] = frac*(calibobj[good].(idevflux))[ii,*] + $
               (1-frac)*(calibobj[good].(iexpflux))[ii,*]
             sdss_ivar[ii,good] = (calibobj[good].(idevflux_ivar))[ii,*]/(frac^2+(frac eq 0))*(frac ne 0) + $
               (calibobj[good].(iexpflux_ivar))[ii,*]/((1-frac)^2+((1-frac) eq 0))*((1-frac) ne 0)
          endif       
       endfor
    endif else begin
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
    endelse
    iextinction=tag_indx(calibobj[0], 'extinction')
    if(iextinction eq -1) then begin
        red_fac = [5.155, 3.793, 2.751, 2.086, 1.479 ]
; check for an alternate ra,dec tag
        ratag = tag_indx(calibobj,'ra')
        dectag = tag_indx(calibobj,'dec')
        if (ratag eq -1) or (dectag eq -1) then begin
           ratag = tag_indx(calibobj,'racen')
           dectag = tag_indx(calibobj,'deccen')
           if (ratag eq -1) or (dectag eq -1) then begin
              splog, 'RA,DEC tags missing!'
              return
           endif
        endif
        euler,calibobj.(ratag),calibobj.(dectag),ll,bb,1
        extinction= red_fac # dust_getval(ll, bb, /interp, /noloop)
    endif else begin
        extinction=calibobj.extinction
    endelse
    if(n_elements(extinction[*,0]) ne 5) then begin
        splog, 'ERROR: calibobj.extinction MUST have five elements!'
        return
    endif
    sdss_nmgy=sdss_nmgy*10.D^(0.4*extinction)
    sdss_ivar=sdss_ivar*10.D^(-0.8*extinction)
    maggies=sdss_nmgy*1.e-9
    ivar=sdss_ivar*1.e+18
    k_abfix, maggies, ivar
    k_minerror, maggies, ivar
    maggies = float(maggies)
    ivar = float(ivar)
    return
endif

end
