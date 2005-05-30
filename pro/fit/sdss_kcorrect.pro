;+
; NAME:
;   sdss_kcorrect
; PURPOSE:
;   calculate K-corrections for standard SDSS input
; CALLING SEQUENCE:
;   kcorrect= sdss_kcorrect(redshift [, nmgy=, ivar=, mag=, err=, $
;                           calibobj=, tsobj=, flux=, band_shift=,$
;                           chi2=, rmaggies=, omaggies=, vname=, $
;                           oivar=, mass=, mtol=, absmag=, amivar=, $
;                           omega0=, omegal0= ])
; INPUTS:
;   redshift - [N] redshifts
;   calibobj - [N] photoop-style structure, containing:
;                  .PETROFLUX[5]
;                  .PETROFLUX_IVAR[5]
;                  .MODELFLUX[5]
;                  .MODELFLUX_IVAR[5]
;                  .PSFFLUX[5]
;                  .PSFFLUX_IVAR[5]
;                  .EXTINCTION[5]
;   tsobj - [N] opdb-style structure, containing:
;                  .PETROCOUNTS[5]
;                  .PETROCOUNTSERR[5]
;                  .COUNTS_MODEL[5]
;                  .COUNTS_MODELERR[5]
;                  .PSFCOUNTS[5]
;                  .PSFCOUNTSERR[5]
;                  .REDDENINg[5]
;   nmgy, ivar - [5, N] nanomaggies, Galactic-reddening corrected, and inverse
;                variance of same
;   mag, err - [5, N] asinh magnitudes, Galactic-reddening corrected and
;              errors of same
; OPTIONAL INPUTS:
;   flux - use this version of the fluxes ('PETRO', 'MODEL', or 'PSF')
;          [defaults to 'PETRO'] if tsobj or calibobj keywords are
;          used 
;   band_shift    - blueshift of bandpasses to apply (to get ^{z}b
;                   type bands) [default 0.]
;   vname - name of fit to use (defaults to 'default')
;   omega0, omegal0 - cosmological parameters for calculating distance
;                     moduli [default 0.3, 0.7]
; OUTPUTS:
;   kcorrect - [5, ngals] K-corrections in ugriz satisfying
;                m = M + DM(z) + K(z)
;              based on the best fit sum of templates
;   mtol - [5, ngals] mass-to-light ratios from model in each band
;   mass - [ngals] total mass from model in each band
;   absmag - [5, ngals] absolute magnitude (for missing data, substitutes
;            model fit)
;   amivar - [5, ngals] inverse variance of absolute magnitude (for
;            missing data = 0)
; OPTIONAL OUTPUTS:
;   coeffs - coefficients of fit
;   chi2 - chi^2 of fit
;   rmaggies - [5, N] reconstructed maggies from the fit (ugriz)
;   omaggies, oivar - [5, N] maggies and inverse variances used for fit
;                           (after extinction, AB correction, etc)  (ugriz)
; COMMENTS:
;   This is a simple wrapper on kcorrect.pro which is almost always
;   just what you want. It keeps a version of rmatrix and zvals in
;   memory to save time, recalculating them each time you change
;   band_shift.
;
;   You must specify nmgy,ivar OR mag,err OR calibobj OR tsobj. If
;   nmgy or mag, make sure they are AB calibrated and Galactic
;   extinction corrected.
;
;   Uses sdss_to_maggies to convert tsobj or calibobj structure to
;   AB, Galactic extinction corrected maggies. Passes optional
;   argument "flux" to sdss_to_maggies.
;
;   For v4_0b templates and later, coefficients are in units of:
;     1 solar mass / (D/10pc)^2
;   That is, sum the coefficients and multiply by (D/10pc)^2 to get
;   masses. (In fact, for Omega0=0.3 and OmegaL0=0.7, this is what the
;   "mass" keyword returns).
; EXAMPLE:
;   For using with photoop system:
; 
;    ra=136.
;    dec=20.
;    obj= sdss_findobj(ra, dec, rerun=137, childobj=calibobj)
;    findspec, ra, dec, slist=slist
;    readspec, slist.plate, slist.fiberid, mjd=slist.mjd, zans=zans
;    kc= sdss_kcorrect(zans.z, calibobj=calibobj) 
;  
;   For reading tsobj structures:
;
;    tsobj=mrdfits('tsObj-01336-3-0456.fit',1,row=100)
;    findspec, tsobj.ra, tsobj.dec, slist=slist
;    readspec, slist.plate, slist.fiberid, mjd=slist.mjd, zans=zans
;    kc= sdss_kcorrect(zans.z, tsobj=tsobj) 
; 
; REVISION HISTORY:
;   07-Apr-2005  Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function sdss_kcorrect, redshift, nmgy=nmgy, ivar=ivar, mag=mag, err=err, $
                        calibobj=calibobj, tsobj=tsobj, flux=flux, $
                        band_shift=in_band_shift, chi2=chi2, coeffs=coeffs, $
                        rmaggies=rmaggies, omaggies=omaggies, $
                        oivar=oivar, vname=vname, mass=mass, mtol=mtol, $
                        absmag=absmag, amivar=amivar, omega0=omega0, $
                        omegal0=omegal0

common com_sdss_kcorrect, rmatrix, zvals, band_shift

if(n_params() lt 1 OR $
   (((keyword_set(nmgy) eq 0 OR keyword_set(ivar) eq 0)) AND $
    ((keyword_set(mag) eq 0 OR keyword_set(err) eq 0)) AND $
    (n_tags(calibobj) eq 0) AND $
    (n_tags(tsobj) eq 0))) $
  then begin
    doc_library, 'sdss_kcorrect'
    return, -1
endif 

;; need to reset rmatrix if band_shift changes
if(n_elements(in_band_shift) gt 0) then begin
    if(n_elements(band_shift) ne 0) then begin
        if(band_shift ne in_band_shift) then begin
            rmatrix=0
            zvals=0
        endif
        band_shift=in_band_shift
    endif else begin
        band_shift=in_band_shift
    endelse 
endif else begin
    if(n_elements(band_shift) eq 0) then $
      band_shift=0.
endelse

sdss_to_maggies, mgy, mgy_ivar, calibobj=calibobj, tsobj=tsobj, flux=flux

;; call kcorrect
kcorrect, mgy, mgy_ivar, redshift, kcorrect, band_shift=band_shift, $
  rmatrix=rmatrix, zvals=zvals, coeffs=coeffs, rmaggies=rmaggies, $
  vname=vname, mass=mass, mtol=mtol, absmag=absmag, amivar=amivar, $
  omega0=omega0, omegal0=omegal0

if(arg_present(omaggies)) then $
  omaggies=mgy
if(arg_present(oivar)) then $
  oivar=mgy_ivar

return, kcorrect

end
