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
;                  .REDDENING[5]
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
; OPTIONAL KEYWORDS:
;   lrg - do "luminous red galaxy" fit; this means changing vname to
;         'lrg1', using "model" fluxes, and ignoring the u-band
;         (setting ivar[0,*]=0); this uses a single template
;         appropriate for the SDSS Luminous Red Galaxy sample
; OUTPUTS:
;   kcorrect - [5, ngals] K-corrections in ugriz satisfying
;                m = M + DM(z) + K(z)
;              based on the best fit sum of templates
;   mtol - [5, ngals] current stellar mass-to-light ratios from model
;          in each band
;   mass - [ngals] total current stellar mass from model 
;   mets - [ngals] average metallicity in current stars 
;   intsfh - [ngals] total integrated star formation history
;   absmag - [5, ngals] absolute magnitude (for missing data, substitutes
;            model fit). (evolution correction *not* applied)
;   amivar - [5, ngals] inverse variance of absolute magnitude (for
;            missing data = 0)
; OPTIONAL OUTPUTS:
;   coeffs - [Nt, ngals] coefficients of fit
;   chi2 - chi^2 of fit
;   rmaggies - [5, ngals] reconstructed maggies from the fit (ugriz)
;   omaggies, oivar - [5, ngals] maggies and inverse variances used for fit
;                           (after extinction, AB correction, etc)  (ugriz)
;   b300 - [ngals] star-formation within last 300Myrs relative to average
;          star-formation rate
;   b1000 - [ngals] star-formation within last 1Gyrs relative to average
;           star-formation rate
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
;   argument "flux" to sdss_to_maggies to indicate which type of flux
;   to use. Note that if you get magnitudes like petroMag or modelMag
;   from the Catalog Archive Servers, these numbers are exactly like
;   the petroCounts and counts_model numbers in the tsObj structures.
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
                        oivar=oivar, vname=in_vname, mass=mass, mtol=mtol, $
                        absmag=absmag, amivar=amivar, omega0=omega0, $
                        omegal0=omegal0, lrg=lrg, mets=mets, b300=b300, $
                        b1000=b1000, intsfh=intsfh, silent=silent

common com_sdss_kcorrect, rmatrix, zvals, band_shift, vname, ermatrix

if(n_params() lt 1 OR $
   (((keyword_set(nmgy) eq 0 OR keyword_set(ivar) eq 0)) AND $
    ((keyword_set(mag) eq 0 OR keyword_set(err) eq 0)) AND $
    (n_tags(calibobj) eq 0) AND $
    (n_elements(coeffs) eq 0) AND $
    (n_tags(tsobj) eq 0))) $
  then begin
    doc_library, 'sdss_kcorrect'
    return, -1
endif 

if(keyword_set(lrg)) then $
  flux='model'

if(n_elements(in_vname) gt 0) then begin
    use_vname=in_vname
endif else begin
    if(keyword_set(lrg)) then $
      use_vname='lrg1' $
    else $
      use_vname='default'
endelse
if(n_elements(vname) gt 0) then begin
    if(vname ne use_vname) then begin
        rmatrix=0
        ermatrix=0
        zvals=0
    endif
endif
vname=use_vname

;; need to reset rmatrix if band_shift changes
if(n_elements(in_band_shift) gt 0) then $
  use_band_shift=in_band_shift $
else $
  use_band_shift=0. 
if(n_elements(band_shift) gt 0) then begin
    if(band_shift ne use_band_shift) then begin
        rmatrix=0
        ermatrix=0
        zvals=0
    endif
endif
band_shift=use_band_shift

if(keyword_set(mag) AND keyword_set(err)) then begin
    mgy=(10.D)^(-(0.4D)*(mag))
    mags_ivar=1./err^2
    mgy_ivar= mags_ivar/(0.4*alog(10.)*mgy)^2.
endif
if(keyword_set(nmgy) AND keyword_set(ivar)) then begin
    mgy=1.e-9*nmgy
    mgy_ivar=1.e+18*ivar
endif
if(n_tags(tsobj) gt 0 OR n_tags(calibobj) gt 0) then $
  sdss_to_maggies, mgy, mgy_ivar, calibobj=calibobj, tsobj=tsobj, flux=flux

if(keyword_set(lrg)) then $
  if(keyword_set(mgy_ivar)) then $
  mgy_ivar[0,*]=0.

;; call kcorrect
kcorrect, mgy, mgy_ivar, redshift, kcorrect, band_shift=band_shift, $
  rmatrix=rmatrix, zvals=zvals, coeffs=coeffs, rmaggies=rmaggies, $
  vname=vname, mass=mass, mtol=mtol, absmag=absmag, amivar=amivar, $
  omega0=omega0, omegal0=omegal0, chi2=chi2, mets=mets, b300=b300, $
  intsfh=intsfh, b1000=b1000, silent=silent

if(arg_present(omaggies)) then $
  omaggies=mgy
if(arg_present(oivar)) then $
  oivar=mgy_ivar

return, kcorrect

end
