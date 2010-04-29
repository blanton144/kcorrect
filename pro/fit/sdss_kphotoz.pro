;+
; NAME:
;   sdss_kphotoz
; PURPOSE:
;   calculate photometric redshifts from SDSS input (v4 use ONLY for LRGs)
; CALLING SEQUENCE:
;   photoz= sdss_kphotoz([nmgy=, ivar=, mag=, err=, $
;                         calibobj=, tsobj=, flux=, band_shift=,$
;                         chi2=, rmaggies=, omaggies=, vname=, $
;                         oivar=, mass=, mtol=, absmag=, amivar=, $
;                         omega0=, omegal0=, kcorrect= ])
; INPUTS:
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
;   cas - [N] CAS-style structure, with:
;                  .MODELMAG_[X]  (for X = U,G,R,I, and Z)
;                  .MODELMAGERR[X]  (for X = U,G,R,I, and Z)
;                  .PSFMAG_[X]  (for X = U,G,R,I, and Z)
;                  .PSFMAGERR[X]  (for X = U,G,R,I, and Z)
;                  .PETROMAG_[X]  (for X = U,G,R,I, and Z)
;                  .PETROMAGERR[X]  (for X = U,G,R,I, and Z)
;                  .EXTINCTION_[X]  (for X = U,G,R,I, and Z)
;   nmgy, ivar - [5, N] nanomaggies, Galactic-reddening corrected, and inverse
;                variance of same
;   mag, err - [5, N] asinh magnitudes, Galactic-reddening corrected and
;              errors of same
; OPTIONAL INPUTS:
;   flux - use this version of the fluxes ('PETRO', 'MODEL', or 'PSF')
;          [defaults to 'MODEL'] if tsobj or calibobj keywords are
;          used 
;   band_shift    - blueshift of bandpasses to apply for K-corrections
;                   returned (to get ^{z}b type bands) [default 0.] 
;   vname - name of fit to use (defaults to 'default')
;   omega0, omegal0 - cosmological parameters for calculating distance
;                     moduli [default 0.3, 0.7]
; OPTIONAL KEYWORDS:
;   lrg - do "luminous red galaxy" fit; this means changing vname to
;         'lrg1' and ignoring the u-band (setting ivar[0,*]=0);
;         this is the right way to get photometric redshifts for red
;         galaxies; if you take things with good chi2 values from
;         this fit, they will have good photo-z's
; OUTPUTS:
;   photoz - [N] photometric redshifts
;   kcorrect - [5, ngals] K-corrections in ugriz satisfying
;                m = M + DM(z) + K(z)
;              based on the best fit sum of templates
;   mtol - [5, ngals] mass-to-light ratios from model in each band
;   mass - [ngals] total current stellar mass from model in each band
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
;   For the v4 templates, this has not be sufficiently well-tested.
;   Use it ONLY with the /lrg flag described above (for which the
;   template is known to be good).
; 
;   This is a simple wrapper on kphotoz.pro which is almost always
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
; 
;     1 solar mass / (D/10pc)^2 
;
;   That is, sum the coefficients and multiply by (D/10pc)^2 to get
;   TOTAL INTEGRATED STAR FORMATION. (In fact, for Omega0=0.3 and
;   OmegaL0=0.7, this is what the "mass" keyword returns). Note that
;   the total integrated star formation DIFFERS from the current
;   stellar mass --- which is returned in the mass and mtol variables.
; REVISION HISTORY:
;   07-June-2005  Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function sdss_kphotoz, nmgy=nmgy, ivar=ivar, mag=mag, err=err, $
                       calibobj=calibobj, tsobj=tsobj, cas=cas, flux=flux, $
                       band_shift=in_band_shift, chi2=chi2, coeffs=coeffs, $
                       rmaggies=rmaggies, omaggies=omaggies, $
                       oivar=oivar, vname=in_vname, mass=mass, mtol=mtol, $
                       absmag=absmag, amivar=amivar, omega0=omega0, $
                       omegal0=omegal0, kcorrect=kcorrect, lrg=lrg, $
                       noprior=noprior

common com_sdss_kphotoz, rmatrix, zvals, band_shift, vname

if((((keyword_set(nmgy) eq 0 OR keyword_set(ivar) eq 0)) AND $
    ((keyword_set(mag) eq 0 OR keyword_set(err) eq 0)) AND $
    (n_tags(calibobj) eq 0) AND $
    (n_tags(cas) eq 0) AND $
    (n_tags(tsobj) eq 0))) $
  then begin
    doc_library, 'sdss_kphotoz'
    return, -1
endif 

if(NOT keyword_set(flux)) then flux='MODEL'

if(keyword_set(lrg)) then in_vname='lrg1'

if(n_elements(in_vname) gt 0) then begin
    if(n_elements(vname) gt 0) then begin
        if(vname ne in_vname) then begin
            rmatrix=0
            zvals=0
        endif
    endif
    vname=in_vname
endif else begin
    if(n_elements(vname) eq 0) then $
      vname='default'
endelse

if(keyword_set(mag) AND keyword_set(err)) then begin
    mgy=(10.D)^(-(0.4D)*(mag))
    mags_ivar=1./err^2
    mgy_ivar= mags_ivar/(0.4*alog(10.)*mgy)^2.
endif
if(keyword_set(nmgy) AND keyword_set(ivar)) then begin
    mgy=1.e-9*nmgy
    mgy_ivar=1.e+18*ivar
endif
if(n_tags(tsobj) gt 0 OR n_tags(calibobj) gt 0 OR n_tags(cas) gt 0) then $
  sdss_to_maggies, mgy, mgy_ivar, calibobj=calibobj, tsobj=tsobj, flux=flux, cas=cas

if(keyword_set(lrg)) then mgy_ivar[0,*]=0.

;; call kcorrect
kphotoz, mgy, mgy_ivar, photoz, vname=vname, $
  rmatrix=rmatrix, zvals=zvals, coeffs=coeffs, $
  vmatrix=vmatrix, lambda=lambda, noprior=noprior

kcorrect= sdss_kcorrect(photoz, nmgy=1.e+9*mgy, ivar=mgy_ivar*1.e-18, $
                        band_shift=in_band_shift, chi2=chi2, coeffs=coeffs, $
                        rmaggies=rmaggies, vname=vname, mass=mass, mtol=mtol, $
                        absmag=absmag, amivar=amivar, omega0=omega0, $
                        omegal0=omegal0, lrg=lrg)

if(arg_present(omaggies)) then $
  omaggies=mgy
if(arg_present(oivar)) then $
  oivar=mgy_ivar

return, photoz

end
