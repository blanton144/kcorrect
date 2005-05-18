;+
; NAME:
;   galex_kcorrect
; PURPOSE:
;   calculate K-corrections for GALEX+SDSS input
; CALLING SEQUENCE:
;   kcorrect= galex_kcorrect(redshift [, nmgy=, ivar=, mag=, err=, $
;                            galex=, calibobj=, tsobj=, flux=, band_shift=,
;                            chi2=, rmaggies=, omaggies=, oivar=, $
;                            mass=, mtol=, absmag=, amivar=])
; INPUTS:
;   redshift - [N] redshifts
;   galex - [N] GALEX "mcat" style input, with:
;            .ALPHA_J2000
;            .DEC_J2000
;            .NUV_MAG
;            .NUV_MAGERR
;            .FUV_MAG
;            .FUV_MAGERR
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
;   nmgy, ivar - [7, N] nanomaggies, Galactic-reddening corrected, and inverse
;                variance of same
;   mag, err - [7, N] standard Pogson magnitudes, Galactic-reddening
;              corrected, and errors of same
; OPTIONAL INPUTS:
;   flux - use this version of the SDSS fluxes ('PETRO', 'MODEL', or 'PSF')
;          [defaults to 'PETRO'] if tsobj or calibobj keywords are
;          used 
;   band_shift - blueshift of bandpasses to apply (to get ^{z}b
;                type bands) [default 0.]
;   vname - name of fit to use (defaults to 'default')
; OUTPUTS:
;   kcorrect - [7, N] K-corrections in FNugriz satisfying
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
;   rmaggies - [7, N] reconstructed maggies from the fit (FNugriz)
;   omaggies, oivar - [7, N] maggies and inverse variances used for fit
;                     (after extinction, AB correction, etc) (FNugriz)
; COMMENTS:
;   This is a simple wrapper on kcorrect.pro which is designed for
;   users with GALEX data matched to SDSS data. It will deal
;   appropriately if you leave out the SDSS or GALEX data (but not
;   both!).
;
;   Uses galex_to_maggies to convert a GALEX-style catalog to maggies
;   and inverse variances.
;
;   Uses sdss_to_maggies to convert tsobj or calibobj structure to
;   AB, Galactic extinction corrected maggies. Passes optional
;   argument "flux" to sdss_to_maggies.
;
;   You must specify nmgy,ivar OR mag,err OR calibobj OR tsobj OR
;   galex.  If nmgy or mag, make sure they are Galactic extinction
;   corrected and AB calibrated.
;
;   For v4_0b templates and later, coefficients are in units of:
;     1 solar mass / (D/10pc)^2
;   That is, sum the coefficients and multiply by (D/10pc)^2 to get
;   masses. (In fact, for Omega0=0.3 and OmegaL0=0.7, this is what the
;   "mass" keyword returns).
; REVISION HISTORY:
;   07-Apr-2005  Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function galex_kcorrect, redshift, nmgy=nmgy, ivar=ivar, mag=mag, err=err, $
                         calibobj=calibobj, tsobj=tsobj, flux=flux, $
                         band_shift=in_band_shift, chi2=chi2, coeffs=coeffs, $
                         rmaggies=rmaggies, omaggies=omaggies, $
                         oivar=oivar, galex=galex, vname=vname, $
                         mass=mass, mtol=mtol, absmag=absmag, amivar=amivar

common com_galex_kcorrect, rmatrix, zvals, band_shift

if(n_params() lt 1 OR $
   (((keyword_set(nmgy) eq 0 OR keyword_set(ivar) eq 0)) AND $
    ((keyword_set(mag) eq 0 OR keyword_set(err) eq 0)) AND $
    (n_tags(calibobj) eq 0) AND $
    (n_tags(galex) eq 0) AND $
    (n_tags(tsobj) eq 0))) $
  then begin
    doc_library,'galex_kcorrect'
    return, -1
endif 

;; interpret band_shift
if(NOT keyword_set(in_band_shift)) then in_band_shift=0.

;; need to reset rmatrix if band_shift changes
if(n_elements(band_shift) ne 0) then begin
    if(band_shift ne in_band_shift) then begin
       rmatrix=0
       zvals=0
    endif
endif else begin
    band_shift=in_band_shift
endelse 

;; set up maggies, and initialize to mag or nmgy input if available
mgy=fltarr(7, n_elements(redshift))
mgy_ivar=fltarr(7, n_elements(redshift))
if(n_elements(mag) gt 0) then begin
    mgy[*,*]=10.^(-0.4*mag)
    mgy_ivar[*,*]=1./(0.4*alog(10.)*mgy*err)^2.
endif
if(n_elements(nmgy) gt 0) then begin
    mgy[*,*]=nmgy*1.e-9
    mgy_ivar[*,*]=ivar*1.e+18
endif

;; get GALEX stuff from GALEX structure
if(n_tags(galex) gt 0) then begin
    galex_to_maggies, galex, galex_mgy, galex_mgy_ivar
    mgy[0:1,*]=galex_mgy
    mgy_ivar[0:1,*]=galex_mgy_ivar
endif

;; get SDSS stuff
if(n_tags(tsobj) gt 0 OR n_tags(calibobj) gt 0) then begin
    sdss_to_maggies, sdss_mgy, sdss_mgy_ivar, tsobj=tsobj, calibobj=calibobj, $
      flux=flux
    mgy[2:6,*]=sdss_mgy
    mgy_ivar[2:6,*]=sdss_mgy_ivar
endif

;; call kcorrect
filterlist=['galex_FUV.par', 'galex_NUV.par', 'sdss_u0.par', $
            'sdss_g0.par', 'sdss_r0.par', 'sdss_i0.par', 'sdss_z0.par']
kcorrect, mgy, mgy_ivar, redshift, kcorrect, band_shift=band_shift, $
  rmatrix=rmatrix, zvals=zvals, coeffs=coeffs, rmaggies=rmaggies, $
  filterlist=filterlist, vname=vname

if(arg_present(omaggies)) then $
  omaggies=mgy
if(arg_present(oivar)) then $
  oivar=mgy_ivar

return, kcorrect

end
