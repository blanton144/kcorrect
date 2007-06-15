;+
; NAME:
;   twomass_kcorrect
; PURPOSE:
;   calculate K-corrections for 2MASS+SDSS input
; CALLING SEQUENCE:
;   kcorrect= twomass_kcorrect(redshift [, nmgy=, ivar=, mag=, err=, $
;                              twomass=, calibobj=, tsobj=, flux=, band_shift=,
;                              chi2=, rmaggies=, omaggies=, oivar=, $
;                              mass=, mtol=, absmag=, amivar=, $
;                              omega0=, omegal0=])
; INPUTS:
;   redshift - [N] redshifts
;   twomass - [N] 2MASS "mcat" style input, with:
;               .RA
;               .DECL
;               .J_M_EXT
;               .J_MSIG_EXT
;               .J_FLG_EXT
;               .H_M_EXT
;               .H_MSIG_EXT
;               .H_FLG_EXT
;               .K_M_EXT
;               .K_MSIG_EXT
;               .K_FLG_EXT
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
;   nmgy, ivar - [8, N] nanomaggies, Galactic-reddening corrected, and inverse
;                variance of same
;   mag, err - [8, N] standard Pogson magnitudes, Galactic-reddening
;              corrected, and errors of same
; OPTIONAL INPUTS:
;   flux - use this version of the SDSS fluxes ('PETRO', 'MODEL', or 'PSF')
;          [defaults to 'PETRO'] if tsobj or calibobj keywords are
;          used 
;   band_shift - blueshift of bandpasses to apply (to get ^{z}b
;                type bands) [default 0.]
;   vname - name of fit to use (defaults to 'default')
;   omega0, omegal0 - cosmological parameters for calculating distance
;                     moduli [default 0.3, 0.7]
; OUTPUTS:
;   kcorrect - [8, N] K-corrections in ugrizJHK satisfying
;                m = M + DM(z) + K(z)
;              based on the best fit sum of templates
;   mtol - [8, ngals] mass-to-light ratios from model in each band
;   mass - [ngals] total mass from model in each band
;   absmag - [8, ngals] absolute magnitude (for missing data, substitutes
;            model fit)
;   amivar - [8, ngals] inverse variance of absolute magnitude (for
;            missing data = 0)
; OPTIONAL OUTPUTS:
;   coeffs - coefficients of fit
;   chi2 - chi^2 of fit
;   rmaggies - [8, N] reconstructed maggies from the fit (ugrizJHK)
;   omaggies, oivar - [8, N] maggies and inverse variances used for fit
;                     (after extinction, AB correction, etc)  (ugrizJHK)
; COMMENTS:
;   This is a simple wrapper on kcorrect.pro which is designed for
;   users with 2MASS data matched to SDSS data. It will deal
;   appropriately if you leave out the SDSS or 2MASS data (but not
;   both!).
;
;   Uses twomass_to_maggies to convert a 2MASS-style catalog to maggies
;   and inverse variances.
;
;   Uses sdss_to_maggies to convert tsobj or calibobj structure to
;   AB, Galactic extinction corrected maggies. Passes optional
;   argument "flux" to sdss_to_maggies.
;
;   You must specify nmgy,ivar OR mag,err OR calibobj OR tsobj OR
;   twomass.  If nmgy or mag, make sure they are Galactic extinction
;   corrected and AB calibrated.
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
;   07-Apr-2005  Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function twomass_kcorrect, redshift, nmgy=nmgy, ivar=ivar, mag=mag, err=err, $
                           calibobj=calibobj, tsobj=tsobj, flux=flux, $
                           band_shift=in_band_shift, chi2=chi2, $
                           coeffs=coeffs, rmaggies=rmaggies, $
                           omaggies=omaggies, twomass=twomass, $
                           oivar=oivar, galex=galex, vname=in_vname, $
                           mass=mass, mtol=mtol, absmag=absmag, $
                           amivar=amivar, omega0=omega0, omegal0=omegal0, $
                           b300=b300, b1000=b1000, intsfh=intsfh, mets=mets

common com_twomass_kcorrect, rmatrix, zvals, band_shift, vname

if(n_params() lt 1 OR $
   (((keyword_set(nmgy) eq 0 OR keyword_set(ivar) eq 0)) AND $
    ((keyword_set(mag) eq 0 OR keyword_set(err) eq 0)) AND $
    (n_tags(calibobj) eq 0) AND $
    (n_tags(galex) eq 0) AND $
    (n_tags(tsobj) eq 0))) $
  then begin
    doc_library,'twomass_kcorrect'

    return, -1
endif 

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
        zvals=0
    endif
endif
band_shift=use_band_shift

;; set up maggies arrays
mgy=fltarr(8, n_elements(redshift))
mgy_ivar=fltarr(8, n_elements(redshift))
if(n_elements(mag) gt 0) then begin
    mgy[*,*]=10.^(-0.4*mag)
    
    mgy_ivar[*,*]=1./(0.4*alog(10.)*mgy*err)^2.
endif
if(n_elements(nmgy) gt 0) then begin
    mgy[*,*]=nmgy*1.e-9
    mgy_ivar[*,*]=nmgy_ivar*1.e+18
endif

;; get 2MASS stuff from twomass structure
if(n_tags(twomass) gt 0) then begin
    twomass_to_maggies, twomass, twomass_mgy, twomass_mgy_ivar
    mgy[5:7,*]=twomass_mgy
    mgy_ivar[5:7,*]=twomass_mgy_ivar
endif

;; get SDSS stuff
if(n_tags(tsobj) gt 0 OR n_tags(calibobj) gt 0) then begin
    sdss_to_maggies, sdss_mgy, sdss_mgy_ivar, tsobj=tsobj, calibobj=calibobj, $
      flux=flux
    mgy[0:4,*]=sdss_mgy
    mgy_ivar[0:4,*]=sdss_mgy_ivar
endif

;; call kcorrect
filterlist=['sdss_u0.par','sdss_g0.par', 'sdss_r0.par', 'sdss_i0.par', $
            'sdss_z0.par', 'twomass_J.par', 'twomass_H.par', 'twomass_Ks.par']
kcorrect, mgy, mgy_ivar, redshift, kcorrect, band_shift=band_shift, $
  rmatrix=rmatrix, zvals=zvals, coeffs=coeffs, rmaggies=rmaggies, $
  filterlist=filterlist, vname=vname, mass=mass, mtol=mtol, $
  absmag=absmag, amivar=amivar, omega0=omega0, omegal0=omegal0, $
  mets=mets, b300=b300, intsfh=intsfh, b1000=b1000

if(arg_present(omaggies)) then $
  omaggies=mgy
if(arg_present(oivar)) then $
  oivar=mgy_ivar

return, kcorrect

end
