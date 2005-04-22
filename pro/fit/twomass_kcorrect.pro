;+
; NAME:
;   twomass_kcorrect
; PURPOSE:
;   calculate K-corrections for 2MASS+SDSS input
; CALLING SEQUENCE:
;   kcorrect= twomass_kcorrect(redshift [, nmgy=, ivar=, mag=, err=, $
;                              twomass=, calibobj=, tsobj=, flux=, band_shift=,
;                              chi2=, rmaggies=, omaggies=, oivar=])
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
; OUTPUTS:
;   kcorrect - [8, N] K-corrections in ugrizJHK satisfying
;                m = M + DM(z) + K(z)
;              based on the best fit sum of templates
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
; REVISION HISTORY:
;   07-Apr-2005  Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function twomass_kcorrect, redshift, nmgy=nmgy, ivar=ivar, mag=mag, err=err, $
                           calibobj=calibobj, tsobj=tsobj, flux=flux, $
                           band_shift=in_band_shift, chi2=chi2, $
                           coeffs=coeffs, rmaggies=rmaggies, $
                           omaggies=omaggies, oivar=oivar, twomass=twomass, $
                           vname=vname

common com_twomass_kcorrect, rmatrix, zvals, band_shift

if(n_params() lt 1 OR $
   (((keyword_set(nmgy) eq 0 OR keyword_set(ivar) eq 0)) AND $
    ((keyword_set(mag) eq 0 OR keyword_set(err) eq 0)) AND $
    (n_tags(calibobj) eq 0) AND $
    (n_tags(twomass) eq 0) AND $
    (n_tags(tsobj) eq 0))) $
  then begin
    doc_library,'twomass_kcorrect'
    return, -1
endif 

;; need to reset rmatrix if band_shift changes
if(n_elements(band_shift) ne 0) then begin
    if(band_shift ne in_band_shift) then begin
       rmatrix=0
       zvals=0
    endif
endif else begin
    band_shift=in_band_shift
endelse 

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
  filterlist=filterlist, vname=vname

if(arg_present(omaggies)) then $
  omaggies=mgy
if(arg_present(oivar)) then $
  oivar=mgy_ivar

return, kcorrect

end
