;+
; NAME:
;   deep_kcorrect
; PURPOSE:
;   calculate K-corrections for standard DEEP input (from BRI to NUB)
; CALLING SEQUENCE:
;   kcorrect= deep_kcorrect(redshift [, nmgy=, ivar=, mag=, err=, $
;                           zcat=, /sdss, band_shift=, chi2=, rmaggies=, $
;                           omaggies=, vname=, oivar=, mass=, mtol=, $
;                           absmag=, amivar= ])
; INPUTS:
;   redshift - [N] redshifts
;   zcat - [N] DEEP zcat-style structure, containing:
;                  .MAGB
;                  .MAGR
;                  .MAGI
;                  .SFD_EBV
;   nmgy, ivar - [3, N] nanomaggies, Galactic-reddening corrected, and inverse
;                variance of same
;   mag, err - [3, N] Pogson magnitudes, Galactic-reddening corrected and
;              errors of same
; OPTIONAL INPUTS:
;   band_shift    - blueshift of bandpasses to apply (to get ^{z}b
;                   type bands) [default 0.]
;   vname - name of fit to use (defaults to 'default')
; OUTPUTS:
;   kcorrect - [3, N] K-corrections from BRI to NUV, U, and B (or NUV,
;              u, g if /sdss is set) satisfying
;                m_R = M_Q + DM(z) + K_QR(z)
;              based on the best fit sum of templates. 
;              
;   mtol - [3, N] mass-to-light ratios from model in each output band
;   mass - [N] total mass from model in each band
;   absmag - [3, N] absolute magnitude (for missing data, substitutes
;            model fit) in each output band
;   amivar - [3, N] inverse variance of absolute magnitude (for
;            missing data = 0) in each output band
; OPTIONAL OUTPUTS:
;   coeffs - coefficients of fit
;   chi2 - chi^2 of fit
;   rmaggies - [3, N] reconstructed maggies from the fit (BRI)
;   omaggies, oivar - [3, N] maggies and inverse variances used for fit
;                           (after extinction, AB correction, etc)
;                           (BRI)
; COMMENTS:
;   This is a simple wrapper on kcorrect.pro which is almost always
;   just what you want. It keeps a version of rmatrix and zvals in
;   memory to save time, recalculating them each time you change
;   band_shift.
;
;   You must specify nmgy,ivar OR mag,err OR zcat. If
;   nmgy or mag, make sure they are AB calibrated and Galactic
;   extinction corrected.
;
;   Uses deep_to_maggies to convert zcat structure to Galactic
;   extinction corrected maggies with errors.
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
function sdss_kcorrect, redshift, nmgy=nmgy, ivar=ivar, mag=mag, err=err, $
                        calibobj=calibobj, tsobj=tsobj, flux=flux, $
                        band_shift=in_band_shift, chi2=chi2, coeffs=coeffs, $
                        rmaggies=rmaggies, omaggies=omaggies, $
                        oivar=oivar, vname=vname, mass=mass, mtol=mtol, $
                        absmag=absmag, amivar=amivar

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

sdss_to_maggies, mgy, mgy_ivar, calibobj=calibobj, tsobj=tsobj, flux=flux

;; call kcorrect
kcorrect, mgy, mgy_ivar, redshift, kcorrect, band_shift=band_shift, $
  rmatrix=rmatrix, zvals=zvals, coeffs=coeffs, rmaggies=rmaggies, $
  vname=vname, mass=mass, mtol=mtol, absmag=absmag, amivar=amivar

if(arg_present(omaggies)) then $
  omaggies=mgy
if(arg_present(oivar)) then $
  oivar=mgy_ivar

return, kcorrect

end
