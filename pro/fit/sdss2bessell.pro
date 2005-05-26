;+
; NAME:
;   sdss2bessell
; PURPOSE:
;   take SDSS data and return rest-frame UBVRI data
; CALLING SEQUENCE:
;   kcorrect= sdss2bessell(redshift [, nmgy=, ivar=, mag=, err=, $
;                          calibobj=, tsobj=, flux=, chi2=, rmaggies=, $
;                          omaggies=, vname=, oivar=, mass=, mtol=, $
;                          absmag=, amivar=, band_shift=, /vega ])
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
;   vname - name of fit to use (defaults to 'default')
;   band_shift - evaluate UBVRI shifted blue by 1+band_shift
; OPTIONAL KEYWORDS:
;   /vega - output Vega magnitudes, fluxes (STILL TAKES AB INPUTS!)
; OUTPUTS:
;   kcorrect - [5, N] K-corrections from ugriz to UBVRI; e.g.:
;                   M_U = m_u - DM(z) - K_{uU}(z)
;   mtol - [5, N] mass-to-light ratios from model in each band
;   mass - [N] total mass from model in each band
;   absmag - [5, N] absolute magnitude (for missing data, substitutes
;            model fit) in UBVRI
;   amivar - [5, N] inverse variance of absolute magnitude (for
;            missing data = 0) in UBVRI
; OPTIONAL OUTPUTS:
;   coeffs - coefficients of fit
;   chi2 - chi^2 of fit
;   rmaggies - [5, N] reconstructed maggies from the fit (ugriz)
;   omaggies, oivar - [5, N] maggies and inverse variances used for fit
;                           (after extinction, AB correction, etc)  (ugriz)
; COMMENTS:
;   This is a simple wrapper on kcorrect.pro. It keeps a version of
;   rmatrix and zvals in memory to save time.
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
; REVISION HISTORY:
;   07-Apr-2005  Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function sdss2bessell, redshift, nmgy=nmgy, ivar=ivar, mag=mag, err=err, $
                       calibobj=calibobj, tsobj=tsobj, flux=flux, $
                       chi2=chi2, coeffs=coeffs, rmaggies=rmaggies, $
                       omaggies=omaggies, oivar=oivar, vname=vname, $
                       mass=mass, mtol=mtol, absmag=absmag, amivar=amivar, $
                       band_shift=in_band_shift, vega=vega

common com_sdss2bessell, rmatrix, zvals, band_shift

if(n_params() lt 1 OR $
   (((keyword_set(nmgy) eq 0 OR keyword_set(ivar) eq 0)) AND $
    ((keyword_set(mag) eq 0 OR keyword_set(err) eq 0)) AND $
    (n_tags(calibobj) eq 0) AND $
    (n_tags(tsobj) eq 0))) $
  then begin
    doc_library, 'sdss2bessell'
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

kcdum=sdss_kcorrect(redshift,nmgy=nmgy, ivar=ivar, mag=mag, err=err, $
                    calibobj=calibobj, tsobj=tsobj, flux=flux, $
                    chi2=chi2, coeffs=coeffs, rmaggies=rmaggies, $
                    omaggies=omaggies, oivar=oivar, vname=vname, $
                    mass=mass, mtol=mtol, band_shift=band_shift)

; calculate the preliminaries
filterlist=['bessell_U.par', $
            'bessell_B.par', $
            'bessell_V.par', $
            'bessell_R.par', $
            'bessell_I.par']
if(NOT keyword_set(rmatrix) OR NOT keyword_set(zvals)) then begin
    if(NOT keyword_set(vmatrix) OR NOT keyword_set(lambda)) then $
      k_load_vmatrix, vmatrix, lambda, vfile=vfile, lfile=lfile, $
      vpath=vpath, vname=vname
    k_projection_table,rmatrix,vmatrix,lambda,zvals,filterlist, $ 
      zmin=zmin,zmax=zmax,nz=nz,filterpath=filterpath
endif

; Reconstruct the magnitudes as observed and in the rest frame
k_reconstruct_maggies,coeffs,replicate(band_shift,n_elements(redshift)), $
  reconstruct_maggies,rmatrix=rmatrix,zvals=zvals
reconstruct_maggies=reconstruct_maggies/(1.+band_shift)

kcorrect=reconstruct_maggies/rmaggies
kcorrect=2.5*alog10(kcorrect)

if(keyword_set(vega)) then begin
    for i=0L, n_elements(filterlist)-1L do begin
        v2ab=k_vega2ab(filterlist=filterlist[i], /kurucz, $
                       band_shift=band_shift)
        ;; you need to convert FROM AB TO VEGA. the K-correction is
        ;; M=m-DM-K, so adding v2ab to K makes M=m-DM-K-v2ab,
        ;; equivalent to M=m-DM-K+ab2v
        kcorrect[i,*]=kcorrect[i,*]+v2ab[0]
    endfor
endif

smaggies=10.^(-0.4*k_solar_magnitudes(filterlist=filterlist))
mtol=fltarr(n_elements(filterlist), n_elements(redshift))
mm=total(coeffs,1)
for i=0L, n_elements(filterlist)-1L do $
  mtol[i,*]=mm/reconstruct_maggies[i,*]*smaggies[i]

if(arg_present(absmag)) then begin
    absmag=fltarr(n_elements(filterlist), n_elements(redshift))
    amivar=fltarr(n_elements(filterlist), n_elements(redshift))
    for i=0L, n_elements(filterlist)-1L do $
      absmag[i,*]=-2.5*alog10(reconstruct_maggies[i,*])- $
      lf_distmod(redshift, omega0=omega0, omegal0=omegal0)
    if(keyword_set(vega)) then begin
        for i=0L, n_elements(filterlist)-1L do begin
            v2ab=k_vega2ab(filterlist=filterlist[i], /kurucz, $
                           band_shift=band_shift)
            absmag[i,*]=absmag[i,*]-v2ab[0]
        endfor
    endif
    for i=0L, n_elements(filterlist)-1L do begin
        ig=where(oivar[i,*] gt 0. AND omaggies[i,*] gt 0., ng)
        if(ng gt 0) then begin
            absmag[i,ig]=-2.5*alog10(omaggies[i,ig])- $
              lf_distmod(redshift[ig], omega0=omega0, omegal0=omegal0)- $
              kcorrect[i,ig]
            amivar[i,ig]=omaggies[i,ig]^2*oivar[i,ig]* $
              (0.4*alog(10.))^2
        endif 
    endfor
endif

return, kcorrect

end
