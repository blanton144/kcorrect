;+
; NAME:
;   sdss2deep
; PURPOSE:
;   take SDSS data and return DEEP BRI at some redshift
; CALLING SEQUENCE:
;   bri= sdss2deep(sdss_redshift, deep_redshift, [, nmgy=, ivar=, mag=, err=, $
;                  calibobj=, tsobj=, flux=, chi2=, rmaggies=, $
;                  omaggies=, vname=, oivar=, mass=, mtol= ]
; INPUTS:
;   sdss_redshift - [N] redshifts of input 
;   deep_redshift - [N] redshifts of desired output
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
; OUTPUTS:
;   bri - [3, N] apparent magnitudes in BRI (AB)
;   mtol - [5, N] mass-to-light ratios from model in each band
;   mass - [N] total mass from model in each band
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
function sdss2deep, sdss_redshift, deep_redshift, nmgy=nmgy, ivar=ivar, $
                    mag=mag, err=err, calibobj=calibobj, tsobj=tsobj, $
                    flux=flux, chi2=chi2, coeffs=coeffs, rmaggies=rmaggies, $
                    omaggies=omaggies, oivar=oivar, vname=vname, $
                    mass=mass, mtol=mtol

common com_sdss2deep, rmatrix, zvals, band_shift

closest=1L

if(n_params() lt 1 OR $
   (((keyword_set(nmgy) eq 0 OR keyword_set(ivar) eq 0)) AND $
    ((keyword_set(mag) eq 0 OR keyword_set(err) eq 0)) AND $
    (n_tags(calibobj) eq 0) AND $
    (n_tags(tsobj) eq 0))) $
  then begin
    doc_library, 'sdss2deep'
    return, -1
endif 

sdss_filterlist=['sdss_u0.par', $
                 'sdss_g0.par', $
                 'sdss_r0.par', $
                 'sdss_i0.par', $
                 'sdss_z0.par']

kcdum=sdss_kcorrect(sdss_redshift,nmgy=nmgy, ivar=ivar, mag=mag, err=err, $
                    calibobj=calibobj, tsobj=tsobj, flux=flux, $
                    chi2=chi2, coeffs=coeffs, rmaggies=rmaggies, $
                    omaggies=omaggies, oivar=oivar, vname=vname, $
                    mass=mass, mtol=mtol, band_shift=band_shift)

; calculate the preliminaries
filterlist=['deep_B.par', $
            'deep_R.par', $
            'deep_I.par']
if(NOT keyword_set(rmatrix) OR NOT keyword_set(zvals)) then begin
    if(NOT keyword_set(vmatrix) OR NOT keyword_set(lambda)) then $
      k_load_vmatrix, vmatrix, lambda, vfile=vfile, lfile=lfile, $
      vpath=vpath, vname=vname
    k_projection_table,rmatrix,vmatrix,lambda,zvals,filterlist, $ 
      zmin=zmin,zmax=zmax,nz=nz,filterpath=filterpath
endif

; Reconstruct the magnitudes as observed by DEEP
k_reconstruct_maggies,coeffs, deep_redshift, $
  reconstruct_maggies,rmatrix=rmatrix,zvals=zvals

obands=lindgen(n_elements(filterlist))#replicate(1L, n_elements(sdss_redshift))

if(keyword_set(closest)) then begin
    lambda_in=k_lambda_eff(filterlist=sdss_filterlist)
    lambda_out=k_lambda_eff(filterlist=filterlist)
    for i=0L, n_elements(sdss_redshift)-1L do begin
        for j=0L, n_elements(lambda_out)-1L do begin
            dmin=min(abs(lambda_in/(1.+sdss_redshift[i])- $
                         lambda_out[j]/(1.+deep_redshift[i])), imin)
            obands[j, i]= imin
        endfor
    endfor
endif

offset=fltarr(n_elements(filterlist), n_elements(sdss_redshift))
for i=0L, n_elements(sdss_redshift)-1L do $
  for j=0L, n_elements(filterlist)-1L do $
  offset[j,i]=reconstruct_maggies[j,i]/rmaggies[obands[j,i],i]
offset=2.5*alog10(offset)

bri=fltarr(n_elements(filterlist), n_elements(sdss_redshift))
bri_ivar=fltarr(n_elements(filterlist), n_elements(sdss_redshift))
dm_sdss=lf_distmod(sdss_redshift, omega0=omega0, omegal0=omegal0)
dm_deep=lf_distmod(deep_redshift, omega0=omega0, omegal0=omegal0)
for j=0L, n_elements(filterlist)-1L do $
  bri[j,*]=-2.5*alog10(reconstruct_maggies[j,*])-dm_sdss+dm_deep 

for i=0L, n_elements(sdss_filterlist)-1L do begin
    for j=0L, n_elements(filterlist)-1L do begin
        ifrom=where(obands[j,*] eq i, nfrom)
        if(nfrom gt 0) then begin
            ig=where(oivar[i,ifrom] gt 0. AND $
                     omaggies[i,ifrom] gt 0., ng)
            if(ng gt 0) then begin
                ig=ifrom[ig]
                bri[j,ig]=-2.5*alog10(omaggies[i,ig])- $
                  dm_sdss[ig]+dm_deep[ig]-offset[j,ig]
                bri_ivar[j,ig]=omaggies[i,ig]^2*oivar[i,ig]*(0.4*alog(10.))^2
            endif
        endif
    endfor
endfor

return, bri

end
