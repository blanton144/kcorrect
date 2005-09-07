;+
; NAME:
;   sdss2goods
; PURPOSE:
;   take SDSS data and return GOODS BVizJHK at some redshift
; CALLING SEQUENCE:
;   mag= sdss2goods(sdss_redshift, goods_redshift, [, nmgy=, ivar=, $
;                   mag=, err=, calibobj=, tsobj=, flux=, chi2=, $
;                   rmaggies=, omaggies=, vname=, oivar=, mass=, mtol= ]
; INPUTS:
;   sdss_redshift - [N] redshifts of input 
;   goods_redshift - [N] redshifts of desired output
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
;   mag - [7, N] apparent magnitudes in BVizJHK (AB)
;   mass - [N] total current stellar mass from model in each band
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
function sdss2goods, sdss_redshift, goods_redshift, nmgy=nmgy, ivar=ivar, $
                    mag=mag, err=err, calibobj=calibobj, tsobj=tsobj, $
                    flux=flux, chi2=chi2, coeffs=coeffs, rmaggies=rmaggies, $
                    omaggies=omaggies, oivar=oivar, vname=vname, $
                    mass=mass, mtol=mtol

common com_sdss2goods, rmatrix, zvals, band_shift

if(n_params() lt 1 OR $
   (((keyword_set(nmgy) eq 0 OR keyword_set(ivar) eq 0)) AND $
    ((keyword_set(mag) eq 0 OR keyword_set(err) eq 0)) AND $
    (n_tags(calibobj) eq 0) AND $
    (n_tags(tsobj) eq 0))) $
  then begin
    doc_library, 'sdss2goods'
    return, -1
endif 

kcdum=sdss_kcorrect(sdss_redshift,nmgy=nmgy, ivar=ivar, mag=mag, err=err, $
                    calibobj=calibobj, tsobj=tsobj, flux=flux, $
                    chi2=chi2, coeffs=coeffs, rmaggies=rmaggies, $
                    omaggies=omaggies, oivar=oivar, vname=vname, $
                    mass=mass, mtol=mtol, band_shift=band_shift)

; calculate the preliminaries
filterlist=['goods_acs_f435w.par', $
            'goods_acs_f606w.par', $
            'goods_acs_f775w.par', $
            'goods_acs_f850lp.par', $
            'goods_J_isaac_etc.par', $
            'goods_H_isaac_etc.par', $
            'goods_Ks_isaac_etc.par']
if(NOT keyword_set(rmatrix) OR NOT keyword_set(zvals)) then begin
    if(NOT keyword_set(vmatrix) OR NOT keyword_set(lambda)) then $
      k_load_vmatrix, vmatrix, lambda, vfile=vfile, lfile=lfile, $
      vpath=vpath, vname=vname
    k_projection_table,rmatrix,vmatrix,lambda,zvals,filterlist, $ 
      zmin=zmin,zmax=zmax,nz=nz,filterpath=filterpath
endif

; Reconstruct the magnitudes as observed by GOODS
k_reconstruct_maggies,coeffs, goods_redshift, $
  reconstruct_maggies,rmatrix=rmatrix,zvals=zvals

sdss_filterlist=['sdss_u0.par', $
                 'sdss_g0.par', $
                 'sdss_r0.par', $
                 'sdss_i0.par', $
                 'sdss_z0.par']
lambda_in=k_lambda_eff(filterlist=sdss_filterlist)
lambda_out=k_lambda_eff(filterlist=filterlist, band_shift=band_shift)
obands=lonarr(n_elements(filterlist), n_elements(goods_redshift))
offset=fltarr(7, n_elements(goods_redshift))
for i=0L, n_elements(goods_redshift)-1L do begin
    for j=0L, n_elements(lambda_out)-1L do begin
        dmin=min(abs(lambda_in/(1.+sdss_redshift[i])- $
                     lambda_out[j]/(1.+goods_redshift[i])), imin)
        obands[j, i]= imin
        offset[j,i]=reconstruct_maggies[j,i]/rmaggies[obands[j,i],i]
    endfor
endfor
offset=2.5*alog10(offset)

mag=fltarr(n_elements(filterlist), n_elements(sdss_redshift))
mag_ivar=fltarr(n_elements(filterlist), n_elements(sdss_redshift))
for i=0L, n_elements(filterlist)-1L do $
  mag[i,*]=-2.5*alog10(reconstruct_maggies[i,*])- $
  (lf_distmod(sdss_redshift, omega0=omega0, omegal0=omegal0))[0]+ $
  (lf_distmod(goods_redshift, omega0=omega0, omegal0=omegal0))[0]
for i=0L, n_elements(goods_redshift)-1L do begin
    ig=where(oivar[obands[*,i], i] gt 0. AND $
             omaggies[obands[*,i], i] gt 0., ng) 
    if(ng gt 0) then begin
        mag[ig,i]=-2.5*alog10(omaggies[obands[ig,i],i])- $
          (lf_distmod(sdss_redshift[i], omega0=omega0, omegal0=omegal0))[0]+ $
          (lf_distmod(goods_redshift[i], omega0=omega0, omegal0=omegal0))[0]- $
          offset[ig,i]
        mag_ivar[ig,i]=omaggies[obands[ig,i],i]^2*oivar[obands[ig,i],i]* $
          (0.4*alog(10.))^2
    endif 
endfor

return, mag

end
