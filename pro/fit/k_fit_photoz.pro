;+
; NAME:
;   k_fit_photoz
; PURPOSE:
;   Fit photometric redshift using nonnegative sum of templates
; CALLING SEQUENCE:
;   photoz= k_fit_photoz(maggies, maggies_ivar, vmatrix, lambda, $
;      filterlist= [, chi2=, rmatrix=, zvals=, maxiter=, filterpath=, $ 
;      zmin=, zmax=, nz=, band_shift=, /verbose])
; INPUTS:
;   maggies - fluxes in each band for each galaxy
;   maggies_var - inverse variance in each band for each galaxy
;   vmatrix - templates
;   lambda - pixel edges for all templates
;   filterlist - list of filternames
; OPTIONAL INPUTS:
;   filterpath - path to look for filters on
;   maxiter - maximum number of iterations in fit
;   z[min|max] - limits for redshift to use in making lookup table
;   nz - number of redshift rows in lookup table
;   band_shift - shifted bands if desired
;   zprior/lprior - grid of redshift and -2ln(prior) values to apply
;   /verbose - verbose output
; OUTPUTS:
;   photoz - estimated redshift
;   chi2 - chi^2 value for each fit
;   coeffs - coeffs for each fit
;   rmatrix - projection table used
;   zvals - redshift list for projection table used
;   niter - number of iterations for last fit
; REVISION HISTORY:
;   01-May-2003  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function k_fit_photoz, maggies, maggies_ivar, vmatrix, lambda, $
                       filterlist=filterlist, chi2=chi2, lprior=lprior, $
                       zprior=zprior, rmatrix=rmatrix, zvals=zvals, $
                       maxiter=maxiter, filterpath=filterpath, zmin=zmin, $
                       zmax=zmax, nz=nz, band_shift=band_shift, $
                       tolerance=tolerance, verbose=verbose, niter=niter, $
                       coeffs=coeffs
                  
if(n_params() ne 4) then begin
    print,'Syntax - photoz= k_fit_photoz(maggies, maggies_ivar, vmatrix, lambda, $'
    print,'          filterlist= [, chi2=, rmatrix=, zvals=, maxiter=, filterpath=, ]'
    print,'          zmin=, zmax=, nz=, band_shift=, coeffs=, lprior=, zprior=, $' 
    print,'          /verbose ])'
    return,-1
endif 

if(NOT keyword_set(verbose)) then verbose=0L
if(NOT keyword_set(maxiter)) then maxiter=50000
if(NOT keyword_set(tolerance)) then tolerance=1.e-6

k_load_filters,filterlist,filter_n,filter_lambda,filter_pass, $
  filterpath=filterpath

; Set source object name
soname=filepath('libkcorrect.'+kcorrect_so_ext(), $
                root_dir=getenv('KCORRECT_DIR'), subdirectory='lib')

; Create rmatrix if necessary
nk=n_elements(filterlist)
ngalaxy=long(n_elements(maggies))/nk
if (keyword_set(vmatrix) AND keyword_set(filterlist) AND keyword_set(lambda)) $
  then begin
    k_projection_table,rmatrix,vmatrix,lambda,zvals,filterlist, $ 
      zmin=zmin,zmax=zmax,nz=nz,band_shift=band_shift,filterpath=filterpath
endif else begin
    if (NOT keyword_set(rmatrix) and NOT keyword_set(zvals)) then begin
        klog, 'ERROR: need to specify (rmatrix,zvals) or (vmatrix,lambda,filterlist)'
    endif
endelse 
nz=long(n_elements(zvals))
nv=long(n_elements(rmatrix)/(nz*nk)) 	 

nprior=n_elements(zprior)
if(n_elements(lprior) eq 0 or n_elements(zprior) eq 0) then begin
    nprior=2L
    lprior=fltarr(nprior)
    zprior=[0.,1000.]
endif

; Set rmatrix
photoz=fltarr(ngalaxy)
coeffs=fltarr(nv,ngalaxy)
chi2=fltarr(ngalaxy)
niter=0L
retval=call_external(soname, 'idl_k_fit_photoz', float(photoz), $
                     float(coeffs), float(rmatrix), long(nk), long(nv), $
                     float(lprior), float(zprior), long(nprior), $
                     float(zvals), long(nz), float(maggies), $
                     float(maggies_ivar), long(ngalaxy), float(tolerance), $
                     long(maxiter), long(niter),float(chi2),long(verbose)) 

return,photoz

end
;------------------------------------------------------------------------------
