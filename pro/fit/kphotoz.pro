;+
; NAME:
;   kphotoz
; PURPOSE:
;   Given AB maggies estimate the redshift of a galaxy
; CALLING SEQUENCE:
;   kphotoz, maggies, maggies_ivar, redshift, photoz [ , $
;        /magnitude, /stddev, lfile=, $
;        vfile=, vpath=, filterlist=, filterpath=, rmatrix=, $
;        zvals=, lambda=, vmatrix=, /sdssfix, coeffs=, chi2=, $
;        maxiter=, zmin=, zmax=, nz=, /verbose ]
; INPUTS:
;   maggies    - AB maggies of galaxies [N_band, N_gal] (magnitudes if
;                /magnitude set)
;   maggies_ivar - inverse variance in maggies (magnitudes if
;                  /magnitude set, std. dev. if /stddev set)
;   redshift      - redshifts of galaxies [N_gal]
; OPTIONAL INPUTS:
;   magnitude     - set if input and output in -2.5 log_10(maggies)
;   stddev        - maggies_ivar actual contains standard dev.
;   lfile    - wavelength file for vmatrix [default lambda.default.dat]
;   vfile         - vmatrix file [default vmatrix.default.dat]
;   vpath         - path to templates [default $KCORRECT_DIR/data/templates]
;   filterlist    - list of filters [default
;                                    ['sdss_u0.par', 'sdss_g0.par',
;                                     'sdss_r0.par', 'sdss_i0.par',
;                                     'sdss_z0.par']]
;   filterpath    - path to filters [default $KCORRECT_DIR/data/filters]
;   /sdssfix      - uses k_sdssfix to "fix" input SDSS magnitudes and 
;                   standard deviations (treats as if /magnitude and
;                   /stddev are also set)
;   /verbose      - call k_fit_nonneg verbosely
;   maxiter       - maximum number of iterations for fit [default
;                   3000]
; OUTPUTS:
;   photoz     - photometric redshift
;   coeffs     - coefficients fit to each template
;   chi2       - chi^2 of fit
; OPTIONAL INPUT/OUTPUTS:
;   lambda        - wavelengths for templates (to use)/(which were used)
;   vmatrix       - templates (to use)/(which were used)
;   rmatrix       - look up table for bmatrix and filter information 
;                   [N_z, N_dim, N_band]
;   zvals         - look up redshift table for rmatrix [N_z]
; COMMENTS:
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   04-Jun-2003  converted from kcorrect MRB, NYU
;-
;------------------------------------------------------------------------------
pro kphotoz, maggies, maggies_ivar, photoz, $
             magnitude=magnitude, stddev=stddev, $
             lfile=lfile, vfile=vfile, vpath=vpath, $
             filterlist=filterlist, filterpath=filterpath, $
             rmatrix=rmatrix, zvals=zvals, lambda=lambda, $
             vmatrix=vmatrix, sdssfix=sdssfix, coeffs=coeffs, $
             chi2=chi2, maxiter=maxiter, verbose=verbose

; Need at least 6 parameters
if (N_params() LT 3) then begin
    print, 'Syntax - kphotoz, maggies, maggies_ivar, photoz [ , $'
    print, '             /magnitude, /stddev, lfile=, $'
    print, '             vfile=, vpath=, filterlist=, filterpath=, rmatrix=, $'
    print, '             zvals=, lambda=, vmatrix=, coeffs=, /verbose /sdssfix ]'
    return
endif

ndim=(size(maggies))[0]
if(ndim eq 1) then begin
    nk=1
endif else begin
    nk=(size(maggies,/dimens))[0]
endelse
ngalaxy=long(n_elements(maggies))/nk
if(NOT keyword_set(filterlist)) then $
  filterlist=['sdss_u0.par','sdss_g0.par','sdss_r0.par','sdss_i0.par', $
              'sdss_z0.par']
if(NOT keyword_set(maxiter)) then maxiter=3000L

; Fix SDSS mags if desired
use_maggies=maggies
use_maggies_ivar=maggies_ivar
if(keyword_set(sdssfix)) then begin
    magnitude=1
    stddev=1
    tmp_mag=maggies
    tmp_magerr=maggies_ivar
    k_sdssfix,tmp_mag,tmp_magerr,use_maggies,use_maggies_ivar
endif else begin
    if(keyword_set(magnitude)) then begin
        use_maggies=10.^(-0.4*use_maggies)
        if(keyword_set(stddev)) then $
          use_maggies_ivar=1./(use_maggies*0.4*alog(10.)*use_maggies_ivar)^2
    endif else begin
        if(keyword_set(stddev)) then $
          use_maggies_ivar=1./use_maggies_ivar^2
    endelse
endelse 

; Calculate coeffs
if(NOT keyword_set(rmatrix) OR NOT keyword_set(zvals)) then $
  k_load_vmatrix, vmatrix, lambda, vfile=vfile, lfile=lfile, $
  vpath=vpath
photoz=k_fit_photoz(use_maggies,use_maggies_ivar,vmatrix,lambda, $
                    filterlist=filterlist,chi2=chi2, $
                    rmatrix=rmatrix,zvals=zvals,maxiter=maxiter, zmin=zmin, $
                    zmax=zmax, nz=nz, verbose=verbose)

end
