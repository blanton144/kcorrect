;+
; NAME:
;   kcorrect
; PURPOSE:
;   Given a set of AB maggies, returns the K-correction for each
;   band. Allows the user to shift the bandpasses by a factor
;   band_shift. If no band_shift is specified, the K-correction is to
;   z=0.
; CALLING SEQUENCE:
;   kcorrect, maggies, maggies_ivar, redshift, kcorrect [ , $
;        band_shift=, /magnitude, /stddev, lfile=, $
;        vfile=, vpath=, filterlist=, filterpath=, rmatrix=, $
;        zvals=, lambda=, vmatrix=, /sdssfix, coeffs=, chi2=, $
;        maxiter=, zmin=, zmax=, nz=, /verbose ]
;      
; INPUTS:
;   maggies    - [nk, ngals] AB maggies of galaxies (magnitudes if
;                /magnitude or /sdssfix set)
;   maggies_ivar - [nk, ngals] inverse variance in maggies (magnitudes if
;                  /magnitude set, std. dev. if /stddev or /sdssifx set)
;   redshift      - [ngals] redshifts of galaxies 
;
; OPTIONAL INPUTS:
;   band_shift    - blueshift of bandpasses to apply (to get ^{z}b
;                   type bands) [default 0]
;   magnitude     - set if input and output in -2.5 log_10(maggies)
;   stddev        - maggies_ivar actual contains standard dev.
;   lfile         - wavelength file for vmatrix [default lambda.default.dat]
;   vfile         - vmatrix file [default vmatrix.default.dat]
;   vpath         - path to templates [default $KCORRECT_DIR/data/templates]
;   filterlist    - [nk] list of filters [default
;                                         ['sdss_u0.par', 'sdss_g0.par',
;                                          'sdss_r0.par', 'sdss_i0.par',
;                                          'sdss_z0.par']]
;   filterpath    - path to filters [default $KCORRECT_DIR/data/filters]
;   /sdssfix      - uses k_sdssfix to "fix" input SDSS magnitudes and 
;                   standard deviations; this applies the AB fix to
;                   the SDSS magnitudes; it treats as if /magnitude and
;                   /stddev are also set 
;   /verbose      - call k_fit_nonneg verbosely
;   maxiter       - maximum number of iterations for fit [default 3000]
; OUTPUTS:
;   kcorrect   - [nk, ngals] K-corrections satisfying
;                   m = M + DM(z) + K(z)
;                based on the best fit sum of templates
;   chi2       - chi^2 of fit
;
; OPTIONAL INPUT/OUTPUTS:
;   coeffs        - coefficients fit to each template (if maggies
;                   input are nonexistent, just use these input
;                   coeffs)
;   lambda        - [nl+1] wavelengths for templates (to use)/(which were used)
;                   (pixel edges)
;   vmatrix       - [nl, nv] templates (to use)/(which were used)
;   rmatrix       - look up table for bmatrix and filter information 
;                   [nz, nv, nk]
;   zvals         - look up redshift table for rmatrix [N_z]
;
; COMMENTS:
;
;   Be careful when sending SDSS tsObj outputs directly into this program.
;   Eg. occasionally the magnitudes or the errors have crazy values,
;   such as -9999. The normal "garbage in, garbage out" rules
;   apply. However, I have supplied the /sdssfix flag, which does a
;   reasonable job in most cases of identifying problem cases and
;   doing something OK about it. 
;
; EXAMPLES:
;
; BUGS:
;   Many inputs NOT commented, including zmin and zmax.
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   24-Jan-2002  Translated to IDL by Mike Blanton, NYU
;   19-Jul-2002  Major bug fix (pointed out by I. Baldry) MRB, NYU
;   02-Jun-2003  Updated to new v3_0 standards MRB, NYU
;-
;------------------------------------------------------------------------------
pro kcorrect, maggies, maggies_ivar, redshift, kcorrect, $
              band_shift=band_shift, magnitude=magnitude, stddev=stddev, $
              lfile=lfile, vfile=vfile, vpath=vpath, $
              filterlist=filterlist, filterpath=filterpath, $
              rmatrix=rmatrix, zvals=zvals, lambda=lambda, $
              vmatrix=vmatrix, sdssfix=sdssfix, coeffs=coeffs, $
              chi2=chi2, maxiter=maxiter, verbose=verbose, nz=nz

; Need at least 6 parameters
if (N_params() LT 4) then begin
    print, 'Syntax - kcorrect, maggies, maggies_ivar, redshift, kcorrect [ , $'
    print, '             band_shift=, /magnitude, /stddev, lfile=, $'
    print, '             vfile=, vpath=, filterlist=, filterpath=, rmatrix=, $'
    print, '             zvals=, lambda=, vmatrix=, coeffs=, /verbose /sdssfix ]'
    return
endif

ngalaxy=long(n_elements(redshift))
if(NOT keyword_set(filterlist)) then $
  filterlist=['sdss_u0.par','sdss_g0.par','sdss_r0.par','sdss_i0.par', $
              'sdss_z0.par']
if(NOT keyword_set(band_shift)) then band_shift=0.
if(NOT keyword_set(maxiter)) then maxiter=3000L

; calculate the preliminaries
if(NOT keyword_set(rmatrix) OR NOT keyword_set(zvals)) then begin
    k_load_vmatrix, vmatrix, lambda, vfile=vfile, lfile=lfile, $
      vpath=vpath
    k_projection_table,rmatrix,vmatrix,lambda,zvals,filterlist, $ 
      zmin=zmin,zmax=zmax,nz=nz,band_shift=band_shift,filterpath=filterpath
endif

; Calculate the coefficients if needed
if(n_elements(maggies) gt 0) then begin 
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
              use_maggies_ivar= $
              1./(use_maggies*0.4*alog(10.)*use_maggies_ivar)^2
        endif else begin
            if(keyword_set(stddev)) then $
              use_maggies_ivar=1./use_maggies_ivar^2
        endelse
    endelse 
    
; Calculate coeffs
    coeffs=k_fit_nonneg(use_maggies,use_maggies_ivar,redshift=redshift, $
                        chi2=chi2,rmatrix=rmatrix,zvals=zvals, $
                        maxiter=maxiter,zmin=zmin,zmax=zmax,nz=nz, $
                        verbose=verbose)
endif else begin
    if(n_elements(coeffs) eq 0) then begin
        splog,'Must set either maggies or coeffs'
        return
    endif
endelse
    
; Reconstruct the magnitudes as observed and in the rest frame
k_reconstruct_maggies,coeffs,replicate(band_shift,n_elements(redshift)), $
  reconstruct_maggies,rmatrix=rmatrix,zvals=zvals
reconstruct_maggies=reconstruct_maggies/(1.+band_shift)
k_reconstruct_maggies,coeffs,redshift,reconstruct_maggies0, $
  rmatrix=rmatrix,zvals=zvals

; get kcorrection
kcorrect=reconstruct_maggies/reconstruct_maggies0
kcorrect=2.5*alog10(kcorrect)

end
