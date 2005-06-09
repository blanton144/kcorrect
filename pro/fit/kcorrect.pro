;+
; NAME:
;   kcorrect
; PURPOSE:
;   Given a set of AB maggies, returns the K-correction for each band.
; CALLING SEQUENCE:
;   kcorrect, maggies, maggies_ivar, redshift, kcorrect [ , $
;        band_shift=, /magnitude, /stddev, lfile=, absmag=, $
;        vfile=, vpath=, filterlist=, filterpath=, rmatrix=, $
;        zvals=, lambda=, vmatrix=, /sdssfix, /abfix, minerrors=, $
;        coeffs=, chi2=, maxiter=, zmin=, zmax=, nz=, /verbose ]
; INPUTS:
;   maggies    - [nk, ngals] AB maggies of galaxies (magnitudes if
;                /magnitude set, asinh magnitudes if /sdssfix set)
;   maggies_ivar - [nk, ngals] inverse variance in maggies (magnitudes if
;                  /magnitude set, std. dev. if /stddev set; if
;                  /sdssfix set, is std. dev. of asinh magnitudes)
;   redshift      - [ngals] redshifts of galaxies 
; OPTIONAL INPUTS:
;   band_shift    - blueshift of bandpasses to apply (to get ^{z}b
;                   type bands) [default 0]
;   magnitude     - set if input and output in -2.5 log_10(maggies)
;   stddev        - maggies_ivar actual contains standard dev.
;   minerrors     - [nk] add this set of errors (in magnitude units)
;                   in quadrature to all uncertainties
;   vname         - name of fit to us (default 'default')
;   lfile         - wavelength file for vmatrix [default
;                   lambda.[vname].dat]
;   vfile         - vmatrix file [default vmatrix.[vname].dat]
;   vpath         - path to templates [default $KCORRECT_DIR/data/templates]
;   filterlist    - [nk] list of filters [default
;                                         ['sdss_u0.par', 'sdss_g0.par',
;                                          'sdss_r0.par', 'sdss_i0.par',
;                                          'sdss_z0.par']]
;   filterpath    - path to filters [default $KCORRECT_DIR/data/filters]
;   /verbose      - call k_fit_nonneg verbosely
;   maxiter       - maximum number of iterations for fit [default 3000]
;   omega0, omegal0 - cosmological parameters for calculating distance
;                     moduli [default 0.3, 0.7]
; DEPRECATED INPUTS:
;   /abfix        - uses k_abfix to fix input SDSS maggies to AB
;                   (better to just call with 'sdss_kcorrect')
;   /sdssfix      - uses k_sdssfix to "fix" input SDSS asinh magnitudes and 
;                   standard deviations; treats as if /abfix, and
;                   minerrors=[0.05,0.02,0.02,0.02,0.03] are all set
;                   (better to just call with 'sdss_kcorrect')
; OUTPUTS:
;   kcorrect   - [nk, ngals] K-corrections satisfying
;                   m = M + DM(z) + K(z)
;                based on the best fit sum of templates
;   chi2       - chi^2 of fit
;   rmaggies      - reconstructed maggies from the fit
; OPTIONAL INPUT/OUTPUTS:
;   coeffs        - coefficients fit to each template (if maggies
;                   input are nonexistent, just use these input
;                   coeffs)
;   mass          - model mass derived from the coeffs
;   absmag        - absolute magnitude (for missing data, substitutes
;                   model fit)
;   amivar        - absolute magnitude invvar (for missing data = 0)
;   mtol          - model mass-to-light in each *final* bandpass (the
;                   bandpass you are kcorrecting *to*) in SOLAR UNITS!
;   lambda        - [nl+1] wavelengths for templates (to use)/(which were used)
;                   (pixel edges)
;   vmatrix       - [nl, nv] templates (to use)/(which were used)
;   rmatrix       - look up table for bmatrix and filter information 
;                   [nz, nv, nk]
;   zvals         - look up redshift table for rmatrix [N_z]
;   zmin,zmax     - minimum and maximum redshifts for lookup table
;                   (default 0., 2.)
;   nz            - number of redshifts in lookup table (default 1000)
; COMMENTS:
;   For v4_0b templates and later, coefficients are in units of:
;     1 solar mass / (D/10pc)^2
;   That is, sum the coefficients and multiply by (D/10pc)^2 to get
;   masses. (In fact, for Omega0=0.3 and OmegaL0=0.7, this is what the
;   "mass" keyword returns).
;
;   If you just want to do SDSS kcorrections, it is better to use the
;   wrapper 'sdss_kcorrect'. 
; 
;   If you want to do GALEX kcorrections (perhaps including matched
;   SDSS data) use the wrapper 'galex_kcorrect'.
; 
;   If you want to do 2MASS kcorrections (perhaps including matched
;   SDSS data) use the wrapper 'twomass_kcorrect'.
; 
;   If you want to do DEEP kcorrections use the wrapper 'deep_kcorrect'.
;
;   Allows the user to shift the bandpasses by a factor band_shift. 
;   If no band_shift is specified, the K-correction is to z=0.
; 
;   Defaults to SDSS filters.  
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
              chi2=chi2, maxiter=maxiter, verbose=verbose, nz=nz, $
              zmin=zmin, zmax=zmax, abfix=abfix, minerrors=minerrors, $
              rmaggies=rmaggies, mass=mass, mtol=mtol, vname=vname, $
              absmag=absmag, amivar=amivar, omega0=omega0, omegal0=omegal0

; Need at least 6 parameters
if (N_params() LT 4) then begin
    doc_library, 'kcorrect'
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
    if(NOT keyword_set(vmatrix) OR NOT keyword_set(lambda)) then $
      k_load_vmatrix, vmatrix, lambda, vfile=vfile, lfile=lfile, $
      vpath=vpath, vname=vname
    k_projection_table,rmatrix,vmatrix,lambda,zvals,filterlist, $ 
      zmin=zmin,zmax=zmax,nz=nz,filterpath=filterpath
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
        if(keyword_set(abfix)) then $
          k_abfix, use_maggies, use_maggies_ivar
        if(keyword_set(minerrors)) then $
          k_minerror, use_maggies, use_maggies_ivar, minerrors
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
k_reconstruct_maggies,coeffs,redshift,rmaggies,rmatrix=rmatrix,zvals=zvals

if(arg_present(mass)) then $
  mass=total(coeffs,1)*10.^(0.4*lf_distmod(redshift, omega0=omega0, $
                                           omegal0=omegal0))
smaggies=10.^(-0.4*k_solar_magnitudes(filterlist=filterlist, $
                                      band_shift=band_shift))
mtol=fltarr(n_elements(filterlist), n_elements(redshift))
mm=total(coeffs,1)
for i=0L, n_elements(filterlist)-1L do $
  mtol[i,*]=mm/reconstruct_maggies[i,*]*smaggies[i]

; get kcorrection
kcorrect=reconstruct_maggies/rmaggies
kcorrect=2.5*alog10(kcorrect)

if(arg_present(absmag)) then begin
    absmag=fltarr(n_elements(filterlist), n_elements(redshift))
    amivar=fltarr(n_elements(filterlist), n_elements(redshift))
    for i=0L, n_elements(filterlist)-1L do $
      absmag[i,*]=-2.5*alog10(reconstruct_maggies[i,*])- $
      lf_distmod(redshift, omega0=omega0,omegal0=omegal0)
    for i=0L, n_elements(filterlist)-1L do begin
        ig=where(use_maggies_ivar[i,*] gt 0. AND use_maggies[i,*] gt 0., ng)
        if(ng gt 0) then begin
            absmag[i,ig]=-2.5*alog10(use_maggies[i,ig])- $
              lf_distmod(redshift[ig], omega0=omega0, omegal0=omegal0)- $
              kcorrect[i,ig]
            amivar[i,ig]=use_maggies[i,ig]^2*use_maggies_ivar[i,ig]* $
              (0.4*alog(10.))^2
        endif 
    endfor
endif

end
