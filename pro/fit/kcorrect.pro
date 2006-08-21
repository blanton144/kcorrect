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
;                   (implies /stddev set)
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
;   /silent       - shut up
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
;   rmaggies      - [nk, ngals] reconstructed maggies in each bands
;   coeffs        - [nt, ngals] coefficients fit to each template (if maggies
;                   input are nonexistent, just use these input
;                   coeffs)
;   mass          - [ngals] model mass derived from the coeffs (current
;                   surviving stellar mass)
;   intsfh        - [ngals] total integrated star-formation history (number of
;                   stars formed total)
;   mets          - [ngals] average metallicity of currently surviving stars

;   absmag        - [nk, ngals] absolute magnitude (for missing data,
;                   substitutes model fit)

;   amivar        - [nk, ngals] absolute magnitude invvar (for missing
;                   data = 0) 
;   mtol          - [nk, ngals] model mass-to-light in each *final*
;                   bandpass (the bandpass you are kcorrecting *to*)
;                   in SOLAR UNITS!
;   lambda        - [nl+1] wavelengths for templates (to use)/(which were used)
;                   (pixel edges)
;   vmatrix       - [nl, nv] templates (to use)/(which were used)
;   rmatrix       - [nz, nv, nk] look up table for bmatrix and filter
;                   information
;   zvals         - [nz] look up redshift table for rmatrix [N_z]
;   zmin,zmax     - minimum and maximum redshifts for lookup table
;                   (default 0., 2.)
;   nz            - number of redshifts in lookup table (default 1000)
;   b300          - star-formation within last 300Myrs relative to
;                   average star-formation rate
;   b1000         - star-formation within last 1Gyrs relative to
;                   average star-formation rate
; COMMENTS:
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
              absmag=absmag, amivar=amivar, omega0=omega0, omegal0=omegal0, $
              mets=mets, b300=b300, b1000=b1000, intsfh=intsfh, silent=silent

littleh=0.7 ;; for evolution

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
      zmin=zmin,zmax=zmax,nz=nz,filterpath=filterpath, silent=silent
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
            stddev=1
            use_maggies=10.^(-0.4*use_maggies)
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

tspecfile=getenv('KCORRECT_DIR')+'/data/templates/k_nmf_derived.'+ $
  vname+'.fits'
tmass=mrdfits(tspecfile, 16, silent=silent)
tmremain=mrdfits(tspecfile, 17, silent=silent)
tmetallicity=mrdfits(tspecfile, 18, silent=silent)
tmass300=mrdfits(tspecfile, 19, silent=silent)
tmass1000=mrdfits(tspecfile, 20, silent=silent)

mrcoeffs=coeffs*(tmremain#replicate(1., n_elements(redshift)))
if(arg_present(mass)) then $
  mass=total(mrcoeffs,1)*10.^(0.4*lf_distmod(redshift, omega0=omega0, $
                                             omegal0=omegal0))
if(arg_present(intsfh)) then $
  intsfh=total(coeffs,1)*10.^(0.4*lf_distmod(redshift, omega0=omega0, $
                                             omegal0=omegal0))
smaggies=10.^(-0.4*k_solar_magnitudes(filterlist=filterlist, $
                                      band_shift=band_shift, silent=silent))
mtol=fltarr(n_elements(filterlist), n_elements(redshift))
mm=total(mrcoeffs,1)
for i=0L, n_elements(filterlist)-1L do $
  mtol[i,*]=mm/reconstruct_maggies[i,*]*smaggies[i]

; get metallicity
b300=fltarr(n_elements(redshift))
b1000=fltarr(n_elements(redshift))
mets=fltarr(n_elements(redshift))
if((size(coeffs))[0] eq 2) then begin
    for i=0L, n_elements(redshift)-1L do begin
        b300[i]=total(tmass300*coeffs[*,i])/total(tmass*coeffs[*,i])
        b1000[i]=total(tmass1000*coeffs[*,i])/total(tmass*coeffs[*,i])
        tmp_mass=total(tmremain*coeffs[*,i])
        mets[i]=total(tmremain*tmetallicity*coeffs[*,i])/tmp_mass
    endfor
endif else begin
    b300=total(tmass300*coeffs)/total(tmass*coeffs)
    b1000=total(tmass1000*coeffs)/total(tmass*coeffs)
    tmp_mass=total(tmremain*coeffs)
    mets=total(tmremain*tmetallicity*coeffs)/tmp_mass
endelse

; get kcorrection
kcorrect=reconstruct_maggies/rmaggies
kcorrect=2.5*alog10(kcorrect)

if(arg_present(absmag)) then begin
    absmag=fltarr(n_elements(filterlist), n_elements(redshift))
    amivar=fltarr(n_elements(filterlist), n_elements(redshift))
    for i=0L, n_elements(filterlist)-1L do $
      absmag[i,*]=-2.5*alog10(reconstruct_maggies[i,*])- $
      lf_distmod(redshift, omega0=omega0,omegal0=omegal0)
    if(keyword_set(use_maggies) gt 0) then begin
        for i=0L, n_elements(filterlist)-1L do begin
            ig=where(use_maggies_ivar[i,*] gt 0. AND use_maggies[i,*] gt 0., $
                     ng)
            if(ng gt 0) then begin
                absmag[i,ig]=-2.5*alog10(use_maggies[i,ig])- $
                  lf_distmod(redshift[ig], omega0=omega0, omegal0=omegal0)- $
                  kcorrect[i,ig]
                amivar[i,ig]=use_maggies[i,ig]^2*use_maggies_ivar[i,ig]* $
                  (0.4*alog(10.))^2
            endif 
        endfor
    endif
endif

if(arg_present(evolve)) then begin
; calculate the preliminaries
    if(keyword_set(evname) eq 0 AND $
       keyword_set(vname) gt 0) then $
      evname=vname+'late'
    if(NOT keyword_set(ermatrix) OR NOT keyword_set(zvals)) then begin
        if(NOT keyword_set(evmatrix) OR NOT keyword_set(lambda)) then $
          k_load_vmatrix, evmatrix, lambda, vfile=evfile, lfile=lfile, $
          vpath=vpath, vname=evname
        k_projection_table,ermatrix,evmatrix,lambda,zvals,filterlist, $ 
          zmin=zmin,zmax=zmax,nz=nz,filterpath=filterpath, silent=silent
    endif

; now calculate the bandpasses if they were measured 1.0 Gyr ago
    k_reconstruct_maggies,coeffs,replicate(band_shift,n_elements(redshift)), $
      emaggies,rmatrix=ermatrix,zvals=zvals
    emaggies=emaggies/(1.+band_shift)
    ediff=2.5*alog10(emaggies/reconstruct_maggies)

; how far back in redshift was that?
    time0=lf_z2t(band_shift, omega0=omega0, omegal0=omegal0)/littleh
    times=lf_z2t(redshift, omega0=omega0, omegal0=omegal0)/littleh
    eredshift=lf_t2z((times+0.1)*littleh > 0.1, omega0=omega0, omegal0=omegal0)
    zdiff=eredshift-redshift

; what is evolution in mags/redshift
    eper=ediff
    for i=0L, n_elements(filterlist)-1L do $
      eper[i,*]=ediff[i,*]/zdiff
    evolve=eper

; now what is the evolution correction?
    for i=0L, n_elements(filterlist)-1L do $
      evolve[i,*]=eper[i,*]*(band_shift-redshift)
endif

end
