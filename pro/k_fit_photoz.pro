;+
; NAME:
;   k_fit_photoz
;
; PURPOSE:
;   Fit redshifts and coefficients to eigentemplates given photometry
;   and redshift of a galaxy
;
; CALLING SEQUENCE:
;   k_fit_coeffs, galaxy_flux, galaxy_invvar, galaxy_z, coeff, ematrix, zvals,$
;      [filterlist=filterlist, bmatrix=bmatrix, lambda=lambda $
;       | rmatrix=rmatrix]
;
; INPUTS:
;   galaxy_flux   - flux in each band for each galaxy [N_band, N_gal]
;   galaxy_invvar - errors in each band for each galaxy [N_band, N_gal]
;   ematrix       - eigentemplates [N_dim, N_template]
;   filterlist    - list of files with filter information [N_band]
;   bmatrix       - orthogonal templates spanning SED space [N_lambda, N_dim]
;   lambda        - wavelengths for orthogonal templates [N_lambda]
;   rmatrix       - look up table for bmatrix and filter information 
;                   [N_z, N_dim, N_band]
;   zvals         - look up table for rmatrix [N_z]
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;   coeff    - output coefficients [N_template, N_gal]
;   galaxy_photoz      - estimated redshift for each galaxy [N_gal]
;
; OPTIONAL INPUT/OUTPUTS:
;
; COMMENTS:
;   If filterlist and bmatrix are specified, rmatrix is created (and replaced
;   if asked by the user)
;
; EXAMPLES:
;
; BUGS:
;   This is not completed yet and will not work.
;
; PROCEDURES CALLED:
;   Dynamic link to fitCoeffs.c
;
; REVISION HISTORY:
;   13-Jan-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_fit_photoz, galaxy_flux, galaxy_invvar, galaxy_photoz, coeff, ematrix, zvals, filterlist=filterlist, bmatrix=bmatrix, lambda=lambda, rmatrix=rmatrix

; Need at least 6 parameters
if (N_params() LT 6) then begin
    klog, 'Syntax - k_fit_photoz, galaxy_flux, galaxy_invvar, galaxy_photoz, coeff, ' 
    klog, 'ematrix, zvals, [filterlist=filterlist, bmatrix=bmatrix, lambda=lambda'
    klog, '                 | rmatrix=rmatrix]'
    return
endif

; Set source object name
soname=filepath('libkcorrect.so', $
                root_dir=getenv('KCORRECT_DIR'), subdirectory='lib')

; Determine dimensions
nk=long((size(galaxy_flux))[1])
ngalaxy=long(n_elements(galaxy_flux)/nk)
if (keyword_set(bmatrix) AND keyword_set(filterlist) $
    AND keyword_set(lambda)) then begin
    k_create_r,rmatrix,bmatrix,lambda,zvals,filterlist
endif else begin
    if (NOT keyword_set(rmatrix)) then begin
        klog, 'need to specify rmatrix or bmatrix, lambda, and filterlist'
    endif
endelse 
nz=long(n_elements(zvals))
nb=long(n_elements(rmatrix)/(nz*nk)) 	 
nt=long(n_elements(ematrix)/nb)

; Call coefficient software
coeff=dblarr(nt,ngalaxy)
galaxy_photoz=dblarr(ngalaxy)
retval=call_external(soname, 'idl_k_fit_photoz', double(ematrix), $
                     long(nt), double(zvals), long(nz), double(rmatrix), $
                     long(nk), long(nb), double(coeff), $
                     double(galaxy_flux), double(galaxy_invvar), $
                     double(galaxy_photoz), long(ngalaxy) )

end
;------------------------------------------------------------------------------
