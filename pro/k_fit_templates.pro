;+
; NAME:
;   k_fit_templates
;
; PURPOSE:
;   Fit eigentemplates given photometry and redshift of a set of galaxies, as
;   well as their coefficients given a previous set of eigentemplates
;
; CALLING SEQUENCE:
;   k_fit_templates, galaxy_flux, galaxy_invvar, galaxy_z, coeff, $
;      ematrix, zvals, [filterlist=, bmatrix=, lambda= $
;      rmatrix=, dmatrix=, galaxy_clip=]
;
; INPUTS:
;   galaxy_flux   - flux in each band for each galaxy [N_band, N_gal]
;   galaxy_invvar - errors in each band for each galaxy [N_band, N_gal]
;   galaxy_z      - redshift for each galaxy [N_gal]
;   galaxy_clip   - 1 to exclude galaxy from fit, 0 to include
;   coeff         - coefficients [N_template, N_gal]
;   filterlist    - list of files with filter information [N_band]
;   bmatrix       - orthogonal templates spanning SED space [N_lambda, N_dim]
;   lambda        - wavelengths for orthogonal templates [N_lambda]
;   rmatrix       - look up table for bmatrix and filter information 
;                   [N_z, N_dim, N_band]
;   dmatrix       - another lookup table for fitting the templates
;   /initialize_dmatrix  - reset dmatrix
;   zvals         - look up table for rmatrix [N_z]
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;   ematrix       - output eigentemplates [N_dim, N_template]
;
; OPTIONAL INPUT/OUTPUTS:
;
; COMMENTS:
;   If filterlist and bmatrix are specified, rmatrix is created (and replaced
;   if asked by the user). If dmatrix is specified, and initialize_dmatrix
;   is set, dmatrix is set and returned. If dmatrix is specified, but 
;   initialize_dmatrix is not set, the input dmatrix is used without being
;   reset.
;
; EXAMPLES:
;
; BUGS:
;   Will fail if N_dim is unity.
;
; PROCEDURES CALLED:
;   k_create_r
;   Dynamic link to idl_k_fit_templates.c
;
; REVISION HISTORY:
;   05-Jan-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_fit_templates, galaxy_flux, galaxy_invvar, galaxy_z, coeff, ematrix, zvals, filterlist=filterlist, bmatrix=bmatrix, lambda=lambda, rmatrix=rmatrix, dmatrix=dmatrix, galaxy_clip=galaxy_clip

; Need at least 6 parameters
if (N_params() LT 6) then begin
    klog, 'Syntax - k_fit_templates, galaxy_flux, galaxy_invvar, galaxy_z, coeff, $'
    klog, '   ematrix, zvals, [filterlist=, bmatrix=, lambda=, rmatrix=, dmatrix=, $'
    klog, '   galaxy_clip=]'
    return
endif

; Set source object name
soname=filepath('libkcorrect.so', $
                root_dir=getenv('KCORRECT_DIR'), subdirectory='lib')

; Determine dimensions
ngalaxy=long(n_elements(galaxy_z))
nk=long(n_elements(galaxy_flux)/ngalaxy)
nz=long(n_elements(zvals))
if (keyword_set(bmatrix) AND keyword_set(filterlist)) then begin
    k_create_r,rmatrix,bmatrix,lambda,zvals,filterlist
endif else begin
    if (NOT keyword_set(rmatrix)) then begin
        klog, 'need to specify rmatrix or bmatrix and filterlist'
    endif
endelse 
nb=long(n_elements(rmatrix)/(nz*nk)) 	 
nt=long(n_elements(ematrix)/nb)

; Set defaults
if (NOT keyword_set(galaxy_clip)) then begin
    galaxy_clip=lonarr(ngalaxy)
endif
if (NOT keyword_set(dmatrix)) then begin
    dmatrix=dblarr(nb,nk,ngalaxy)
    initialized_dmatrix=0l
endif else begin
    if (keyword_set(initialize_dmatrix)) then begin
        initialized_dmatrix=0l
    endif else begin
        initialized_dmatrix=1l
    endelse
endelse

; Call coefficient software
ematrix=dblarr(nb,nt)
retval=call_external(soname, 'idl_k_fit_templates', double(ematrix), $
                     long(nt), double(zvals), long(nz), double(rmatrix), $
                     long(nk), long(nb), double(coeff), $
                     double(galaxy_flux), double(galaxy_invvar), $
                     double(galaxy_z), long(galaxy_clip), long(ngalaxy), $
                     double(dmatrix), long(initialized_dmatrix))

end
;------------------------------------------------------------------------------
