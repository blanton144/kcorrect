;+
; NAME:
;   k_model_fluxes
;
; PURPOSE:
;   Calculates the observed fluxes of a bunch of galaxies at a
;   redshift, given the coefficients and the templates, as well
;   as the filter information.
;
; CALLING SEQUENCE:
;   k_model_fluxes,
;
; INPUTS:
;   coeff    - coefficients [N_template, N_gal]
;   galaxy_z      - redshift for each galaxy (where you want to
;                   calculate flux in model) [N_gal]
;
; OPTIONAL INPUTS:
;   version   - version of eigentemplates to use
;   vpath   - version of eigentemplates to use
;
; OUTPUTS:
;   model_flux  - flux in each band
;
; OPTIONAL INPUT/OUTPUTS:
;   ematrix       - eigentemplates [N_dim, N_template]
;   filterlist    - list of files with filter information [N_band]
;   bmatrix       - orthogonal templates spanning SED space [N_lambda, N_dim]
;   lambda        - wavelengths for orthogonal templates [N_lambda]
;   rmatrix       - look up table for bmatrix and filter information 
;                   [N_z, N_dim, N_band]
;   zvals         - look up table for rmatrix [N_z]
;
; COMMENTS:
;   Outputs fluxes in maggies. Does not calculate errors.
; 
; EXAMPLES:
;   Given coeffs from a call to k_fit_coeff, create reconstructed fluxes
;   for each galaxy using default templates
; 
;   IDL> k_model_fluxes,coeffs,galaxy_z,model_flux,/default
;  
;   To K-correct all the fluxes to the same redshift (0.1):
;
;   IDL> k_model_fluxes, coeffs,replicate(0.1,n_elements(galaxy_z)), $
;        model_flux, /default
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   05-Jan-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_model_fluxes,coeffs,galaxy_z,model_flux,ematrix=ematrix,zvals=zvals,rmatrix=rmatrix,bmatrix=bmatrix,lambda=lambda,version=version,vpath=vpath,filterpath=filterpath,filterlist=filterlist

; Need at least 3 parameters
if (N_params() LT 3) then begin
    klog, 'Syntax - k_model_fluxes, coeffs, galaxy_z, model_flux, [ematrix=, zvals=, $'
    klog, '    rmatrix=, bmatrix=, lambda=, version=, vpath=, filterpath=]'
    return
endif


if(NOT keyword_set(filterpath)) then $
  filterpath=getenv('KCORRECT_DIR')+'/data/filters'

; Get bmatrix and stuff from files if necessary
if(keyword_set(version)) then begin 
    if(NOT keyword_set(vpath)) then $
      vpath=getenv('KCORRECT_DIR')+'/data/etemplates'
    k_load_ascii_table,ematrix,vpath+'/ematrix.'+version+'.dat'
    k_load_ascii_table,bmatrix,vpath+'/bmatrix.'+version+'.dat'
    k_load_ascii_table,lambda,vpath+'/lambda.'+version+'.dat'
    filtfile=vpath+'/filterlist.'+version+'.dat'
    spawn,'cat '+filtfile+' | wc -l',nfilters
    filterlist=strarr(nfilters[0])
    openr,unit,vpath+'/filterlist.'+version+'.dat',/get_lun
    readf,unit,filterlist
    close,unit
    free_lun,unit
endif

ngalaxy=long(n_elements(galaxy_z))
nt=long(n_elements(coeffs))/ngalaxy
nb=long(n_elements(ematrix))/nt
if (keyword_set(bmatrix) AND keyword_set(filterlist)  $
    AND keyword_set(lambda)) then begin
    nk=long(n_elements(filterlist))
    k_create_r,rmatrix,bmatrix,lambda,zvals,filterlist,filterpath=filterpath
endif else begin
    if (NOT keyword_set(rmatrix) OR NOT keyword_set(zvals)) then begin
        klog, 'need to specify rmatrix or bmatrix and filterlist'
    endif
endelse 
nz=long(n_elements(zvals))
nk=long(n_elements(rmatrix))/(nb*nz)

; Set source object name
soname=filepath('libkcorrect.so', $
                root_dir=getenv('KCORRECT_DIR'), subdirectory='lib')

; Call coefficient software
model_flux=dblarr(nk,ngalaxy)
retval=call_external(soname, 'idl_k_model_fluxes', double(ematrix), $
                     long(nt), double(zvals), long(nz), double(rmatrix), $
                     long(nk), long(nb), double(coeffs), $
                     double(galaxy_z), double(model_flux), long(ngalaxy) )

end
;------------------------------------------------------------------------------
