;+
; NAME:
;   k_fit_coeffs
;
; PURPOSE:
;   Fit coefficients to eigentemplates given photometry and redshift of 
;   a galaxy.
;
; CALLING SEQUENCE:
;   k_fit_coeffs, galaxy_maggies, galaxy_invvar, galaxy_z, coeff, 
;      [ematrix=, zvals=, filterlist=, bmatrix=, lambda=, rmatrix=, /default]
;
; INPUTS:
;   galaxy_maggies   - maggies in each band for each galaxy [N_band, N_gal]
;   galaxy_invvar - errors in each band for each galaxy [N_band, N_gal]
;   galaxy_z      - redshift for each galaxy [N_gal]
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
;
; OPTIONAL INPUT/OUTPUTS:
;   default  - set to use the default templates
;
; COMMENTS:
;   galaxy_maggies is in maggies (f=10.^{-0.4*mag}). galaxy_invvar is in 
;   maggies^{-2} ((f*0.4*ln(10)*magerr)^{-2}).
;
;   If filterlist, bmatrix, and lambda are specified, rmatrix is created 
;   (and replaced if asked by the user)
;
; EXAMPLES:
;   To get the coefficients using the standard templates:
; 
;   IDL> k_fit_coeffs,galaxy_maggies,galaxy_invvar,galaxy_z, coeff, /default
; 
;   Then you can pass "coeff" into k_reconstruct_maggies
;
; BUGS:
;   Will fail if N_dim is unity.
;
; PROCEDURES CALLED:
;   k_load_ascii_table
;   k_create_r
;   Dynamic link to idl_k_fit_coeffs.c in libkcorrect.so
;
; REVISION HISTORY:
;   04-Jan-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_fit_coeffs, galaxy_maggies, galaxy_invvar, galaxy_z, coeff, ematrix=ematrix, zvals=zvals, filterlist=filterlist, bmatrix=bmatrix, lambda=lambda, rmatrix=rmatrix, version=version, vpath=vpath, filterpath=filterpath

; Need at least 6 parameters
if (N_params() LT 4) then begin
    klog, 'Syntax - k_fit_coeffs, galaxy_maggies, galaxy_invvar, galaxy_z, coeff, $'  
    klog, '         [ematrix=, zvals=, filterlist=, bmatrix=, lambda=, rmatrix=, $'
    klog, '          version=, vpath=, filterpath=]'
    return
endif

if(NOT keyword_set(filterpath)) then $j
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
    nk=long(nfilters[0])-1l
    filterlist=strarr(nk)
    openr,unit,vpath+'/filterlist.'+version+'.dat',/get_lun
    readf,unit,nk
    readf,unit,filterlist
    close,unit
    free_lun,unit
endif

; Set source object name
soname=filepath('libkcorrect.so', $
                root_dir=getenv('KCORRECT_DIR'), subdirectory='lib')

; Determine dimensions
ngalaxy=long(n_elements(galaxy_z))
nk=long(n_elements(galaxy_maggies)/ngalaxy)
if (keyword_set(bmatrix) AND keyword_set(filterlist) AND keyword_set(lambda)) $
  then begin
    k_create_r,rmatrix,bmatrix,lambda,zvals,filterlist,filterpath=filterpath
endif else begin
    if (NOT keyword_set(rmatrix)) then begin
        klog, 'need to specify rmatrix or bmatrix and filterlist'
    endif
endelse 
nz=long(n_elements(zvals))
nb=long(n_elements(rmatrix)/(nz*nk)) 	 
nt=long(n_elements(ematrix)/nb)

; Call coefficient software
coeff=dblarr(nt,ngalaxy)
retval=call_external(soname, 'idl_k_fit_coeffs', double(ematrix), $
                     long(nt), double(zvals), long(nz), double(rmatrix), $
                     long(nk), long(nb), double(coeff), $
                     double(galaxy_maggies), double(galaxy_invvar), $
                     double(galaxy_z), long(ngalaxy) )

end
;------------------------------------------------------------------------------
