;+
; NAME:
;   k_create_r
;
; PURPOSE:
;   Create lookup table for calculating SED fits. This tabulates the
;   projection of each basis element b onto each filter k, as a
;   function of redshift. You only have to perform this projection
;   once, since every spectrum we will deal with will be a linear
;   combination of the basis elements b. To get a particular redshift,
;   you interpolate the rmatrix in the z direction.
;
; CALLING SEQUENCE:
;   k_create_r, rmatrix, bmatrix, lambda, zvals, filterlist, $
;      [zmax=, zmin=, nz=]
;
; INPUTS:
;   bmatrix       - orthogonal templates spanning SED space [N_lambda, N_dim]
;   lambda        - wavelengths for orthogonal templates [N_lambda]
;   filterlist    - list of files with filter information [N_band]
;
; OPTIONAL INPUTS:
;   filterpath    - path for filters (default '$KCORRECT_DIR/data/filters')
;
; KEYWORDS:
;   system - Magnitude system: one of the predefined list: ['AB'] or a
;            filename containing g(nu) [default to AB]
;
; OUTPUTS:
;   rmatrix       - look up table for bmatrix and filter information 
;                   [N_z, N_dim, N_band]
;   zvals         - look up table for rmatrix [N_z]
;
; OPTIONAL INPUT/OUTPUTS:
;   zmin, zmax, nz  - settings for setting zvals
;
; COMMENTS:
;   Keep in mind that this only creates the r matrix for a specific
;   redshift range; e.g. you have to change the defaults in order to
;   consider redshifts greater than unity.
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   k_load_filters
;   Dynamic link to idl_k_create_r.c in libkcorrect.so
;
; REVISION HISTORY:
;   05-Jan-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_create_r, rmatrix, bmatrix, lambda, zvals, filterlist, zmin=zmin, $
                zmax=zmax, nz=nz, filterpath=filterpath, system=system

klog,'Creating rmatrix ...'

; Set defaults
if (NOT keyword_set(zmin)) then zmin=1.0e-4
if (NOT keyword_set(zmax)) then zmax=1.0e-0
if (NOT keyword_set(nz)) then nz=1000l
if (NOT keyword_set(system)) then system='AB'
if(NOT keyword_set(filterpath)) then $
  filterpath=getenv('KCORRECT_DIR')+'/data/filters'

; set g(nu)
if(system eq 'AB') then begin
    ; glambda is specified as 10.^{glambdaexp} * glambda
    glambdaexp=0.
    i=lindgen(n_elements(lambda)-1L)
    ip1=i+1
    glambda=2.99792*3.631e-2/(0.5*(lambda[i]+lambda[ip1]))^2
endif else begin
    if(file_test(system)) then begin
        k_load_ascii_table,glambda,system
        glambdaexp=0.
        glambda=glambda*10.^(-glambdaexp)
    endif else begin
        klog,'No such magnitude system: '+system
        return
    endelse
endelse

; Set zvals 
zvals=zmin+(zmax-zmin)*(dindgen(nz)+0.5)/double(nz)

; Get filters
k_load_filters,filterpath+'/'+filterlist,filter_n,filter_lambda,filter_pass

; Set source object name
soname=filepath('libkcorrect.so', $
                root_dir=getenv('KCORRECT_DIR'), subdirectory='lib')

; Set rmatrix
nl=long(n_elements(lambda))-1l
nb=long(n_elements(bmatrix)/nl)
nk=long(n_elements(filterlist))
rmatrix=dblarr(nz,nb,nk)
retval=call_external(soname, 'idl_k_create_r', double(rmatrix), $
                     long(nk), long(nb), double(bmatrix), $
                     double(lambda), long(nl), double(zvals), long(nz), $
                     long(filter_n), double(filter_lambda), $
                     double(filter_pass), long(max(filter_n)), $
                     double(glambda), double(glambdaexp))

klog,'Done.'

end
;------------------------------------------------------------------------------
