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
pro k_create_r, rmatrix, bmatrix, lambda, zvals, filterlist, zmin=zmin, zmax=zmax, nz=nz, filterpath=filterpath

klog,'Creating rmatrix ...'

; Set defaults
if (n_elements(zmin) eq 0) then zmin=1.0e-4
if (n_elements(zmax) eq 0) then zmax=1.0e-0
if (NOT keyword_set(nz)) then nz=1000l
if(NOT keyword_set(filterpath)) then $
  filterpath=getenv('KCORRECT_DIR')+'/data/filters'

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
                     double(filter_pass), long(max(filter_n)))

klog,'Done.'

end
;------------------------------------------------------------------------------
