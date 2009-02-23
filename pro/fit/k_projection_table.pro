;+
; NAME:
;   k_projection_table
; PURPOSE:
;   Create lookup table for calculating maggies corresponding to spectra
; CALLING SEQUENCE:
;   k_projection_table, rmatrix, vmatrix, lambda, zvals, filterlist, $
;      [zmax=, zmin=, nz=, filterpath=, band_shift= ]
; INPUTS:
;   vmatrix       - templates spanning SED space [N_lambda, N_dim]
;   lambda        - wavelengths for orthogonal templates [N_lambda]
;   filterlist    - list of files with filter information [N_band]
; OPTIONAL INPUTS:
;   filterpath    - path for filters (default '$KCORRECT_DIR/data/filters')
;   band_shift    - shift the bands blueward by factor (1+band_shift)
; KEYWORDS:
;   verbose       - be loud
;   silent        - be silent (overrides verbose)
; OUTPUTS:
;   rmatrix       - look up table for vmatrix and filter information 
;                   [N_z, N_dim, N_band]
; OPTIONAL INPUT/OUTPUTS:
;   zvals         - look up table for rmatrix [N_z]
;   zmin, zmax, nz  - settings for setting zvals (default 0., 2., 1000)
; COMMENTS:
;   This tabulates the projection of each basis element v onto each
;   filter k, as a function of redshift. You only have to perform this
;   projection once, since every spectrum we will deal with will be a
;   linear combination of the basis elements b. To get a particular
;   redshift, you interpolate the rmatrix in the z direction.  Keep in
;   mind that this only creates the r matrix for a specific redshift
;   range; e.g. you have to change the defaults in order to consider
;   redshifts greater than unity.
; PROCEDURES CALLED:
;   k_load_filters
;   Dynamic link to idl_k_projection_table.c in libkcorrect.so
; REVISION HISTORY:
;   05-Jan-2002  Translated to IDL by Mike Blanton, NYU
;   22-Feb-2009  ZVALS now an optional input
;-
;------------------------------------------------------------------------------
pro k_projection_table, rmatrix, vmatrix, lambda, zvals, filterlist, $
                        zmin=zmin,zmax=zmax, nz=nz, filterpath=filterpath, $
                        band_shift=band_shift, verbose=verbose, silent=silent

if(n_params() ne 5) then begin
    print,'Syntax - k_projection_table, rmatrix, vmatrix, lambda, zvals, filterlist '
    print,'         [, zmin=, zmax=, nz=, filterpath=, band_shift=]'
    return
endif

if(NOT keyword_set(band_shift)) then band_shift=0.

if(keyword_set(silent)) then verbose=0
if(keyword_set(verbose)) then klog,'Creating rmatrix ...'

; Set defaults
if (n_elements(zmin) eq 0) then zmin=0.0
if (n_elements(zmax) eq 0) then zmax=2.0
if (NOT keyword_set(nz)) then nz=1000l

; Set zvals; the "not keyword_set()" business is for backwards
; compatibility with, e.g., deep_kcorrect
if (n_elements(zvals) eq 0L) or (keyword_set(zvals) eq 0) then $
  zvals=zmin+(zmax-zmin)*(findgen(nz)+0.5)/float(nz) $
else $
  nz = n_elements(zvals)

; Get filters
k_load_filters,filterlist,filter_n,filter_lambda,filter_pass, $
  filterpath=filterpath

; Set source object name
soname=filepath('libkcorrect.'+kcorrect_so_ext(), $
                root_dir=getenv('KCORRECT_DIR'), subdirectory='lib')

; Set rmatrix
nl=long(n_elements(lambda))-1l
nv=long(n_elements(vmatrix)/nl)
nk=long(n_elements(filterlist))
rmatrix=fltarr(nz,nv,nk)
retval=call_external(soname, 'idl_k_projection_table', float(rmatrix), $
                     long(nk), long(nv), float(vmatrix), $
                     float(lambda), long(nl), float(zvals), long(nz), $
                     long(filter_n), float(filter_lambda), $
                     float(filter_pass), float(band_shift), $
                     long(max(filter_n)))

if(keyword_set(verbose)) then klog,'Done.'

end
;------------------------------------------------------------------------------
