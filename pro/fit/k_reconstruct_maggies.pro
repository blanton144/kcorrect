;+
; NAME:
;   k_reconstruct_maggies
; PURPOSE:
;   Reconstruct AB galaxy maggies given an observed redshift 
; CALLING SEQUENCE:
;   k_reconstruct_maggies,coeffs,redshift,reconstruct_maggies [, $
;       band_shift=, zvals=, rmatrix=, vmatrix=, lambda=, $
;       filterpath=, filterlist=, zmin=, zmax=, nz= ]
; INPUTS:
;   coeff    - coefficients [N_template, N_gal]
;   redshift      - redshift for each galaxy (where you want to
;                   calculate maggies in model) [N_gal]
; OPTIONAL INPUTS:
; OUTPUTS:
; OPTIONAL INPUT/OUTPUTS:
;   filterlist    - list of files with filter information [N_band]
;   vmatrix       - templates spanning SED space [N_lambda, N_dim]
;   lambda        - wavelengths for orthogonal templates [N_lambda+1]
;   rmatrix       - look up table for bmatrix and filter information 
;                   [N_z, N_dim, N_band]
;   zvals         - look up table for rmatrix [N_z]
;   band_shift    - shift to apply to bandpasses (default 0.)
; COMMENTS:
;   Outputs maggies in maggies. Does not calculate errors. 
; 
;   Reconstruct AB galaxy maggies given an observed redshift and a
;   shift to apply to the bandpasses (band_shift). To reconstruct the
;   observed galaxy maggies:
; 
;      k_reconstruct_maggies,coeff,redshift,reconstruct_maggies
; 
;   To construct what would be observed if the galaxy were observed 
;   at z=0. through a bandpass blueshifted by z=0.1:
; 
;      k_reconstruct_maggies,coeff,replicate(0.,ngals),reconstruct_maggies, $
;         band_shift=replicate(0.1,ngals)
; EXAMPLES:
;   Given coeffs from a call to k_fit_coeff, create reconstructed maggies
;   for each galaxy using default templates
; 
;   IDL> k_reconstruct_maggies,coeffs,redshift,reconstruct_maggies,/default
;  
;   To K-correct all the maggies to the same redshift (0.1):
;
;   IDL> k_reconstruct_maggies, coeffs,replicate(0.1,n_elements(redshift)), $
;        reconstruct_maggies, /default
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   05-Jan-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_reconstruct_maggies,coeffs,redshift,reconstruct_maggies, $
                          band_shift=band_shift,zvals=zvals, $
                          rmatrix=rmatrix,vmatrix=vmatrix,lambda=lambda, $
                          filterpath=filterpath, filterlist=filterlist, $
                          zmin=zmin,zmax=zmax,nz=nz

; Need at least 3 parameters
if (N_params() LT 3) then begin
    print, 'Syntax - k_reconstruct_maggies, coeffs, redshift, reconstruct_maggies $'
    print, '    [band_shift=, zvals=, rmatrix=, vmatrix=, lambda=, filterlist=, $'
    print, '      filterpath=]'
    return
endif

ngalaxy=long(n_elements(redshift))
nv=long(n_elements(coeffs))/ngalaxy
if (keyword_set(vmatrix) AND keyword_set(filterlist) AND keyword_set(lambda)) $
  then begin
    k_projection_table, rmatrix, vmatrix, lambda, zvals, filterlist, $
      zmin=zmin,zmax=zmax, nz=nz, filterpath=filterpath, band_shift=band_shift
endif else if(NOT keyword_set(rmatrix) OR NOT keyword_set(zvals)) then begin
    klog, 'need to specify rmatrix and zvals or vmatrix, lambda and filterlist'
    return
endif
nz=long(n_elements(zvals))
nk=long(n_elements(rmatrix))/(nv*nz)

; Set source object name
soname=filepath('libkcorrect.so', $
                root_dir=getenv('KCORRECT_DIR'), subdirectory='lib')

; Call coefficient software
reconstruct_maggies=fltarr(nk,ngalaxy)
retval=call_external(soname, 'idl_k_reconstruct_maggies', float(zvals), $
                     long(nz), float(rmatrix), long(nk), long(nv), $
                     float(coeffs), float(redshift), $
                     float(reconstruct_maggies), long(ngalaxy) )

end
;------------------------------------------------------------------------------