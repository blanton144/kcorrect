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
;   coeffs        - [nv, ngals] coefficients 
;   redshift      - [ngals] redshift for each galaxy (where you want to
;                   calculate maggies in model) 
; KEYWORDS:
;   silent        - shut up
; OUTPUTS:
;   reconstruct_maggies - [nk, ngals] maggies of each
; OPTIONAL INPUT/OUTPUTS:
;   filterlist    - [nk] list of files with filter information 
;   filterpath    - path for filterlist (default $KCORRECT_DIR/data/filters)
;   vmatrix       - [nl,nk] templates spanning SED space 
;   lambda        - [nl+1] wavelengths for templates 
;   rmatrix       - [nz, nv, nk] look up table for bmatrix and filter 
;                   information 
;   zvals         - [nz] look up table for rmatrix 
;   zmin, zmax    - redshifts limits of lookup table (default 0., 2)
;   nz            - number of redshifts in lookup table (default 1000)
;   band_shift    - shift to apply to bandpasses (default 0.)
; COMMENTS:
;   Reconstruct AB galaxy maggies given an observed redshift and a
;   shift to apply to the bandpasses (band_shift) under the assumption
;   that the bolometric flux is conserved. To reconstruct the
;   observed galaxy maggies:
; 
;      k_reconstruct_maggies,coeffs,redshift,reconstruct_maggies
; 
;   To construct what would be observed if the galaxies were observed 
;   at z=0. through a bandpass blueshifted by z=0.1:
; 
;      k_reconstruct_maggies,coeffs,replicate(0.,ngals),reconstruct_maggies, $
;         band_shift=replicate(0.1,ngals)
; EXAMPLES:
;   Given coeffs from a call to k_fit_coeff or kcorrect, create 
;   reconstructed maggies for each galaxy using default templates
; 
;   IDL> k_reconstruct_maggies,coeffs,redshift,reconstruct_maggies
;  
;   To K-correct all the maggies to the same redshift (0.1):
;
;   IDL> k_reconstruct_maggies, coeffs,replicate(0.1,n_elements(redshift)), $
;        reconstruct_maggies
; REVISION HISTORY:
;   05-Jan-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_reconstruct_maggies,coeffs,redshift,reconstruct_maggies, $
                          band_shift=band_shift,zvals=zvals, $
                          rmatrix=rmatrix,vmatrix=vmatrix,lambda=lambda, $
                          filterpath=filterpath, filterlist=filterlist, $
                          zmin=zmin,zmax=zmax,nz=nz,silent=silent

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
      zmin=zmin,zmax=zmax, nz=nz, filterpath=filterpath, $
      band_shift=band_shift, silent=silent
endif else if(NOT keyword_set(rmatrix) OR NOT keyword_set(zvals)) then begin
    klog, 'need to specify rmatrix and zvals or vmatrix, lambda and filterlist'
    return
endif
nz=long(n_elements(zvals))
nk=long(n_elements(rmatrix))/(nv*nz)

; Set source object name
soname=filepath('libkcorrect.'+kcorrect_so_ext(), $
                root_dir=getenv('KCORRECT_DIR'), subdirectory='lib')

; Call coefficient software
reconstruct_maggies=fltarr(nk,ngalaxy)
retval=call_external(soname, 'idl_k_reconstruct_maggies', float(zvals), $
                     long(nz), float(rmatrix), long(nk), long(nv), $
                     float(coeffs), float(redshift), $
                     float(reconstruct_maggies), long(ngalaxy) )

end
;------------------------------------------------------------------------------
