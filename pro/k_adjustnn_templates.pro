;+
; NAME:
;   k_fitnn_coeff
;
; PURPOSE:
;   Given a set of maggies, redshifts, and spectra to fit to, find the 
;   nonnegative combination of the inputs which produces the outputs.
;
; CALLING SEQUENCE:
;   k_fit_nn, galaxy_maggies, galaxy_invvar, galaxy_z, coeff, $
;      [ematrix=, zvals=, filterlist=, bmatrix=, lambda=, rmatrix=]
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
; KEYWORDS:
;   covar - galaxy_invvar is actually a set of covariance matrices
;
; OUTPUTS:
;   coeff    - output coefficients [N_template, N_gal]
;
; OPTIONAL INPUT/OUTPUTS:
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
;   IDL> k_fit_nn,galaxy_maggies,galaxy_invvar,galaxy_z, coeff
; 
;   Then you can pass "coeff" into k_reconstruct_nn
;
; BUGS:
;   Will fail if N_dim is unity.
;
; PROCEDURES CALLED:
;   k_load_ascii_table
;   k_create_r
;   Dynamic link to idl_k_fit_nn.c in libkcorrect.so
;
; REVISION HISTORY:
;   04-Jancb2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_adjustnn_templates, galaxy_maggies, galaxy_invvar, galaxy_z, coeff, $
                          ematrix=ematrix, bmatrix=bmatrix, lambda=lambda, $
                          zvals=zvals,filterlist=filterlist,rmatrix=rmatrix, $
                          covar=covar, qarun=qarun, qastop=qastop

; Need at least 6 parameters
if (N_params() LT 4) then begin
    print, 'Syntax - k_fitnn_templates, galaxy_maggies, galaxy_invvar, galaxy_z, coeff, $'  
    print, '        [ematrix=, bmatrix=, lambda=, zvals=, filterlist=, rmatrix=, $'
    print, '         version=, vpath=, filterpath=, /covar]'
    return
endif

if(NOT keyword_set(filterpath)) then $j
  filterpath=getenv('KCORRECT_DIR')+'/data/filters'

; Set source object name
soname=filepath('libkcorrect.so', root_dir=getenv('KCORRECT_DIR'), $
                subdirectory='lib')

; Determine dimensions
ngalaxy=long(n_elements(galaxy_z))
nk=long(n_elements(galaxy_maggies)/ngalaxy)
if (keyword_set(bmatrix) AND keyword_set(filterlist) $
    AND keyword_set(lambda) AND keyword_set(ematrix)) $
  then begin
    k_create_r,rmatrix,bmatrix,lambda,zvals,filterlist,filterpath=filterpath
endif else begin
    if (NOT keyword_set(rmatrix)) then begin
        klog, 'need to specify rmatrix or bmatrix and filterlist'
    endif
endelse 
nz=long(n_elements(zvals))
nb=long(n_elements(rmatrix)/(nz*nk)) 	 
nt=long(n_elements(ematrix))/nb

; set inverse covariance matrices
galaxy_invcovar=dblarr(nk,nk,ngalaxy)
if(NOT keyword_set(covar)) then begin 
    for k=0, nk-1 do $
      galaxy_invcovar[k,k,*]=galaxy_invvar[k,*]
endif else begin
    for i=0L, ngalaxy-1L do $
      galaxy_invcovar[*,*,i]=invert(galaxy_invvar[*,*,i])
endelse

; construct inner matrices
innerc=dblarr(nb,nb,ngalaxy)
innerd=dblarr(nb,ngalaxy)
for i=0L, ngalaxy-1L do begin
    rlocal=dblarr(nk,nb)
    for k=0L, nk-1L do $
      for l=0L, nb-1L do $
      rlocal[k,l]=interpolate(rmatrix[*,l,k],double(nz)* $
                              (galaxy_z[i]-zvals[0])/(zvals[nz-1]-zvals[0]))
    innerc[*,*,i]=rlocal##galaxy_invcovar[*,*,i]##(transpose(rlocal))
		innerd[*,i]=rlocal##galaxy_invcovar[*,*,i]##galaxy_maggies[*,i]
endfor

; construct chi^2 matrices: x.C.x + d 
cmatrix=dblarr(nt*nb,nt*nb)
dmatrix=dblarr(nt*nb)
for j=0, nt-1 do begin
    for l=0, nb-1 do begin
        for jp=0, nt-1 do begin
            for lp=0, nb-1 do begin
                cmatrix[j*nb+l,jp*nb+lp]=total(coeff[j,*]*coeff[jp,*]* $
																							 innerc[l,lp,*],/double)
            endfor
        endfor
				dmatrix[j*nb+l]=-total(coeff[j,*]*innerd[l,*],/double)
    endfor
endfor

; now get derivative of chi2
dchi2=reform(ematrix,n_elements(ematrix))##cmatrix+dmatrix

; now perturb the templates, but only a little
dchi2_mag=sqrt(total(dchi2^2,/double))
ematrix_mag=sqrt(total(ematrix^2,/double))
scale=0.01*ematrix_mag/dchi2_mag
ematrix=ematrix+scale*dchi2
    
end
;------------------------------------------------------------------------------
