;+
; NAME:
;   k_reconstruct_nn
;
; PURPOSE:
;   Given a set of coefficients from the fit in k_fit_nn,
;   reconstruct the magnitudes in any band 
;
; CALLING SEQUENCE:
;
; INPUTS:
;   coeff    - coefficients [N_template, N_gal]
;
; OPTIONAL INPUTS:
;
; KEYWORDS:
; OUTPUTS:
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
;
; BUGS:
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
pro k_reconstructnn_maggies, coeff, galaxy_z, reconstruct_maggies, $
                      band_shift=band_shift, ematrix=ematrix, bmatrix=bmatrix, lambda=lambda, $
                      zvals=zvals, filterlist=filterlist, rmatrix=rmatrix, $
                      version=version, vpath=vpath, filterpath=filterpath

; Need at least 6 parameters
if (N_params() LT 3) then begin
    print, 'Syntax - k_reconstruct_nn, coeff, galaxy_z, reconstruct_maggies $'  
    print, '        [, band_shift=, bmatrix=, lambda=lambda, zvals=, filterlist=, ' 
    print, '         rmatrix=, version=, vpath=, filterpath=]'  
    return
endif

if(NOT keyword_set(filterpath)) then $j
  filterpath=getenv('KCORRECT_DIR')+'/data/filters'

; Get bmatrix and stuff from files if necessary
if(keyword_set(version)) then begin 
    if(NOT keyword_set(vpath)) then $
      vpath=getenv('KCORRECT_DIR')+'/data/etemplates'
    k_load_ascii_table,bmatrix,vpath+'/bmatrix.'+version+'.dat'
    k_load_ascii_table,ematrix,vpath+'/ematrix.'+version+'.dat'
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
nt=long(n_elements(ematrix)/(nb))
ngalaxy=long(n_elements(coeff)/nt)

; now reconstruct for each galaxy
reconstruct_maggies=dblarr(nk,ngalaxy)
for i=0L, ngalaxy-1L do begin
; construct rjk(z)
    rlocal=dblarr(nk,nb)
    for k=0L, nk-1L do $
      for j=0L, nb-1L do $
      rlocal[k,j]=interpolate(rmatrix[*,j,k],double(nz)* $
                              (galaxy_z[i]-zvals[0])/(zvals[nz-1]-zvals[0]))
  				
; now evaluate it
  reconstruct_maggies[*,i]=(ematrix##rlocal)#coeff[*,i]
endfor

end
;------------------------------------------------------------------------------
