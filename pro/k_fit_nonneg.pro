;+
; NAME:
;   k_fit_nonneg
; PURPOSE:
;   Given a set of spectra, find the best combination of 
;   them which explain a given set of maggies and uncertainties,
;   with nonegative coefficients
; CALLING SEQUENCE:
;   coeffs = k_fit_nonneg(maggies, maggies_err, vmatrix, $
;   lambda
; INPUTS:
;   maggies - maggies to fit to
;   maggies_err - uncertainties in maggies
;   vmatrix - model spectra to use
;   lambda - model wavelength scale to use
; OPTIONAL INPUTS:
;   filterlist - bands in the kcorrect product to use
;   redshift - observed at higher redshift
;   band_shift - observed in shifted bandpasses
;   maxiter - maximum number of iterations in fit
; OUTPUTS:
; OPTIONAL OUTPUTS:
;   rmatrix - conversion matrix from spectra to maggies
;   zvals - redshift scale for rmatrix
;   chi2 - chi^2 of fit
; COMMENTS:
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   18-Jan-2003  Written MRB (NYU) -- based on lumden/pro/fit_sfh_maggies.pro
;-
;------------------------------------------------------------------------------
function k_fit_nonneg, maggies, maggies_err, vmatrix, $
                       lambda,redshift=redshift, $
                       band_shift=band_shift, version=version, vpath=vpath, $
                       filterlist=filterlist, filterpath=filterpath, $
                       maxiter=maxiter, rmatrix=rmatrix,zvals=zvals,chi2=chi2

if(n_elements(band_shift) eq 0) then band_shift=0.
if(n_elements(redshift) eq 0) then redshift=0.
if(NOT keyword_set(filterpath)) then $j
  filterpath=getenv('KCORRECT_DIR')+'/data/filters'

; Get vmatrix and stuff from files if necessary
if(keyword_set(version)) then begin 
    if(NOT keyword_set(vpath)) then $
      vpath=getenv('KCORRECT_DIR')+'/data/etemplates'
    k_load_ascii_table,vmatrix,vpath+'/vmatrix.'+version+'.dat'
    k_load_ascii_table,lambda,vpath+'/lambda.'+version+'.dat'
    if(NOT keyword_set(filterlist)) then begin
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
endif

; Determine dimensions and create conversion rmatrix if necessary
ngalaxy=long(n_elements(redshift))
nk=long(n_elements(maggies)/ngalaxy)
if (keyword_set(vmatrix) AND keyword_set(filterlist) AND keyword_set(lambda)) $
  then begin
    k_create_r,rmatrix,vmatrix,lambda,zvals,filterlist,filterpath=filterpath
endif else begin
    if (NOT keyword_set(rmatrix)) then begin
        klog, 'need to specify rmatrix or bmatrix and filterlist'
    endif
endelse 
nz=long(n_elements(zvals))
nv=long(n_elements(rmatrix)/(nz*nk)) 	 

; Do fit (this will later be integrated into the C-code, which will be
; faster and more portable)
coeffs=dblarr(nv,ngalaxy)
for i=0L, ngalaxy-1L do begin
    splog,i
    invcovar=dblarr(nv,nv)
    bb=dblarr(nv)
    for j=0L, nk-1L do begin
        tmp_rmatrix=dblarr(nv)
        tmp_interppos=double(n_elements(zvals))*(redshift[i]-zvals[0])/ $
          (zvals[n_elements(zvals)-1L]-zvals[0])
        for v=0L, nv-1L do $
          tmp_rmatrix[v]=interpolate(rmatrix[*,v,j], tmp_interppos)
        invcovar[*,*]=invcovar[*,*]+ $
          tmp_rmatrix#transpose(tmp_rmatrix)/maggies_err[j,i]^2
        bb[*]=bb[*]-tmp_rmatrix*maggies[j,i]/maggies_err[j,i]^2
    endfor
    offset=0.5*total((maggies[*,i]/maggies_err[*,i])^2,/double)

    start=replicate(1.e-3,nv)+1.e-3*randomu(seed,nv)
    coeffs[*,i]=nonneg_mult_update_solve(start,invcovar,bb,value=chi2, $
                                         offset=offset,/matrix,tol=tol, $
                                         niter=niter,maxiter=maxiter, $
                                         chi2tol=0.5,skip=100)
    splog,chi2
endfor

return, coeffs

end
