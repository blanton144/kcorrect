;+
; NAME:
;   k_fit_nonneg_photoz
; PURPOSE:
;   Fit redshifts and coefficients to eigentemplates given photometry
; CALLING SEQUENCE:
;   k_fit_nonneg_photoz, maggies, maggies_err, vmatrix, lambda,
;     redshift=redshift, filterlist=filterlist, filterpath=filterpath, 
;     rmatrix=rmatrix, zvals=zvals, chi2=chi2
; INPUTS:
; OPTIONAL INPUTS:
; OUTPUTS:
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   13-Feb-2003  by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_fit_nonneg_photoz, maggies, maggies_err, vmatrix, $
                         lambda,redshift=redshift, $
                         band_shift=band_shift, version=version, vpath=vpath, $
                         filterlist=filterlist, filterpath=filterpath, $
                         maxiter=maxiter, rmatrix=rmatrix,zvals=zvals, $
                         zmin=zmin,zmax=zmax, chi2=chi2

if(n_elements(zmin) eq 0) then zmin=0.001
if(n_elements(zmax) eq 0) then zmax=1.
if(n_elements(nz) eq 0) then nz=300

k_create_r,rmatrix,vmatrix,lambda,zvals,filterlist,filterpath=filterpath, $
  zmin=zmin, zmax=zmax, nz=nz

nbands=n_elements(filterlist)
ngals=n_elements(maggies)/nbands
redshift=dblarr(ngals)
chi2=dblarr(ngals)
chi2vals=dblarr(nz)
use_maggies=reform(maggies,nbands,ngals)
use_maggies_err=reform(maggies_err,nbands,ngals)
for i=0, ngals-1L do begin
    splog,'i='+string(i)
    for j=0L, nz-1L do begin
        ztest=zvals[j]
        coeffs=k_fit_nonneg(use_maggies[*,i], use_maggies_err[*,i], $
                            redshift=ztest,filterlist=filterlist, $
                            filterpath=filterpath,maxiter=maxiter, $
                            rmatrix=rmatrix,zvals=zvals, $
                            chi2=chi2test,/quiet)
        chi2vals[j]=chi2test
    endfor
    chi2[i]=min(chi2vals,indx)
    redshift[i]=zvals[indx]
    splog,'chi2= '+string(chi2[i])
    splog,'redshift= '+string(redshift[i])
endfor

end
;------------------------------------------------------------------------------
