;+
; NAME:
;   k_tweak_templates
; PURPOSE:
;   Given a set of templates, tweak em to better fit
; CALLING SEQUENCE:
;   k_tweak_templates, maggies, maggies_ivar, coeffs, vmatrix,
;      lambda, filterlist=filterlist
; INPUTS:
;   maggies - maggies for each galaxy
;   maggies_ivar - inverse variance of maggie values
;   coeffs  - coefficients 
;   vmatrix - array of template spectra
;   lambda  - lambda pixel boundaries
; OPTIONAL INPUTS:
;   filterlist - list of filters to use
; OUTPUTS:
;   
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
;   Minimizes chi^2 of the difference between the maggies and the 
;   reconstructed maggies by multiplying the template spectra by 
;   a low-order polynomial and running with that.
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   24-Apr-2003  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function k_tt_func, params

common k_tt_com, k_tt_maggies, k_tt_maggies_ivar, k_tt_coeffs, $
  k_tt_vmatrix, k_tt_lambda, k_tt_filterlist, k_tt_redshift

nk=n_elements(k_tt_filterlist)
nt=n_elements(k_tt_vmatrix)/(n_elements(k_tt_lambda)-1L)
maggies_factor=exp(params[0:nk-1])
alpha=params[nk:nk+nt-1]
beta=params[nk+nt:nk+nt+nt-1]
gamma=params[nk+2*nt:nk+2*nt+nt-1]

nv=n_elements(k_tt_vmatrix)/(n_elements(k_tt_lambda)-1L)

; multiply in tweak
vmatrix_tweaked=k_tt_vmatrix
logl=alog(0.5*(lambda[0:n_elements(lambda)-2]+lambda[1:n_elements(lambda)-1]))_
logl0=logl[n_elements(ll)/2]
difflogl=logl-logl0
for i=0L, nt-1L do begin
    vmatrix_tweaked[*,i]=vmatrix_tweaked* $
      exp(difflogl*alpha[i]+ $
          difflogl^2*beta[i]+ $
          difflogl^3*gamma[i])
endfor

k_reconstruct_maggies,k_tt_coeffs,k_tt_redshift,maggies_model, $
  lambda=k_tt_lambda,bmatrix=vmatrix_tweaked,ematrix=identity(nv), $
  filterlist=k_tt_filterlist

; tweak zero points
maggies_obs=k_tt_maggies
maggies_obs_ivar=k_tt_maggies
for i=0, n_elements(k_tt_filterlist)-1L do begin
    maggies_obs[i,*]=k_tt_maggies[i,*]*maggies_factor[i]
    maggies_obs_ivar[i,*]=k_tt_maggies_ivar[i,*]/maggies_factor[i]^2
endfor

dev=(maggies_obs-maggies_model)*sqrt(maggies_obs_ivar)
dev=reform(dev,n_elements(dev))
help,dev
help,maggies_obs
help,maggies_obs_ivar
help,maggies_model
return, dev

end
; 
pro k_tweak_templates, maggies, maggies_ivar, redshift, coeffs, vmatrix, $
                       lambda, filterlist=filterlist, $
                       maggies_factor=maggies_factor

common k_tt_com

if(NOT keyword_set(filterlist)) then $
  filterlist=['sdss_u0','sdss_g0','sdss_r0','sdss_i0','sdss_z0']+'.dat'

k_tt_maggies=maggies
k_tt_maggies_ivar=maggies_ivar
k_tt_redshift=redshift
k_tt_coeffs=coeffs
k_tt_vmatrix=vmatrix
k_tt_lambda=lambda
k_tt_filterlist=filterlist

nt=n_elements(k_tt_vmatrix)/(n_elements(k_tt_lambda)-1L)
start=replicate(0.,n_elements(filterlist)+3*nt)
params=mpfit('k_tt_func',start)

maggies_factor=exp(params[0:n_elements(k_tt_filterlist)-1L])

end
;------------------------------------------------------------------------------
