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
function k_tt_tweak_vmatrix, lambda, vmatrix, alpha, beta, gamma, delta, $
                             epsilon
vmatrix_tweaked=vmatrix
nv=n_elements(vmatrix)/(n_elements(lambda)-1L)
logl=alog(0.5*(lambda[0:n_elements(lambda)-2]+ $
               lambda[1:n_elements(lambda)-1]))
difflogl=((logl)-alog(7000.))
for i=0L, nv-1L do begin
    vmatrix_tweaked[*,i]=vmatrix_tweaked[*,i]* $
      exp(difflogl*alpha+ $
          difflogl^2*beta+ $
          difflogl^3*gamma+ $
          difflogl^4*delta+ $
          difflogl^5*epsilon)
endfor
return,vmatrix_tweaked
end
; 
pro k_tt_params_interpret, params, maggies_factor, alpha, beta, gamma, $
                           delta,epsilon

common k_tt_com, k_tt_maggies, k_tt_maggies_ivar, k_tt_coeffs, $
  k_tt_vmatrix, k_tt_lambda, k_tt_filterlist, k_tt_redshift

nk=n_elements(k_tt_filterlist)
nv=n_elements(k_tt_vmatrix)/(n_elements(k_tt_lambda)-1L)
maggies_factor=dblarr(nk)
maggies_factor[0:1]=exp(params[0:1])
maggies_factor[3:4]=exp(params[2:3])
maggies_factor[2]=exp(0.)
alpha=params[nk-1+0]
beta=params[nk-1+1]
gamma=params[nk-1+2]
delta=params[nk-1+3]
epsilon=params[nk-1+4]
end
;
function k_tt_func, params

common k_tt_com

nk=n_elements(k_tt_filterlist)
nv=n_elements(k_tt_vmatrix)/(n_elements(k_tt_lambda)-1L)
k_tt_params_interpret, params, maggies_factor, alpha, beta, gamma, $
  delta, epsilon

; multiply in tweak
vmatrix_tweaked= $
  k_tt_tweak_vmatrix(k_tt_lambda, k_tt_vmatrix, alpha, beta, gamma, $
                     delta, epsilon)

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
                       maggies_factor=maggies_factor, alpha=alpha, $
                       beta=beta, gamma=gamma, delta=delta, epsilon=epsilon, $
                       vmatrix_tweaked=vmatrix_tweaked

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

nv=n_elements(k_tt_vmatrix)/(n_elements(k_tt_lambda)-1L)
start=replicate(0.,n_elements(filterlist)-1L+5)
params=mpfit('k_tt_func',start)

k_tt_params_interpret, params, maggies_factor, alpha, beta, gamma, delta, $
  epsilon
vmatrix_tweaked= $
  k_tt_tweak_vmatrix(lambda, vmatrix, alpha, beta, gamma, delta, epsilon)

end
;------------------------------------------------------------------------------
