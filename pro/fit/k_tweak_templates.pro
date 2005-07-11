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
; COMMENTS:
;   Minimizes chi^2 of the difference between the maggies and the 
;   reconstructed maggies by multiplying the template spectra by 
;   a low-order polynomial and running with that.
; REVISION HISTORY:
;   24-Apr-2003  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function k_tt_tweak_vmatrix, lambda, vmatrix, tweakpars

vmatrix_tweaked=vmatrix
nl=n_elements(lambda)-1L
nv=n_elements(vmatrix)/nl
ntweak=n_elements(tweakpars)/nv
logl=(alog(0.5*(lambda[0:nl-1]+lambda[1:nl])))
tweakcenters=logl[0]+(logl[nl-1]-logl[0])*(dindgen(ntweak)-1.)/ $
  double(ntweak-3L)
for i=0L, nv-1L do begin
    tweakval=dblarr(nl)
    for j=0L, ntweak-1L do begin
        xx=(logl-tweakcenters[j])/(tweakcenters[1]-tweakcenters[0])
        tweakval=tweakval+tweakpars[j,i]*k_bspline2(xx)
    endfor
    vmatrix_tweaked[*,i]=vmatrix_tweaked[*,i]*exp(tweakval)
endfor
return,vmatrix_tweaked

end
; 
pro k_tt_params_interpret, params, maggies_factor, tweakpars

common k_tt_common, k_tt_maggies, k_tt_maggies_ivar, k_tt_coeffs, $
  k_tt_vmatrix, k_tt_lambda, k_tt_filterlist, k_tt_redshift, $
  k_tt_tweak_maggies, k_tt_smoothness

nk=n_elements(k_tt_filterlist)
nv=n_elements(k_tt_vmatrix)/(n_elements(k_tt_lambda)-1L)
maggies_factor=dblarr(nk)+1.D
nmaggies=0L
if(keyword_set(k_tt_tweak_maggies)) then begin
    maggies_factor[0:1]=exp(params[0:1])
    maggies_factor[3:7]=exp(params[2:6])
    maggies_factor[2]=exp(0.)
    nmaggies=7L
endif
ntweak=(n_elements(params)-(nmaggies))/nv
tweakpars=dblarr(ntweak,nv)
tweakpars[*,*]=params[nmaggies:n_elements(params)-1L]

end
;
function k_tt_func, params

common k_tt_common

nk=n_elements(k_tt_filterlist)
nv=n_elements(k_tt_vmatrix)/(n_elements(k_tt_lambda)-1L)
k_tt_params_interpret, params, maggies_factor, tweakpars

; multiply in tweak
vmatrix_tweaked= $
  k_tt_tweak_vmatrix(k_tt_lambda, k_tt_vmatrix, tweakpars)

k_reconstruct_maggies,k_tt_coeffs,k_tt_redshift,maggies_model, $
  lambda=k_tt_lambda,vmatrix=vmatrix_tweaked, filterlist=k_tt_filterlist, $
  nz=100L
  
; tweak zero points
maggies_obs=k_tt_maggies
maggies_obs_ivar=k_tt_maggies_ivar
for i=0, n_elements(k_tt_filterlist)-1L do begin
    maggies_obs[i,*]=k_tt_maggies[i,*]*maggies_factor[i]
    maggies_obs_ivar[i,*]=k_tt_maggies_ivar[i,*]/maggies_factor[i]^2
endfor

; add tweakpars constraints
ntweak=n_elements(tweakpars)/nv
devpars=dblarr(nv*ntweak)
for i=0, nv-1L do begin
    devpars[i*ntweak+0]=tweakpars[ntweak/2L,i]/k_tt_smoothness
    devpars[i*ntweak+lindgen(ntweak-1L)+1L]= $
      (tweakpars[1:ntweak-1L,i]-tweakpars[0:ntweak-2,i])/k_tt_smoothness
endfor

dev=(maggies_obs-maggies_model)*sqrt(maggies_obs_ivar)
dev=reform(dev,n_elements(dev))
dev=[dev,devpars]
help,dev
return, dev

end
;
pro k_tt_iterproc, myfunc, params, iter, fnorm, FUNCTARGS=fcnargs, $
                   PARINFO=parinfo, QUIET=quiet

common k_tt_common

nl=n_elements(k_tt_lambda)-1L
nk=n_elements(k_tt_filterlist)
nv=n_elements(k_tt_vmatrix)/(n_elements(k_tt_lambda)-1L)
k_tt_params_interpret, params, maggies_factor, tweakpars

!P.MULTI=[0,1,nv]
test=dblarr(nl,nv)+1.
test_tweaked= $
  k_tt_tweak_vmatrix(k_tt_lambda, test, tweakpars)
;window,0
for i=0, nv-1 do begin
    plot, k_tt_lambda,test_tweaked[*,i],/xlog
    oplot, k_tt_lambda,test[*,i]
endfor
;window,1
dev=k_tt_func(params)
;!P.MULTI=[0,1,nk]
;ng=n_elements(k_tt_redshift)
;for k=0L, nk-1L do begin
;    plot,k_tt_redshift,dev[nk*lindgen(ng)+k],psym=3,yra=[-3.,3.]
;endfor

splog,'***********************************************'
splog,'chi2 = '+string(fnorm)
splog,'***********************************************'

end
; 
pro k_tweak_templates, maggies, maggies_ivar, redshift, coeffs, vmatrix, $
                       lambda, filterlist=filterlist, $
                       maggies_factor=maggies_factor, tweakpars=tweakpars, $
                       vmatrix_tweaked=vmatrix_tweaked

common k_tt_common

if(NOT keyword_set(filterlist)) then $
  filterlist=['sdss_u0','sdss_g0','sdss_r0','sdss_i0','sdss_z0']+'.par'

k_tt_maggies=maggies
k_tt_maggies_ivar=maggies_ivar
k_tt_redshift=redshift
k_tt_coeffs=coeffs
k_tt_vmatrix=vmatrix
k_tt_lambda=lambda
k_tt_filterlist=filterlist
k_tt_tweak_maggies=0
if(keyword_set(k_tt_tweak_maggies)) then $
  nmaggies=n_elements(filterlist)-1L $
else $
  nmaggies=0L
ntweak=80
k_tt_smoothness=1.00

nv=n_elements(k_tt_vmatrix)/(n_elements(k_tt_lambda)-1L)
start=0.0*randomn(seed,nmaggies+ntweak*nv,/double)
pi=replicate({fixed:0, limited:[0,0], limits:[0.D,0.D], step:0.01}, $
             n_elements(start))
params=mpfit('k_tt_func',start,iterproc='k_tt_iterproc',ftol=1.d-16, $
             parinfo=pi)

k_tt_params_interpret, params, maggies_factor, tweakpars
vmatrix_tweaked= $
  k_tt_tweak_vmatrix(lambda, vmatrix, tweakpars)

end
;------------------------------------------------------------------------------
