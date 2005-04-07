;+
; NAME:
;   k_lambda_eff
; PURPOSE:
;   Get effective wavelengths of filters in Angstroms
; CALLING SEQUENCE:
;   efflam=k_lambda_eff([filterlist=, filterpath=, band_shift=])
; OPTIONAL INPUTS:
;   filterlist - list of filter names (default
;                ['sdss_u0.par','sdss_g0.par', 'sdss_r0.par',
;                'sdss_i0.par', 'sdss_z0.par'])
;   filterpath - path for filter files; default to kcorrect repository
;   band_shift - blueshift to apply to bandpass; default to 0.0
; COMMENTS:
;   Calculates the effective wavelength using the Schneider et al 1983
;   defn (as quoted in Fukugita et al 1996). Returns results in
;   Angstroms.
; PROCEDURES CALLED:
;   k_load_filters
; REVISION HISTORY:
;   17-Jan-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function k_lambda_eff,filterlist=filterlist, $
                      filterpath=filterpath, $
                      band_shift=band_shift

if(NOT keyword_set(filterlist)) then $
  filterlist=['sdss_u0.par','sdss_g0.par','sdss_r0.par','sdss_i0.par', $
              'sdss_z0.par']
if(NOT keyword_set(band_shift)) then band_shift=0.

; load in the filters
k_load_filters,filterlist,filter_n,filter_lambda,filter_pass
filter_lambda=filter_lambda/(1.+band_shift)

lambda_eff=fltarr(n_elements(filterlist))
for k=0, n_elements(filterlist)-1 do begin
  dloglambda=fltarr(filter_n[k])
  dloglambda[0]=alog(filter_lambda[1,k])-alog(filter_lambda[0,k])
  dloglambda[filter_n[k]-1]=alog(filter_lambda[filter_n[k]-1,k])- $
    alog(filter_lambda[filter_n[k]-2,k])
  dloglambda[1:filter_n[k]-2]=0.5*(alog(filter_lambda[2:filter_n[k]-1,k])- $
                                   alog(filter_lambda[0:filter_n[k]-3,k]))
  denom=total(filter_pass[0:filter_n[k]-1,k]*dloglambda,/double)
  numer=total(filter_pass[0:filter_n[k]-1,k]*dloglambda* $
              alog(filter_lambda[0:filter_n[k]-1,k]),/double)
  lambda_eff[k]=exp(numer/denom)
endfor

return,lambda_eff

end
;------------------------------------------------------------------------------

