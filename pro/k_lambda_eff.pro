;+
; NAME:
;   k_lambda_eff
; PURPOSE:
;   Get effective wavelengths of filters in Angstroms
; CALLING SEQUENCE:
;   efflam=k_lambda_eff([filterlist=, filterpath=, band_shift=])
; INPUTS:
; OPTIONAL INPUTS:
;   filterlist - list of filter names (default ['sdss_u0','sdss_g0','sdss_r0',
;                'sdss_i0','sdss_z0'])
;   filterpath - path for filter files (default
;                getenv('KCORRECT_DIR')+'/data/filters')
;   band_shift - blueshift to apply to bandpass
; OUTPUTS:
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
;   Calculates the effective wavelength using the Schneider et al 1983
;   defn (as quoted in Fukugita et al 1996). Returns results in
;   Angstroms.
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
;   k_load_filters
; REVISION HISTORY:
;   17-Jan-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function k_lambda_eff,filterlist=filterlist, $
                      filterpath=filterpath, $
                      band_shift=band_shift

if(NOT keyword_set(filterpath)) then $
  filterpath=getenv('KCORRECT_DIR')+'/data/filters'
if(NOT keyword_set(filterlist)) then $
  filterlist=['sdss_u0','sdss_g0','sdss_r0','sdss_i0','sdss_z0']
if(NOT keyword_set(band_shift)) then band_shift=0.

; load in the filters
fullfilterlist=filterpath+'/'+filterlist+'.dat'
k_load_filters,fullfilterlist,filter_n,filter_lambda,filter_pass
filter_lambda=filter_lambda/(1.+band_shift)

lambda_eff=dblarr(n_elements(fullfilterlist))
for k=0, n_elements(fullfilterlist)-1 do begin
  dloglambda=dblarr(filter_n[k])
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

