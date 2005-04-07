;+
; NAME:
;   lf_vals
; PURPOSE:
;   return values of the luminosity function
; INPUTS:
;   nphi             number of "bins"
;   sample_absmmin   minimum abs mag for the sample
;   sample_absmmax   maximum abs mag for the sample
;   sigabsmag        width of gaussian
;   phi              amplitude of bin 
; OPTIONAL INPUTS:
;   subsample        how much to subsample for plot
; KEYWORDS:
; OUTPUTS:
; OPTIONAL OUTPUTS:
; BUGS:
; DEPENDENCIES:
; REVISION HISTORY:
;   2002-5-22  written - Blanton
;-
pro lf_vals,absmk,phi,sigabsmag,sample_absmmin,sample_absmmax, $
            amvals,phivals,subsample=subsample,factor=factor, $
            invsig=invsig

; defaults
if(NOT keyword_set(subsample)) then subsample=10L

; settings
pi=3.14159265358979D+0

; set amvals
nphi=n_elements(phi)
amvals=sample_absmmin+(dindgen(subsample*nphi)+0.5)/ $
  double(subsample*nphi)*(sample_absmmax-sample_absmmin)

; set phivals
phivals=dblarr(subsample*nphi)
factor=1.d/(sqrt(2*pi)*sigabsmag)
invsig=0.5d/(sigabsmag^2)
for i=0L, nphi-1L do begin
    phivals=phivals+phi[i]*factor*exp(-invsig*(absmk[i]-amvals)^2)
endfor

end
