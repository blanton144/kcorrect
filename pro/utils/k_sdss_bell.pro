;+
; NAME:
;   k_sdss_bell
; PURPOSE:
;   Bell and de Jong stellar masses given SDSS absolute magnitudes
; CALLING SEQUENCE:
;   smass= k_sdss_bell(absmag)
; INPUTS:
;   absmag - [5, N] SDSS ugriz absolute magnitudes
; OUTPUTS:
;   smass - stellar mass in solar masses
; REVISION HISTORY:
;   2005-04-07 MRB, NYU
;-
;------------------------------------------------------------------------------
function k_sdss_bell, absmag

gmr=absmag[1,*]-absmag[2,*]

ii=where(abs(gmr) gt 1.e+10 OR gmr ne gmr, nii)
if(nii gt 0) then gmr[ii]=0.75

absmv=absmag[2,*]+0.2876+0.2246*(gmr-0.7778)-0.02 ;; vega
bmv=1.24*(gmr)-0.09  ;; vega
mtolv=10.^(-0.734+1.404*bmv) ;; from b&dj

vsolar=(k_solar_magnitudes(filterlist='bessell_V.par'))[0]
lvsolar=10.^(-0.4*(absmv-vsolar))
msolar=mtolv*lvsolar

return, msolar

end
