;+
; NAME:
;   lf_schechter_mden
; PURPOSE:
; INPUTS:
; OPTIONAL INPUTS:
; KEYWORDS:
; OUTPUTS:
; OPTIONAL OUTPUTS:
; BUGS:
; DEPENDENCIES:
; REVISION HISTORY:
;   2002-7-04  written - Blanton
;-
function lf_schechter_mden_plus,schechter_plus,limits=limits,bin=bin

if(NOT keyword_set(bin)) then bin=0.02
if(NOT keyword_set(limits)) then $
  limits=[schechter_plus.mstar-5.,schechter_plus.mstar+15.]

nbin=long((limits[1]-limits[0])/bin)
am=limits[0]+(limits[1]-limits[0])*(dindgen(nbin)+0.5)/double(nbin)
dm=am[1]-am[0]
lum=10^(-0.4*(am+20.))
sch=lf_schechter_plus(am,schechter_plus.phistar,schechter_plus.mstar, $
                      schechter_plus.alpha,schechter_plus.phiplus, $
                      schechter_plus.alphaplus)
lumden=total(dm*lum*sch)
mden=-2.5*alog10(lumden)-20.

return,mden

end
