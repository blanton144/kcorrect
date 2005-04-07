;+
; NAME:
;   lf_schechter_mden_plus
; PURPOSE:
;   calculate total luminosity density from schechter function
; USAGE:
;   mden= lf_schechter_mden_plus(schechter_plus [, limits=, bin=])
; INPUTS:
;   schechter_plus - structure containing double schechter pars 
;                    (.PHISTAR, .MSTAR, .ALPHA, .PHISTAR_PLUS, .ALPHA_PLUS)
; OPTIONAL INPUTS:
;   limits - [2] limits of integration (default to [MSTAR-5., MSTAR+10.])
;   bin - binning of integration (default to 0.02)
; OUTPUTS:
;   mden - total luminosity density in magnitude units
; COMMENTS:
;   The double Schechter function is the sum of two Schechter
;   functions with the same MSTAR but different PHISTAR and ALPHA
;   values.
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
