;+
; NAME:
;   lf_schechter_nden_plus
; PURPOSE:
;   calculate total density from schechter function
; USAGE:
;   nden= lf_schechter_nden_plus(schechter_plus [, limits=, bin=])
; INPUTS:
;   schechter_plus - structure containing double schechter pars 
;                    (.PHISTAR, .MSTAR, .ALPHA, .PHISTAR_PLUS, .ALPHA_PLUS)
; OPTIONAL INPUTS:
;   limits - [2] limits of integration (default to [MSTAR-5., MSTAR+10.])
;   bin - binning of integration (default to 0.02)
; OUTPUTS:
;   nden - total density 
; COMMENTS:
;   The double Schechter function is the sum of two Schechter
;   functions with the same MSTAR but different PHISTAR and ALPHA
;   values.
; REVISION HISTORY:
;   2002-7-04  written - Blanton
;-
function lf_schechter_nden_plus,schechter_plus,limits=limits,bin=bin

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
nden=total(dm*sch)

return,nden

end
