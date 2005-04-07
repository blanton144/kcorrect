;+
; NAME:
;   lf_schechter_nden
; PURPOSE:
;   calculate total number density from schechter function
; USAGE:
;   mden= lf_schechter_nden(schechter [, limits=, bin=])
; INPUTS:
;   schechter - structure containing schechter pars 
;               (.PHISTAR, .MSTAR, .ALPHA)
; OPTIONAL INPUTS:
;   limits - [2] limits of integration (default to [MSTAR-5., MSTAR+10.])
;   bin - binning of integration (default to 0.02)
; OUTPUTS:
;   nden - total number density 
; REVISION HISTORY:
;   2002-7-04  written - Blanton
;-
function lf_schechter_nden,schechter,limits=limits,bin=bin

if(NOT keyword_set(bin)) then bin=0.02
if(NOT keyword_set(limits)) then $
  limits=[schechter.mstar-5.,schechter.mstar+10.]

nbin=long((limits[1]-limits[0])/bin)
am=limits[0]+(limits[1]-limits[0])*(dindgen(nbin)+0.5)/double(nbin)
dm=am[1]-am[0]
lum=10^(-0.4*(am+20.))
sch=lf_schechter(am,schechter.phistar,schechter.mstar,schechter.alpha)
nden=total(dm*sch)

return,nden

end
