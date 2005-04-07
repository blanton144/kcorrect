;+
; NAME:
;   lf_fit_schechter
; PURPOSE:
;   fit schechter function to set of points
; USAGE:
;   lf_fit_schechter, absmag, phi, phierr, schechter [, mden=]
; INPUTS:
;   absmag - absolute magnitude at center of each bin
;   phi - luminosity function at each bin
;   phierr - error in l.f. at each bin
; OUTPUTS:
;   schechter - structure with
;                  .PHISTAR 
;                  .MSTAR 
;                  .ALPHA 
;                  .PHISTAR_ERR
;                  .MSTAR_ERR
;                  .ALPHA_ERR
; OPTIONAL OUTPUTS:
;   mden - total luminosity density in magnitudes of fit
; REVISION HISTORY:
;   2002-7-04  written - Blanton
;-
function lf_fit_schechter_func,params,x=x,y=y,err=err

phimodel=lf_schechter(x,params[0],params[1],params[2])
dev=((y-phimodel)/err)
return,dev

end
;
pro lf_fit_schechter,am,phi,phistddev,schechter,mden=mden

if(n_tags(schechter) eq 0) then $
  schechter={phistar:1.d-2, mstar:-20.d, alpha:-1.d, $
             phistar_err:0.d, mstar_err:0.d, alpha_err:0.d}
start=dblarr(3)
start[0]=schechter.phistar
start[1]=schechter.mstar
start[2]=schechter.alpha

functargs={x:am,y:phi,err:phistddev}
params=mpfit('lf_fit_schechter_func',start,functargs=functargs,perror=perror)

schechter.phistar=params[0]
schechter.mstar=params[1]
schechter.alpha=params[2]
schechter.phistar_err=perror[0]
schechter.mstar_err=perror[1]
schechter.alpha_err=perror[2]
mden=lf_schechter_mden(schechter)

help,/st,schechter

end
