;+
; NAME:
;   lf_fit_schechter_plus
; PURPOSE:
;   fit "double" schechter function to set of points
; USAGE:
;   lf_fit_schechter_plus, absmag, phi, phierr, schechter_plus [, mden=]
; INPUTS:
;   absmag - absolute magnitude at center of each bin
;   phi - luminosity function at each bin
;   phierr - error in l.f. at each bin
; OUTPUTS:
;   schechter - structure with
;                  .MSTAR 
;                  .PHISTAR 
;                  .ALPHA 
;                  .PHIPLUS 
;                  .ALPHAPLUS 
;                  .MSTAR_ERR
;                  .PHISTAR_ERR
;                  .ALPHA_ERR
;                  .PHIPLUS_ERR
;                  .ALPHAPLUS_ERR
; OPTIONAL OUTPUTS:
;   mden - total luminosity density in magnitudes of fit
; COMMENTS:
;   The version of the double schechter function we use here is the
;   sum of two Schecter functions with the same exponential cutoff but
;   different normalizations and faint-end slopes. Returns the steeper
;   faint end slope term in PHIPLUS and ALPHAPLUS. 
; REVISION HISTORY:
;   2002-7-04  written - Blanton
;-
function lf_fit_schechter_plus_func,params,x=x,y=y,err=err

phimodel=lf_schechter_plus(x,exp(params[0]),params[1],params[2], $
                           exp(params[3]), params[4])
dev=((y-phimodel)/err)
return,dev

end
;
pro lf_fit_schechter_plus,am,phi,phistddev,schechter_plus,mden=mden

if(n_tags(schechter_plus) eq 0) then $
  schechter_plus={phistar:1.d-2, mstar:-20.d, alpha:-1.d, phiplus:1.d-3, $
                  alphaplus:-1.d, phistar_err:1.d-2, mstar_err:-20.d, $
                  alpha_err:-1.d, phiplus_err:-16.d, alphaplus_err:-1.3d}
start=dblarr(5)
start[0]=alog(schechter_plus.phistar)
start[1]=schechter_plus.mstar
start[2]=schechter_plus.alpha
start[3]=alog(schechter_plus.phiplus)
start[4]=schechter_plus.alphaplus

parinfo=replicate({limited:[0,0], limits:[0.,0.]}, 5)
parinfo[0].limited[0]=1
parinfo[3].limited[0]=1
parinfo[0].limits[0]=-8.
parinfo[3].limits[0]=-8.
functargs={x:am,y:phi,err:phistddev}
params=mpfit('lf_fit_schechter_plus_func',start,functargs=functargs, $
             perror=perror, parinfo=parinfo)

if(params[4] gt params[2]) then begin
    tmpparams=params[0]
    params[0]=params[3]
    params[3]=tmpparams
    tmpparams=params[2]
    params[2]=params[4]
    params[4]=tmpparams
    tmpperror=perror[0]
    perror[0]=perror[3]
    perror[3]=tmpperror
    tmpperror=perror[2]
    perror[2]=perror[4]
    perror[4]=tmpperror
endif

schechter_plus.phistar=exp(params[0])
schechter_plus.mstar=params[1]
schechter_plus.alpha=params[2]
schechter_plus.phiplus=exp(params[3])
schechter_plus.alphaplus=params[4]
schechter_plus.phistar_err=perror[0]*schechter_plus.phistar
schechter_plus.mstar_err=perror[1]
schechter_plus.alpha_err=perror[2]
schechter_plus.phiplus_err=perror[3]*schechter_plus.phiplus
schechter_plus.alphaplus_err=perror[4]
mden=lf_schechter_mden_plus(schechter_plus)

help,/st,schechter_plus

end
