;+
; NAME:
;   k_solar_magnitudes
;
; PURPOSE:
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL INPUT/OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   17-Jan-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function k_solar_magnitudes,band_shift=band_shift,filterpath=filterpath, $
                            filterlist=filterlist, $
                            solarname=solarname

if(n_elements(band_shift) eq 0) then band_shift=0.
if(NOT keyword_set(solarname)) then solarname='lcbsun.ori'
if(NOT keyword_set(filterpath)) then $
  filterpath=getenv('KCORRECT_DIR')+'/data/filters'
if(NOT keyword_set(filterlist)) then $
  filterlist=['sdss_u0','sdss_g0','sdss_r0','sdss_i0','sdss_z0']

; read in the sun and put it at 10pc
k_read_basel,lambda,flux,getenv('KCORRECT_DIR')+'/data/basel/'+solarname
nspectra=n_elements(flux)/n_elements(lambda)
lambda=lambda*10.
pctocm=3.086d+18
solarradtocm=6.960d+10
cspeed=2.99792d+18 ; ang per sec
flux=3.14159*4.*flux*cspeed/ $
  (lambda#replicate(1,nspectra))^2 ; from "hnu" to flambda
flux=flux*(solarradtocm/(10.*pctocm))^2  ; from solar radius to 10 pc
lambda=lambda ; redshift it
flux=flux ; redshift it

; get answer in maggies
maggies=k_project_filters(lambda,flux,filterlist=filterlist, $
                          filterpath=filterpath,band_shift=band_shift)
solarmags=-2.5*alog10(maggies)

return,solarmags

end
;------------------------------------------------------------------------------

