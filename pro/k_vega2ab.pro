;+
; NAME:
;   k_vega2ab
;
; PURPOSE:
;   Calculate conversion to apply to AB magnitudes to obtain
;   Vega magnitudes.
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
function k_vega2ab,filterpath=filterpath, $
                   filterlist=filterlist, $
                   vega_z=vega_z, $
                   veganame=veganame, $
                   qa=qa,kurucz=kurucz,hayes=hayes

if(n_elements(ab_z) eq 0) then ab_z=0.
if(n_elements(vega_z) eq 0) then vega_z=0.
if(NOT keyword_set(veganame)) then veganame='lcbvega.ori'
if(NOT keyword_set(filterpath)) then $
  filterpath=getenv('KCORRECT_DIR')+'/data/filters'
if(NOT keyword_set(filterlist)) then $
  filterlist=['sdss_u0','sdss_g0','sdss_r0','sdss_i0','sdss_z0']

; read in the vega 
k_read_basel,lambda,flux,getenv('KCORRECT_DIR')+'/data/basel/'+veganame
nspectra=1L
lambda=lambda*10.d
cspeed=2.99792d+18 ; ang per sec
flux=3.14159*4.*flux*cspeed/(lambda)^2 ; from "hnu" to flambda
flux=flux/1.5d+16

; get AB maggies of Vega -- this is the conversion
maggies=k_project_filters(lambda,flux,filterlist=filterlist, $
                          filterpath=filterpath,band_shift=vega_z,qa=qa)
vega2ab=-2.5*alog10(maggies)

return,vega2ab

end
;------------------------------------------------------------------------------

