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
; We can compare the results of this conversion with the results of
; Bessell (1993) in IAU Colloq. 136: Stellar Photometry - Current
; Techniques and Future Developments for the Bessell bandpasses:
;
; IDL> $date
; Tue Aug 20 13:31:12 EDT 2002
; IDL> print,3631.*10.^(-0.4*k_vega2ab(filterlist=['bessell_B','bessell_V','bessell_R','bessell_I'],/kurucz))/[4000.,3600.,3060.,2420.] 
;   K_READ_BASEL: 1 block(s) of spectra
;       0.98953211      0.98817475      0.97688385      0.98599360
; IDL> print,3631.*10.^(-0.4*k_vega2ab(filterlist=['bessell_B','bessell_V','bessell_R','bessell_I'],/hayes))/[4000.,3600.,3060.,2420.] 
;        1.0136841      0.99102008      0.97584356      0.99039233
;
; The first version uses the Kurucz models from the BaSeL
; distribution; the second the Hayes spectrophotometry, as transcribed
; by Hogg (see details in the data/filters/hoggraw/hayes directory). 
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
if(keyword_set(kurucz)) then begin
    k_read_basel,lambda,flux,getenv('KCORRECT_DIR')+'/data/basel/'+veganame
    nspectra=1L
    lambda=lambda*10.d
    cspeed=2.99792d+18          ; ang per sec
    flux=3.14159*4.*flux*cspeed/(lambda)^2 ; from "hnu" to flambda
    flux=flux*6.5043898d-17     ; normalize (to fit Hayes)
endif 

if(keyword_set(hayes)) then begin
    readcol,getenv('KCORRECT_DIR')+'/data/filters/hoggraw/hayes/hayes.txt', $
      lambda,lflux,comment='#',/silent
    flux=10.^(-0.4*lflux)*4.65d-9
endif 

if(NOT keyword_set(flux)) then begin
    klog,'ERROR: Must specify a particular determination of Vega flux'
    return,-100000.
end

; get AB maggies of Vega -- this is the conversion
maggies=k_project_filters(lambda,flux,filterlist=filterlist, $
                          filterpath=filterpath,band_shift=vega_z,qa=qa)
vega2ab=-2.5*alog10(maggies)

return,vega2ab

end
;------------------------------------------------------------------------------

