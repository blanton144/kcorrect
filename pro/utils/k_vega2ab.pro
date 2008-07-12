;+
; NAME:
;   k_vega2ab
; PURPOSE:
;   Calculate conversion such that m_AB = m_Vega + k_vega2ab()
; CALLING SEQUENCE:
;   vega2ab=k_vega2ab([ filterlist=, filterpath=, band_shift=, $
;                       /kurucz, /hayes])
; INPUTS:
; OPTIONAL INPUTS:
;   filterlist - list of .par files with filters ([default: 
;                ['sdss_u0.par', 'sdss_g0.par', 'sdss_r0.par',
;                'sdss_i0.par', 'sdss_z0.par'])
;   filterpath - path where to find them
;   band_shift - blueward shift of bandpass (factor 1.+band_shift)
;   /kurucz - use the Kurucz 1991 model spectrum 
;   /hayes - use the Hayes spectrophotometry
; OUTPUTS:
;   vega2ab - term to add to a Vega magnitude in given band to convert
;             to an AB magnitude
; COMMENTS:
;   We can compare the results of this conversion with the results of
;   Bessell (1993) in IAU Colloq. 136: Stellar Photometry - Current
;   Techniques and Future Developments for the Bessell bandpasses:
;
;    IDL> print,3631.*10.^(-0.4*k_vega2ab(filterlist=['bessell_B.par','bessell_V.par','bessell_R.par','bessell_I.par'],/kurucz))/[4000.,3600.,3060.,2420.]
;    K_READ_BASEL: 1 block(s) of spectra
;     K_PROJECTION_TABLE: Creating rmatrix ...
;     K_PROJECTION_TABLE: Done.
;         0.98735007      0.98666174      0.97275494      0.98666968
;   IDL> print,3631.*10.^(-0.4*k_vega2ab(filterlist=['bessell_B.par','bessell_V.par','bessell_R.par','bessell_I.par'],/hayes))/[4000.,3600.,3060.,2420.] 
;     K_PROJECTION_TABLE: Creating rmatrix ...
;     K_PROJECTION_TABLE: Done.
;          1.0158529      0.98916833      0.97277257      0.98920956
;   IDL> $date
;   Tue Apr 29 20:31:21 EDT 2003
;
;   The first version uses the Kurucz 1991 models from the BaSeL
;   distribution (ADC: J/A+AS/125/229); the second the Hayes 
;   spectrophotometry, as transcribed by Hogg (see details in the 
;   data/filters/hoggraw/hayes directory). 
; PROCEDURES CALLED:
;    k_load_filters
;    k_project_filters
;    readcol (idlutils)
; REVISION HISTORY:
;   17-Jan-2002  Translated to IDL by Mike Blanton, NYU
;   10-Jul-2008  Moustakas - added SILENT keyword
;-
;------------------------------------------------------------------------------
function k_vega2ab,filterpath=filterpath, $
                   filterlist=filterlist, $
                   band_shift=band_shift, $
                   kurucz=kurucz,hayes=hayes,silent=silent

if(NOT keyword_set(filterlist)) then $
  filterlist=['sdss_u0.par','sdss_g0.par','sdss_r0.par','sdss_i0.par', $
              'sdss_z0.par']

; read in the vega 
if(keyword_set(kurucz)) then begin
    veganame='lcbvega.ori'
    k_read_basel,lambda,flux,getenv('KCORRECT_DIR')+'/data/basel/'+veganame,silent=silent
    nspectra=1L
    lambda=lambda*10.e
    cspeed=2.99792e+18          ; ang per sec
    flux=!DPI*4.*flux*cspeed/(lambda)^2 ; from "hnu" to flambda
    flux=flux*6.5043898e-17     ; normalize (to fit Hayes)
endif 

if(keyword_set(hayes)) then begin
    readcol,getenv('KCORRECT_DIR')+'/data/filters/hoggraw/hayes/hayes.txt', $
      lambda,lflux,comment='#',/silent
    flux=10.^(-0.4*lflux)*4.65e-9
endif 

if(NOT keyword_set(flux)) then begin
    klog,'ERROR: Must specify a particular determination of Vega flux'
    return,-100000.
end

; get AB maggies of Vega -- this is the conversion
lambda_edges=k_lambda_to_edges(lambda)
maggies=k_project_filters(lambda_edges,flux,filterlist=filterlist, $
                          filterpath=filterpath,band_shift=band_shift,silent=silent)
vega2ab=reform(-2.5*alog10(maggies),n_elements(filterlist))

if(n_elements(vega2ab) eq 1) then vega2ab=vega2ab[0]
return,vega2ab

end
;------------------------------------------------------------------------------

