;+
; NAME:
;   k_solar_magnitudes
; PURPOSE:
;   calculate the solar magnitudes for bandpasses
; CALLING SEQUENCE:
;   solar_magnitudes= k_solar_magnitudes([band_shift=, filterlist=, $
;                       filterpath=, solarname=]) 
; OPTIONAL INPUTS:
;   filterlist - list of filters (default
;                ['sdss_u0.par','sdss_g0.par','sdss_r0.par','sdss_i0.par', $
;                 'sdss_z0.par'])
;   filterpath - path in which to look for filters
;   band_shift - shift to apply to band passes
;   solarname - name of solar model
;   /silent - shut up
; OUTPUTS:
;   solar_magnitudes - absolute magnitude of sun in specified bands
; EXAMPLES:
;   IDL> print,k_solar_magnitudes()
;    K_READ_BASEL: 1 block(s) of spectra
;     K_PROJECTION_TABLE: Creating rmatrix ...
;     K_PROJECTION_TABLE: Done.
;          6.3800855       5.1210220       4.6432452       4.5322263       4.5107040
;   IDL> print,k_solar_magnitudes(band_shift=0.1)
;    K_READ_BASEL: 1 block(s) of spectra
;     K_PROJECTION_TABLE: Creating rmatrix ...
;     K_PROJECTION_TABLE: Done.
;          6.7792312       5.4336204       4.7553857       4.5749466       4.5155073
; REVISION HISTORY:
;   17-Jan-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function k_solar_magnitudes,band_shift=band_shift,filterpath=filterpath, $
                            filterlist=filterlist, solarname=solarname, $
                            silent=silent

if(n_elements(band_shift) eq 0) then band_shift=0.
if(NOT keyword_set(solarname)) then solarname='lcbsun.ori'
if(NOT keyword_set(filterlist)) then $
  filterlist=['sdss_u0.par','sdss_g0.par','sdss_r0.par','sdss_i0.par', $
              'sdss_z0.par']

; read in the sun and put it at 10pc
k_read_basel,lambda,flux,getenv('KCORRECT_DIR')+'/data/basel/'+solarname, $
  silent=silent
nspectra=n_elements(flux)/n_elements(lambda)
lambda=lambda*10.
pctocm=3.086e+18
solarradtocm=6.960e+10
cspeed=2.99792e+18 ; ang per sec
flux=!DPI*4.*flux*cspeed/ $
  (lambda#replicate(1,nspectra))^2 ; from "hnu" to flambda
flux=flux*(solarradtocm/(10.*pctocm))^2  ; from solar radius to 10 pc

; get answer in maggies
lambda_edges=k_lambda_to_edges(lambda)
maggies=k_project_filters(lambda_edges,flux,filterlist=filterlist, $
                          filterpath=filterpath,band_shift=band_shift, $
                          silent=silent)
solarmags=reform(-2.5*alog10(maggies),n_elements(maggies))

if(n_elements(solarmags) eq 1) then solarmags=solarmags[0]
return,solarmags

end
;------------------------------------------------------------------------------

