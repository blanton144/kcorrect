;+
; NAME:
;   witt_ext
; PURPOSE:
;   Calculate optical depth from the Witt theoretical results in kcorrect/data
; CALLING SEQUENCE:
;   att=witt_ext(dust_str, tauv, lambda)
; INPUTS:
;   dust_str - structure with string elements:
;                   geometry - 'dusty', 'cloudy' or 'shell'
;                   dust - 'SMC' or 'MW' 
;                   structure - 'c' (clumpy) or 'h' (homogeneous)
;   tauv - optical depth in V
;   lambda - wavelengths (in Angstroms) where A is desired
; OUTPUTS:
;   att - extinction at each lambda
; EXAMPLES:
;   A spectrum spec should have att applied as follows:
;     dusty_spec = spec * exp(-att)
; REVISION HISTORY:
;   04-Jan-2002  Translated to IDL by Mike Blanton, NYU
;   29-Apr-2003  Common block added for speed
;-
;------------------------------------------------------------------------------
function witt_ext, dust_str, tauv, lambda

common witt_ext_common, w_lambda, w_tauv, w_tau_att, w_geometry, w_dust, $
  w_structure

if(tauv eq 0.) then $
  return,lambda*0.

; dust_str should be {geometry, dust, structure}
geometry=dust_str.geometry
dust=dust_str.dust
structure=dust_str.structure

if(n_elements(w_geometry) gt 0) then $
  if(w_geometry ne geometry) then w_lambda=0
if(n_elements(w_dust) gt 0) then $
  if(w_dust ne dust) then w_lambda=0
if(n_elements(w_structure) gt 0) then $
  if(w_structure ne structure) then w_lambda=0
if(n_elements(w_lambda) le 1) then begin
    rootdir=getenv('KCORRECT_DIR')+'/data/dustmodels'
    k_read_ascii_table,w_lambda,rootdir+'/witt.'+geometry+'_'+dust+ $
      '_'+structure+'.lambda.dat'
    k_read_ascii_table,w_tauv,rootdir+'/witt.'+geometry+'_'+dust+ $
      '_'+structure+'.tauv.dat'
    k_read_ascii_table,w_tau_att,rootdir+'/witt.'+geometry+'_'+dust+ $
      '_'+structure+'.tau_att.dat'
    w_geometry=geometry
    w_dust=dust
    w_structure=structure
endif

int_tau_att=fltarr(n_elements(w_lambda))
for i=0L, n_elements(w_lambda)-1L do $
  int_tau_att[i]=interpol(w_tau_att[i,*],w_tauv,tauv)
int_tau_att=[int_tau_att[0], int_tau_att[0], int_tau_att, $
             int_tau_att[n_elements(int_tau_att)-1L], $
             int_tau_att[n_elements(int_tau_att)-1L]]
int_lambda=[0., w_lambda[0], w_lambda, $
            w_lambda[n_elements(w_lambda)-1], 1.e+9]
tau_att=interpol(int_tau_att,int_lambda,lambda)

return,tau_att

end
