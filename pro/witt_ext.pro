;+
; NAME:
;   witt_ebv
;
; PURPOSE:
;   Calculate E(B-V) from the Witt theoretical results
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
;   k_fit_photoz
;
; REVISION HISTORY:
;   04-Jan-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function witt_ext, dust_str, tauv, lambda, w_lambda=w_lambda, $
                   w_tauv=w_tauv,w_tau_att=w_tau_att

if(tauv eq 0.) then $
  return,lambda*0.

; dust_str should be {geometry, dust, structure}
geometry=dust_str.geometry
dust=dust_str.dust
structure=dust_str.structure

if(n_elements(w_lambda) le 1) then begin
    rootdir=getenv('KCORRECT_DIR')+'/data/dustmodels'
    k_load_ascii_table,w_lambda,rootdir+'/witt.'+geometry+'_'+dust+ $
      '_'+structure+'.lambda.dat'
    k_load_ascii_table,w_tauv,rootdir+'/witt.'+geometry+'_'+dust+ $
      '_'+structure+'.tauv.dat'
    
    k_load_ascii_table,w_tau_att,rootdir+'/witt.'+geometry+'_'+dust+ $
      '_'+structure+'.tau_att.dat'
endif

int_tau_att=dblarr(n_elements(w_lambda))
for i=0L, n_elements(w_lambda)-1L do begin
    int_tau_att[i]=interpol(w_tau_att[i,*],w_tauv,tauv)
endfor
tau_att=interpol(int_tau_att,w_lambda,lambda)

return,tau_att

end
