;+
; NAME:
;   lf_z2t
; PURPOSE:
;   transform redshift to age of universe
; CALLING SEQUENCE:
;   t=lf_z2t(z [, omega0=, omegal0= ])
; INPUTS:
;   z -  redshiift
; OPTIONAL INPUTS:
;   omega0 - matter density (default 0.3)
;   omegal0 - cosmological constant (default 0.7)
; OUTPUTS:
;   t - age of universe (h^{-1} Gyrs)
; REVISION HISTORY:
;   2005-7-22  written - Blanton
;-
function lf_z2t, z, omega0=omega0, omegal0=omegal0

if(n_elements(omega0) eq 0) then omega0=0.3
if(n_elements(omegal0) eq 0) then omegal0=0.7

soname=filepath('libkcorrect.'+kcorrect_so_ext(), $
                root_dir=getenv('KCORRECT_DIR'), subdirectory='lib')
ngals=n_elements(z)
t=fltarr(ngals)
retval=call_external(soname, 'idl_z2t', float(z), float(omega0), $
                     float(omegal0), float(t),long(ngals))
return,t

end
