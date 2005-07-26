;+
; NAME:
;   lf_t2z
; PURPOSE:
;   transform to redshift from age of universe
; CALLING SEQUENCE:
;   z=lf_t2z(t [, omega0=, omegal0= ])
; INPUTS:
;   t - age of universe (h^{-1} Gyrs)
; OPTIONAL INPUTS:
;   omega0 - matter density (default 0.3)
;   omegal0 - cosmological constant (default 0.7)
; OUTPUTS:
;   z -  redshift
; REVISION HISTORY:
;   2005-7-22  written - Blanton
;-
function lf_t2z, t, omega0=omega0, omegal0=omegal0

if(n_elements(omega0) eq 0) then omega0=0.3
if(n_elements(omegal0) eq 0) then omegal0=0.7

soname=filepath('libkcorrect.'+kcorrect_so_ext(), $
                root_dir=getenv('KCORRECT_DIR'), subdirectory='lib')
ngals=n_elements(t)
z=fltarr(ngals)
retval=call_external(soname, 'idl_t2z', float(t), float(omega0), $
                     float(omegal0), float(z),long(ngals))
return,z

end
