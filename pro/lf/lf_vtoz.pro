;+
; NAME:
;   lf_vtoz
; PURPOSE:
;   transform volume to redshift
; CALLING SEQUENCE:
;   z=lf_vtoz(comvol [, omega0=, omegal0= ])
; INPUTS:
;   comvol -  R^3 in equation V=4*PI*R^3/3.
; OPTIONAL INPUTS:
;   omega0 - matter density (default 0.3)
;   omegal0 - cosmological constant (default 0.7)
; OUTPUTS:
;   z - redshift
; REVISION HISTORY:
;   2002-5-22  written - Blanton
;-
function lf_vtoz, v, omega0=omega0, omegal0=omegal0

if(n_elements(omega0) eq 0) then omega0=0.3
if(n_elements(omegal0) eq 0) then omegal0=0.7

soname=filepath('libkcorrect.'+kcorrect_so_ext(), $
                root_dir=getenv('KCORRECT_DIR'), subdirectory='lib')
ngals=n_elements(v)
z=fltarr(ngals)
retval=call_external(soname, 'idl_Vtoz', float(v), float(omega0), $
                     float(omegal0), float(z),long(ngals))
return,z

end
