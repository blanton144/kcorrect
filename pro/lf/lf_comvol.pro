;+
; NAME:
;   lf_comvol
; PURPOSE:
;   Calculate the comoving volume out to a certain redshift
; CALLING SEQUENCE:
;   comvol=lf_comvol(z [, omega0=, omegal0= ])
; INPUTS:
;   z - redshift
; OPTIONAL INPUTS:
;   omega0 - matter density (default 0.3)
;   omegal0 - vacuum energy density (default 0.7)
; COMMENTS:
;   Uses the ztransform routines
;   returns R^3 in equation V=4*PI*R^3/3.
;   (in units of h^{-3} Mpc^3)
; REVISION HISTORY:
;   2002-5-22  written - Blanton
;-
function lf_comvol, z, omega0=omega0, omegal0=omegal0

if(n_elements(omega0) eq 0) then omega0=0.3
if(n_elements(omegal0) eq 0) then omegal0=0.7

soname=filepath('libkcorrect.'+kcorrect_so_ext(), $
                root_dir=getenv('KCORRECT_DIR'), subdirectory='lib')
ngals=n_elements(z)
v=fltarr(ngals)
retval=call_external(soname, 'idl_ztoV', float(z), float(omega0), $
                     float(omegal0), float(v),long(ngals))
return,v

end
