;+
; NAME:
;   lf_distmod
; PURPOSE:
;   Calculate distance modulus
; USAGE:
;   dm=lf_distmod(z [, omega0=, omegal0=])
; INPUTS:
;   z               [N] redshifts
; OPTIONAL INPUTS:
;   omega0           omega_matter to use (default: 0.3)
;   omegal0          omega_lambda to use (default: 0.7)
; OUTPUTS:
;   dm - distance modulus
; REVISION HISTORY:
;   2002-11-17  written - Blanton
;-
function lf_distmod,z,omega0=omega0,omegal0=omegal0

if(n_elements(omega0) eq 0) then omega0=0.3
if(n_elements(omegal0) eq 0) then omegal0=0.7

soname=filepath('libkcorrect.'+kcorrect_so_ext(), $
                root_dir=getenv('KCORRECT_DIR'), subdirectory='lib')
ngals=n_elements(z)
dm=fltarr(ngals)
retval=call_external(soname, 'idl_z2dm', float(z), float(omega0), $
                     float(omegal0), float(dm), long(ngals))

if(n_elements(dm) eq 1) then dm=dm[0]

return,dm

end
