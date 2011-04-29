;+
; NAME:
;   lf_simple_vmax
; PURPOSE:
;   calculate a simple vmax, without K-corrections or area
; USAGE:
;   vmax= lf_simple_vmax(absm, mlimit, zlimit [, omega0=, omegal0= ])
; INPUTS:
;   absm - [N] absolute magnitudes
;   mlimit - [2] magnitude limits in observed frame band of absm
;   zlimit - [2] redshift limits
; OPTIONAL INPUTS:
;   omega0           omega_matter to use (default: 0.3)
;   omegal0          omega_lambda to use (default: 0.7)
; OUTPUTS:
;   vmax - [N] structure with: .ZMIN, .ZMAX, .VMAX
; COMMENTS:
;   vmax is returned in h^{-3} Mpc^3 comoving
; REVISION HISTORY:
;   2011-04-26  MRB NYU
;-
function lf_simple_vmax, absm, mlimit, zlimit, omega0=omega0, omegal0=omegal0

; settings
pi=3.14159265358979D
dh=2.99792D+5/100.D

if(n_elements(omega0) eq 0) then omega0=0.3
if(n_elements(omegal0) eq 0) then omegal0=0.7

; for each object...
ngals= n_elements(absm)
vmax= replicate({vmax:0., zmin:0., zmax:0.}, ngals)

nz= 10000L
zgrid= findgen(nz)/float(nz-1L)*(zlimit[1]-zlimit[0])+zlimit[0]
dm= lf_distmod(zgrid, omega0=omega0, omegal0=omegal0)

for i=0L,ngals-1L do begin
    mag= absm[i]+dm
    vmax[i].zmin= ((interpol(zgrid, mag, mlimit[0]))>zlimit[0])<zlimit[1]
    vmax[i].zmax= ((interpol(zgrid, mag, mlimit[1]))>zlimit[0])<zlimit[1]
    vmax[i].vmax= (lf_comvol(vmax[i].zmax, omega0=omega0, omegal0=omegal0)- $
                   lf_comvol(vmax[i].zmin, omega0=omega0, omegal0=omegal0))*4.*!DPI/3.
endfor

return, vmax

end
