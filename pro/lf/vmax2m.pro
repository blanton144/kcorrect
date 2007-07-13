;+
; NAME:
;   vmax2m
; USAGE:
;   mass= vmax2m(vmax)
; INPUTS:
;   vmax - maximum circular velocity (in km/s)
; OUTPUTS:
;   mass - mass of halo (in Msolar/h)
; COMMENTS:
;   Inverts m2vmax (good between 10^3 and 10^18 solar masses)
; REVISION HISTORY:
;   2002-11-01  written - Blanton
;-
function vmax2m, in_vmax, littleh=littleh, omega0=omega0

common com_vmax2m, lmgrid, lvgrid, glittleh, gomega0

if(NOT keyword_set(littleh)) then littleh=0.7
if(NOT keyword_set(omega0)) then omega0=0.24
if(NOT keyword_set(glittleh)) then glittleh=littleh
if(NOT keyword_set(gomega0)) then gomega0=omega0

if(n_elements(lmgrid) eq 0 OR $
   glittleh ne littleh OR $
   gomega0 ne omega0) then begin
    gomega0=omega0
    glittleh=littleh
    mmin=1.d+3
    mmax=1.d+18
    ngrid=10000L
    lmgrid=alog10(mmin)+(alog10(mmax)-alog10(mmin))*(findgen(ngrid)+0.5)/ $
      float(ngrid)
    lvgrid=alog10(m2vmax(10.^lmgrid, littleh=littleh, omega0=omega0))
endif

return, 10.^interpol(lmgrid, lvgrid, alog10(in_vmax))

end
