;+
; NAME:
;   k_evolve
; PURPOSE:
;   returns "evolved" absolute magnitude using simple evolution model
; CALLING SEQUENCE:
;   absm_evolved= k_evolve(absm, z, q0, q1, qz0)
; INPUTS:
;   absm - absolute magnitude
;   z - redshift
;   q0, q1, qz0 - parameters in formula:
;     absm_evolved = absm+(q0*(1.+q1*(z-qz0)))*(z-qz0)
; REVISION HISTORY:
;   2005-3-17  written - Blanton
;-
function k_evolve, absm, z, q0, q1, qz0
return, absm+(q0*(1.+q1*(z-qz0)))*(z-qz0)
end
