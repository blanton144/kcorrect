;+
; NAME:
;   lf_evolve
; PURPOSE:
;   evolution in magnitudes, given evolution parameters and redshift
; USAGE:
;   evol=lf_evolve(z [, qevolve=, qz0= ])
; INPUTS:
;   z           redshift
; OPTIONAL INPUTS:
;   qevolve     q parameter in q*(z-z0)
;   qz0         z0 parameter in q*(z-z0)
; KEYWORDS:
; OUTPUTS:
; OPTIONAL INPUTS/OUTPUTS:
; BUGS:
; DEPENDENCIES:
;   idlutils
; REVISION HISTORY:
;   2002-11-17  written - Blanton
;-
function lf_evolve,z, qevolve=qevolve, qz0=qz0

evol=0.
if(n_elements(qevolve) gt 0 and n_elements(qz0) gt 0) then $
  evol=evol+qevolve*(z-qz0)

return,evol

end
