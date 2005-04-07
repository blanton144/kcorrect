;+
; NAME:
;   lf_erf
; PURPOSE:
;   calculate error function
; USAGE:
;   y= lf_erf(x)
; INPUTS:
;   x - input value
; OUTPUTS:
;   y - erf(x)
; COMMENTS:
;   Calculates as 1.-lf_erfc()
; REVISION HISTORY:
;   2003-10-20  written - Blanton
;-
function lf_erf,x

val=1.-lf_erfc(x)
indx=where(x lt 0.,count)
if(count gt 0) then val[indx]=-val[indx]
return,val

end
