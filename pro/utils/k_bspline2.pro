;+
; NAME:
;   k_bspline2
; PURPOSE:
;   return value of a second-order B-spline basis function centered at zero
; CALLING SEQUENCE:
;   vals= k_bspline2(x)
; INPUTS:
;   x - [N] positions
; OUTPUTS:
;   vals - [N] value of basis function at positions
; COMMENTS:
;   A very specialized function. 
; REVISION HISTORY:
;   26-Feb-2003  Written by M. Blanton, NYU
;-  
;------------------------------------------------------------------------------
function k_bspline2, x

val=fltarr(n_elements(x))
indx=where(x gt -1.5 and x le -0.5,count)
if(count gt 0) then val[indx]=0.5*(x[indx]+1.5)^2 
indx=where(x ge 0.5 and x lt 1.5,count)
if(count gt 0) then val[indx]=0.5*(x[indx]-1.5)^2 
indx=where(x gt -0.5 and x lt 0.5,count)
if(count gt 0) then val[indx]=-(x[indx]^2-0.75) 
return, val

end
;------------------------------------------------------------------------------
