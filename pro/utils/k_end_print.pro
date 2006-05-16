;+
; NAME:
;   k_end_print
; PURPOSE:
;   end a postscript file and return plotting device to X
; CALLING SEQUENCE:
;   k_end_print [, pold=, xold=, yold= ]
; OPTIONAL INPUTS:
;   pold, xold, yold - if these exist, resets !P, !X, !Y to these values
; COMMENTS:
;   A very specialized function. 
; REVISION HISTORY:
;   26-Feb-2003  Written by M. Blanton, NYU
;-  
;------------------------------------------------------------------------------
pro k_end_print,pold=pold,xold=xold,yold=yold

common com_k_print, bangp, bangx, bangy

device,/close

if(keyword_set(bangp)) then !P=bangp
if(keyword_set(bangx)) then !X=bangx
if(keyword_set(bangy)) then !Y=bangy
if(keyword_set(pold)) then !P=pold
if(keyword_set(xold)) then !X=xold
if(keyword_set(yold)) then !Y=yold

set_plot,'x'

end
