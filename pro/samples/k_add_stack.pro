;+
; NAME:
;   k_add_stack
; PURPOSE:
;   adds stack of observations for multiple bands/objects
; CALLING SEQUENCE:
;   k_add_stack, stack, stack_ivar, stacked, stacked_ivar
; INPUTS:
;   stack - [Nbands, Ngalaxies, Nobs] maggies in each band for each
;           galaxy for each observation
;   stack_ivar - [Nbands, Ngalaxies, Nobs] inverse variance for each
;                observation
; OUTPUTS:
;   stacked - [Nbands, Ngalaxies] average maggies
;   stacked_ivar - [Nbands, Ngalaxies] inverse variance of average
; COMMENTS:
;   Weights average by inverse variance.
; REVISION HISTORY:
;   18-Jan-2003  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_add_stack, stack, stack_ivar, stacked, stacked_ivar

nk=(size(stack,/dim))[0]
ng=(size(stack,/dim))[1]
no=(size(stack,/dim))[2]

stacked=fltarr(nk,ng)
stacked_ivar=fltarr(nk,ng)
for i=0L, ng-1L do begin
    weights=stack_ivar[2,i,*]
    for k=0, nk-1 do begin
        ii=where(weights gt 0. and stack_ivar[k,i,*] gt 0., iicount)
        if(iicount gt 0) then begin
            stacked[k,i]=total(weights[ii]*stack[k,i,ii],/double)/ $
              total(weights[ii],/double)
            stacked_ivar[k,i]=(total(weights[ii],/double))^2/ $
              total(weights[ii]^2/stack_ivar[k,i,ii],/double)
        endif
    endfor
endfor

end
