;+
; NAME:
;   k_spec_pegase
; PURPOSE:
;   Given a time result for PEGASE, return a sensible
;   spectrum. Currently adds the lines in an incredibly 
;   stupid manner.
; CALLING SEQUENCE:
; INPUTS:
; OPTIONAL INPUTS:
; KEYWORDS:
; OUTPUTS:
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   25-Jul-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_spec_pegase,peg,spec,lambda,nl=nl,lmin=lmin,lmax=lmax, $
                  linewidth=linewidth,nolines=nolines,nocontinuum=nocontinuum

if(n_elements(nl) eq 0) then nl=5000L
if(n_elements(lmin) eq 0) then lmin=double(min(peg.cont))
if(n_elements(lmax) eq 0) then lmax=double(max(peg.cont))
if(n_elements(linewidth) eq 0) then linewidth=10.

lambda=lmin+(lmax-lmin)*dindgen(nl+1)/double(nl)
dl=lambda[1]-lambda[0]
interp_lambda=0.5*(lambda[0:nl-1]+lambda[1:nl])

if(NOT keyword_set(nocontinuum)) then begin
    spec=interpol(peg.contlum,double(peg.cont),interp_lambda)
endif else begin
    spec=dblarr(n_elements(interp_lambda))
endelse
    
if(NOT keyword_set(nolines)) then begin
    for i=0L, n_elements(peg.lines)-1L do begin
        indx=where(interp_lambda gt peg.lines[i]-3.*linewidth and $
                   interp_lambda lt peg.lines[i]+3.*linewidth,count)
        if(count gt 0) then begin
            tmp_line=exp(-0.5*((interp_lambda[indx]-double(peg.lines[i]))/ $
                               linewidth)^2)/(sqrt(2.*!DPI)*linewidth)
            spec[indx]=spec[indx]+peg.linelum[i]*tmp_line
        endif
    endfor 
endif

end
;------------------------------------------------------------------------------
