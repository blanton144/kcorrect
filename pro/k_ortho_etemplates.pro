;+
; NAME:
;   k_ortho_etemplates
;
; PURPOSE:
;   Orthogonalize eigentemplates for SED fitting
;
; CALLING SEQUENCE:
;   k_ortho_etemplates, ematrix, bflux
;
; INPUTS:
;   ematrix   - eigentemplate value
;   bflux     - flux in each SED template
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL INPUT/OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   05-Jan-2002  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_ortho_etemplates, ematrix, bflux

; Need at least 2 parameters
if (N_params() LT 2) then begin
    klog, 'Syntax - k_ortho_etemplates, ematrix, bflux'
    return 
endif

nb=long(n_elements(bflux))
nt=long(n_elements(ematrix)/nb)

; Perform initial orthogonalization
for t = 0l, nt-1l do begin
  ; orthogonalize
  for tp = 0l, t-1l do begin
	  dot=total(ematrix[*,t]*ematrix[*,tp],/double)
		ematrix[*,t]=ematrix[*,t]-dot*ematrix[*,tp]
  endfor 
	
	; normalize
	dot=total(ematrix[*,t]^2,/double)
	dot=1./sqrt(dot)
	ematrix[*,t]=ematrix[*,t]*dot
endfor 

; Calculate flux in each eigenspectrum
eflux=bflux#ematrix

; Set the initial vector to carry all flux
tmpe=ematrix#transpose(eflux)
ematrix[*,0]=tmpe

; Perform final orthogonalization
for t = 0l, nt-1l do begin
  ; orthogonalize
  for tp = 0l, t-1l do begin
	  dot=total(ematrix[*,t]*ematrix[*,tp],/double)
		ematrix[*,t]=ematrix[*,t]-dot*ematrix[*,tp]
  endfor 
	
	; normalize
	dot=total(ematrix[*,t]^2,/double)
	dot=1./sqrt(dot)
	ematrix[*,t]=ematrix[*,t]*dot
endfor 

; Now normalize the ematrix so that unit coefficient 
; means unit visual flux
eflux=bflux#ematrix
ematrix=ematrix/eflux[0]

end 
