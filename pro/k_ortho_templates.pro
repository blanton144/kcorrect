;+
; NAME:
;   k_ortho_templates
;
; PURPOSE:
;   Orthogonalize SED templates
;
; CALLING SEQUENCE:
;   k_ortho_templates, vmatrix, lambda, bmatrix, bflux
;
; INPUTS:
;   vmatrix   - template value
;   lambda    - template wavelength
;
; OPTIONAL INPUTS:
;   sublmin, sublmax   - limits of subrange to calculate flux for
;
; OUTPUTS:
;   bmatrix   - orthogonalized templates
;   bflux     - flux in specified range 
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
;   04-Jan-2002  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_ortho_templates, vmatrix, lambda, bmatrix, bflux, sublmin=sublmin, sublmax=sublmax

; Need at least 3 parameters
if (N_params() LT 4) then begin
    klog, 'Syntax - k_ortho_templates, vmatrix, lambda, bmatrix, bflux,'
    klog, '     [sublmin=, sublmax=]'
    return
endif

; Set defaults
if (NOT keyword_set(sublmin)) then sublmin=3500.d
if (NOT keyword_set(sublmax)) then sublmax=7500.d
nl=long(n_elements(lambda))-1l
nb=long(n_elements(vmatrix)/nl)

; Set subrange 
subindx=where(lambda ge sublmin and lambda lt sublmax)
subindxp1=subindx+1l

; Orthonormalize
bmatrix=vmatrix
bflux=dblarr(nb)
for b = 0l, nb-1l do begin
  ; orthogonalize
  for bp = 0l, b-1l do begin
	  dot=total(0.5*(lambda[subindxp1]^2-lambda[subindx]^2)* $
							bmatrix[subindx,b]*bmatrix[subindx,bp],/double)
		bmatrix[*,b]=bmatrix[*,b]-dot*bmatrix[*,bp]
  endfor 
	
	; normalize
	dot=total(0.5*(lambda[subindxp1]^2-lambda[subindx]^2)*bmatrix[subindx,b]^2, $
						/double)
	dot=1./sqrt(dot)
	bmatrix[*,b]=bmatrix[*,b]*dot

	; calculate flux
	bflux[b]=total(0.5*(lambda[subindxp1]^2-lambda[subindx]^2) $
								 *bmatrix[subindx,b],/double)
endfor 

end
