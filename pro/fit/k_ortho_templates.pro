;+
; NAME:
;   k_ortho_templates
; PURPOSE:
;   Orthogonalize SED templates, scale to constant flux
; CALLING SEQUENCE:
;   k_ortho_templates, coeffs, vmatrix, lambda, bcoeffs, bmatrix, $
;      bflux [, bdotv=, sublmin=, sublmax= ] $
; INPUTS:
;   coeffs    - coefficients (so they can be transformed with the vmatrix)
;   vmatrix   - template value
;   lambda    - template wavelength
; OPTIONAL INPUTS:
;   sublmin, sublmax   - limits of subrange to calculate flux for
; OUTPUTS:
;   bcoeffs   - new coeffs for orthogonalized templates
;   bmatrix   - orthogonalized templates
;   bflux     - flux in specified range 
;   bdotv     - transformation between new and old coords
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   04-Jan-2002  Written by Mike Blanton, NYU
;   06-May-2003  Altered for new system by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_orthogonalize, lambda, inmatrix, outmatrix, outflux, subindx, subindxp1

nl=long(n_elements(lambda))-1L
nv=long(n_elements(inmatrix)/nl)

dinmatrix=double(inmatrix)
doutmatrix=dblarr(nl,nv)

outflux=fltarr(nv)
for b=0L, nv-1L do begin
    doutmatrix[*,b]=dinmatrix[*,b]

;   orthogonalize
    for bp = 0l, b-1l do begin
        dot=total((lambda[subindxp1]-lambda[subindx])* $
                  doutmatrix[subindx,b]*doutmatrix[subindx,bp],/double)
        doutmatrix[*,b]=doutmatrix[*,b]-dot*doutmatrix[*,bp]
    endfor 
    
;   normalize
    dot=total((lambda[subindxp1]-lambda[subindx])*doutmatrix[subindx,b]^2, $
              /double)
    dot=1./sqrt(dot)
    doutmatrix[*,b]=doutmatrix[*,b]*dot
    
;   calculate flux
    outflux[b]=$
      total((lambda[subindxp1]-lambda[subindx])*doutmatrix[subindx,b],/double)
endfor 

outmatrix=float(doutmatrix)

end
;
pro k_ortho_templates, coeffs, vmatrix, lambda, bcoeffs, bmatrix, bflux, $
                       sublmin=sublmin, sublmax=sublmax, bdotv=bdotv, $
                       bdotb=bdotb, vdotv=vdotv

; Need at least 3 parameters
if (N_params() LT 4) then begin
    print, 'Syntax - k_ortho_templates, coeffs, vmatrix, lambda, bcoeffs, bmatrix, $' 
    print, '     bflux [, sublmin=, sublmax=, bdotv=, bdotb= ]'
    return
endif

; Set defaults
if (NOT keyword_set(sublmin)) then sublmin=1500.d
if (NOT keyword_set(sublmax)) then sublmax=30000.d
nl=long(n_elements(lambda))-1L
nv=long(n_elements(vmatrix)/nl)

; Set subrange 
subindx=where(lambda ge sublmin and lambda lt sublmax)
subindxp1=subindx+1L

; orthogonalize
bmatrix=fltarr(nl,nv)
k_orthogonalize, lambda, vmatrix, bmatrix, bflux, subindx, subindxp1

; set the first vector to carry all the flux
tmpb=bmatrix#bflux
bmatrix[*,0]=tmpb
k_orthogonalize, lambda, bmatrix, bmatrix, bflux, subindx, subindxp1

; calculate transformation of coordinates and apply it
bdotv=fltarr(nv,nv)
vdotv=fltarr(nv,nv)
bdotb=fltarr(nv,nv)
for b = 0l, nv-1l do $
  for bp = 0l, nv-1l do $
  bdotv[b,bp]=total((lambda[subindxp1]-lambda[subindx])* $
                    bmatrix[subindx,b]*vmatrix[subindx,bp],/double)* $
  (sublmax-sublmin)
for b = 0l, nv-1l do $
  for bp = 0l, nv-1l do $
  vdotv[b,bp]=total((lambda[subindxp1]-lambda[subindx])* $
                    vmatrix[subindx,b]*vmatrix[subindx,bp],/double)* $
  (sublmax-sublmin)
for b = 0l, nv-1l do $
  for bp = 0l, nv-1l do $
  bdotb[b,bp]=total((lambda[subindxp1]-lambda[subindx])* $
                    bmatrix[subindx,b]*bmatrix[subindx,bp],/double)* $
  (sublmax-sublmin)
bcoeffs=bdotv#coeffs

end
