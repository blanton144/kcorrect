;+
; NAME:
;   lf_schechter
; PURPOSE:
;   given schechter parameters and absolute mags, return vals
; CALLING SEQUENCE:
;   vals= lf_schechter(absmag, phistar, mstar, alpha)
;     OR
;   vals= lf_schechter(absmag, schechter)
; INPUTS:
;   absmag - [N] set of absolute magnitudes
;   AND:
;    phistar - phi* parameter in schechter defn
;    mstar - M* parameter in schechter defn
;    alpha - alpha parameter ('faint end slope') in schechter defn
;   OR:
;   schechter - structure with phistar, mstar, alpha entries
; REVISION HISTORY:
;   20-Oct-2003  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function lf_schechter,m,phistar,mstar,alpha

if(n_tags(phistar) gt 0) then begin
    use_phistar=phistar.phistar
    use_mstar=phistar.mstar
    use_alpha=phistar.alpha
endif else begin
    use_phistar=phistar
    use_mstar=mstar
    use_alpha=alpha
endelse

val=0.4*alog(10.)*use_phistar*10.^(-0.4*(m-use_mstar)*(use_alpha+1.))* $
  exp(-10.^(-0.4*(m-use_mstar)))
return,val

end
