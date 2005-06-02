;+
; NAME:
;   lf_schechter_plus
; PURPOSE:
;   given schechter parameters and absolute mags, return vals
; CALLING SEQUENCE:
;   vals= lf_schechter(absmag, phistar, mstar, alpha, phiplus, alphaplus)
;     OR
;   vals= lf_schechter(absmag, schechter)
; INPUTS:
;   absmag - [N] set of absolute magnitudes
;   AND:
;    phistar - phi* parameter in 1st schechter fn
;    mstar - M* parameter in both schechter fns
;    alpha - alpha parameter ('faint end slope') in 1st schechter fn
;    phiplus - phi* parameter in 2nd schechter defn
;    alphaplus - alpha parameter ('faint end slope') in 2nd schechter fn
;   OR:
;    schechter_plus - structure containing double schechter pars 
;                     (.PHISTAR, .MSTAR, .ALPHA, .PHIPLUS, .ALPHA_PLUS)
; COMMENTS:
;   The double Schechter function is the sum of two Schechter
;   functions with the same MSTAR but different PHISTAR and ALPHA
;   values.
; REVISION HISTORY:
;   20-Oct-2003  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function lf_schechter_plus,m,phistar,mstar,alpha,phiplus,alphaplus

if(n_tags(phistar) gt 0) then begin
    use_phistar=phistar.phistar
    use_mstar=phistar.mstar
    use_alpha=phistar.alpha
    use_phiplus=phistar.phiplus
    use_alphaplus=phistar.alphaplus
endif else begin
    use_phistar=phistar
    use_mstar=mstar
    use_alpha=alpha
    use_phiplus=phiplus
    use_alphaplus=alphaplus
endelse

val=0.4*alog(10.)* $
  use_phistar* $
  10.^(-0.4*(m-use_mstar)*(use_alpha+1.))*exp(-10.^(-0.4*(m-use_mstar)))+ $
  use_phiplus* $
  10.^(-0.4*(m-use_mstar)*(use_alphaplus+1.))*exp(-10.^(-0.4*(m-use_mstar)))
return, val

end
