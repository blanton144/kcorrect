;+
; NAME:
;   lf_normalize
; PURPOSE:
;   normalize npgauss function so integral is unity
; INPUTS:
;   absmk - centers of gaussians
;   phi - heights of gaussians (changed on output)
;   sample_absmmin, sample_absmmax - abs. mag. limits
;   sigabsmag - width of gaussian (recommended to be
;               0.7*(sample_absmmax-sample_absmmin)/nphi)
; COMMENTS:
;   Changes phi to make integral unity.
; REVISION HISTORY:
;   2002-5-22  written - Blanton
;-
pro lf_normalize,absmk,phi,sample_absmmin,sample_absmmax,sigabsmag

; normalize it appropriately
totalphi=total(phi*0.5* $
               (errorf((sample_absmmax-absmk)/(sqrt(2.)*sigabsmag))- $
                errorf((sample_absmmin-absmk)/(sqrt(2.)*sigabsmag))),/double)
phi=phi/totalphi

end
