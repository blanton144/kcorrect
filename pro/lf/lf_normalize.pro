;+
; NAME:
;   lf_normalize
; PURPOSE:
;   For each galaxy, calculate the probability of it being in the
;   observable range
; INPUTS:
;   zz               [N] redshifts
;   absmag           [N] absolute magnitudes
;   absmerr          [N] absolute magnitude errors
;   absmmin          [N] minimum abs mag for each object
;   absmmax          [N] maximum abs mag for each object
;   nphi             number of "bins"
;   sample_absmmin   minimum abs mag for the sample
;   sample_absmmax   maximum abs mag for the sample
;   sigabsmag        width of gaussian (recommended to be
;                    0.7*(sample_absmmax-sample_absmmin)/nphi)
; OPTIONAL INPUTS:
;   zzero            zeropoint in redshift (default 0.1)
; KEYWORDS:
; OUTPUTS:
;   phi              amplitude of bin
;   q                evolution of lum
;   p                evolution of num
; OPTIONAL OUTPUTS:
; BUGS:
;   TODO: 
;         finish construction of likelihood
;         create test set of data
;         test the data ...
; DEPENDENCIES:
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
