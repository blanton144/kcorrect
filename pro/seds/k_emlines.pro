;+
; NAME:
;   k_emlines
; COMMENTS:
;   returns set of generic emission line templates
; CALLING SEQUENCE:
;   k_emlines, loglam, templates
; INPUTS:
;   loglam - [N] desired wavelength grid
;   templates - [N, Nem] templates 
; COMMENTS:
;   Templates are gaussians with 300 km/s dispersion
; REVISION HISTORY:
;   13-Jun-2006  Michael Blanton (NYU)
;-
;------------------------------------------------------------------------------
pro k_emlines, loglam, templates, vdisp=vdisp, bluest=bluest

if (NOT keyword_set(vdisp)) then vdisp=300.
sigma=vdisp/(2.99792e+5*alog(10.))  ;; smoothing sigma in log lambda

nper=[1, $
      2, $
      2, $
      1]

tcen=[ [ 6564.61, 0.00], $
       [ 4862.683, 6564.61], $
       [ 4960.29, 5008.24], $
       [ 3727.09, 0.00] ]

rcen=[ [1.0, 0.], $
       [1.0, 2.7], $
       [1.0, 3.0], $
       [1.0, 0.]]

bluest=tcen[0,*]

templates=fltarr(n_elements(loglam), n_elements(nper))
for i=0L, n_elements(nper)-1L do begin
    for j=0L, nper[i]-1L do begin
        templates[*, i]= templates[*,i]+ $
          rcen[j,i]*exp(-0.5*(loglam-alog10(tcen[j,i]))^2/(sigma)^2)/ $
          (sqrt(2.*!DPI)*sigma)/total(rcen[*,i])
    endfor
    tmax=max(templates[*,i])
    templates[*,i]=templates[*,i]/tmax
endfor

end
