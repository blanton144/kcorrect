;+
; NAME:
;   k_mkspec_qso
; PURPOSE:
;   Given a wavelength scale, interpolate idlspec2d QSO template onto it
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
;   k_load_ascii_table
;   read_peg (in eplusa)
; REVISION HISTORY:
;   25-Jul-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_mkspec_qso, qsovmatrix, lambda

interp_lambda=0.5*(lambda[0L:n_elements(lambda)-2L]+ $
                   lambda[1L:n_elements(lambda)-1L])
qsotemplate= $
  mrdfits(getenv('IDLSPEC2D_DIR')+'/templates/spEigenQSO-52223.fits',0,hdr)
hdrstruct=hdr2struct(hdr)
lambdatemplate=10.^(hdrstruct.coeff0+dindgen(n_elements(qsotemplate[*,0]))* $
                    hdrstruct.coeff1)
qsovmatrix=interpol(qsotemplate[*,0],lambdatemplate,interp_lambda)

end
;------------------------------------------------------------------------------
