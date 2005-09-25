;+
; NAME:
;   k_binspec
; PURPOSE:
;   bins spectrum by integrating over pixels to convert flux density to flux
; CALLING SEQUENCE:
;   newspec= k_rebin(lambda, spec, newlambda)
; INPUTS:
;   lambda - [N] wavelength
;   spec - [N] flux density (in same units as wavelength)
;   newlambda - [Nout] centers of new pixels
; OUTPUTS:
;   newspec - [Nout] total flux in each pixel
; COMMENTS: 
;   Uses k_lambda_to_edges to infer pixel edges from centers.
; REVISION HISTORY:
;   21-Sept-2005  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function k_binspec, lambda, spec, newlambda
                  
edges=k_lambda_to_edges(newlambda)

; Set source object name
soname=filepath('libkcorrect.'+kcorrect_so_ext(), $
                root_dir=getenv('KCORRECT_DIR'), subdirectory='lib')

nl=n_elements(lambda)
nnewl=n_elements(newlambda)
newspec=fltarr(nnewl)
retval=call_external(soname, 'idl_k_binspec', float(lambda), $
                     float(spec), float(edges), float(newspec), $
                     long(nl), long(nnewl))

return, newspec

end
;------------------------------------------------------------------------------
