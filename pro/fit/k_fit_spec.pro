;+
; NAME:
;   k_fit_spec
; PURPOSE:
;   fit a spectrum to the sum of templates 
; CALLING SEQUENCE:
;   k_fit_spec, flux, ivar, coeffs [, vname=, vdisp= ]
; INPUTS:
;   flux - [nl] fluxes (aligned with wavelength grid of models)
;   ivar - [nl] inverse variances (aligned with wavelength grid of models)
; OPTIONAL INPUTS:
;   vname - name of fit to use (default 'default')
;   vdisp - velocity dispersion to smooth fit to
; OUTPUTS:
;   coeffs - [nt] coefficients fit to each template 
;   templates - [nl, nt] template fluxes
; REVISION HISTORY:
;   21-Apr-2005  MRB, NYU
;-
;------------------------------------------------------------------------------
pro k_fit_spec, flux, ivar, coeffs, vname=vname, vdisp=vdisp, $
                templates=templates

k_reconstruct_spec, dum, loglam, /init, vname=vname, nt=nt

nl=n_elements(loglam)
templates=fltarr(nl, nt)
for i=0L, nt-1L do begin
    coeffs=fltarr(nt)
    coeffs[i]=1.
    k_reconstruct_spec, coeffs, loglam, tmpflux, vname=vname, vdisp=vdisp
    templates[*,i]=tmpflux
endfor

if(NOT keyword_set(verbose)) then verbose=1L
if(NOT keyword_set(maxiter)) then maxiter=500000
if(NOT keyword_set(tolerance)) then tolerance=1.e-6

; Set source object name
soname=filepath('libkcorrect.'+kcorrect_so_ext(), $
                root_dir=getenv('KCORRECT_DIR'), subdirectory='lib')

coeffs=fltarr(nt)
chi2=fltarr(1)
niter=0L
retval=call_external(soname, 'idl_k_fit_spec', float(coeffs), $
                     float(flux), float(ivar), float(templates), $
                     long(nt), long(nl),float(tolerance), long(maxiter), $
                     long(niter),float(chi2),long(verbose)) 


end
