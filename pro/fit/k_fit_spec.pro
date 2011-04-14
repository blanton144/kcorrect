;+
; NAME:
;   k_fit_spec
; PURPOSE:
;   fit a spectrum to the sum of templates 
; CALLING SEQUENCE:
;   k_fit_spec, flux, ivar, coeffs [, vname=, vdisp=, templates=, $
;     lambda=, oflux=, oivar=, olambda= ]
; INPUTS:
;   flux - [nl] fluxes (aligned with wavelength grid of models)
;   ivar - [nl] inverse variances (aligned with wavelength grid of models)
; OPTIONAL INPUT/OUTPUTS:
;   templates - [nl] fluxes of templates to fit
;   lambda - [nl] wavelenths 
; OPTIONAL INPUTS:
;   oflux - [nm] aligned version of flux
;   oivar - [nm] aligned version of ivar
;   olambda - [nm] model wavelength grid
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
                templates=templates, lambda=lambda, oflux=oflux, $
                oivar=oivar, olambda=olambda, nolines=nolines, $
                linear=linear, chi2=chi2, tolerance=tolerance, $
                maxiter=maxiter, verbose=verbose

if(NOT keyword_set(ivar)) then begin
    inz=where(flux ne 0., nnz)
    minflux=min(abs(flux[inz]))
    ivar=1./((abs(flux)>minflux)*0.1)^2
endif

if(n_elements(olambda) eq 0) then begin
    k_reconstruct_spec, dum, loglam, /init, vname=vname, nt=nt, nolines=nolines
    olambda=10.^loglam
endif else begin
    loglam= alog10(olambda)
    nl=n_elements(loglam)
    nt=n_elements(templates)/n_elements(olambda)
endelse

if(keyword_set(lambda)) then begin
    nf=n_elements(flux)
    dl=lambda[1]-lambda[0]
    oflux=interpol([flux[0], flux[0], flux, flux[nf-1], flux[nf-1]], $
                   alog10([lambda[0]-2.*dl, lambda[0], lambda, $
                           lambda[nf-1], lambda[nf-1]+2.*dl]), loglam)
    oivar=interpol([0., 0., ivar, 0., 0.], $
                   alog10([lambda[0]-2.*dl, lambda[0], lambda, $
                           lambda[nf-1], lambda[nf-1]+2.*dl]), loglam)
endif else begin
    oflux=flux
    oivar=ivar
endelse

if(keyword_set(templates) eq 0 OR $
   keyword_set(vname) gt 0) then begin
    nl=n_elements(loglam)
    templates=fltarr(nl, nt)
    for i=0L, nt-1L do begin
        coeffs=fltarr(nt)
        coeffs[i]=1.
        k_reconstruct_spec, coeffs, loglam, tmpflux, vname=vname, $
          vdisp=vdisp, nolines=nolines
        templates[*,i]=tmpflux
    endfor
endif

if(n_elements(verbose) eq 0) then verbose=1L
if(NOT keyword_set(maxiter)) then maxiter=50000L
if(NOT keyword_set(tolerance)) then tolerance=1.e-13

; Set source object name
soname=filepath('libkcorrect.'+kcorrect_so_ext(), $
                root_dir=getenv('KCORRECT_DIR'), subdirectory='lib')

coeffs=fltarr(nt)+1.
chi2=fltarr(1)
niter=0L

if(NOT keyword_set(linear)) then begin
    retval=call_external(soname, 'idl_k_fit_spec', float(coeffs), $
                         float(oflux), float(oivar), float(templates), $
                         long(nt), long(nl),float(tolerance), long(maxiter), $
                         long(niter),float(chi2),long(verbose)) 
endif else begin
    retval=call_external(soname, 'idl_k_fit_spec_linear', float(coeffs), $
                         float(oflux), float(oivar), float(templates), $
                         long(nt), long(nl), float(chi2),long(verbose)) 
endelse


end
