;+
; NAME:
;   fit_sdss_fiber
; PURPOSE:
;   fit an sdss fiber spectro to a nonnegative sum of templates
; CALLING SEQUENCE:
;   fit_sdss_fiber [, plate, fiberid, mjd=, slist=, model=, coeffs=, $
;          vname=, /usevdisp ]
; INPUTS:
;   plate, fiberid, mjd - SDSS id numbers 
;     OR
;   slist - structure containing above id numbers in .PLATE, .FIBERID,
;           and .MJD
; OPTIONAL INPUTS:
;   vname - name of fit to use (default 'default')
; OPTIONAL KEYWORDS:
;   usevdisp - use correct vdisp for fit 
; OUTPUTS:
;   coeffs - coefficients fit to each template 
;   model - model spectrum
; REVISION HISTORY:
;   21-Apr-2005  MRB, NYU
;-
;------------------------------------------------------------------------------
pro fit_sdss_fiber, in_plate, in_fiberid, mjd=in_mjd, slist=slist, $
                    model=model, coeffs=coeffs, vname=vname, $
                    usevdisp=usevdisp, flux=flux, ivar=ivar, loglam=loglam

if(n_tags(slist) gt 0) then begin
    plate=slist.plate
    fiberid=slist.fiberid
    mjd=slist.mjd
endif else begin
    plate=in_plate
    fiberid=in_fiberid
    mjd=in_mjd
endelse

k_reconstruct_spec, dum, loglam, /init, vname=vname
readspec, plate, fiberid, mjd=mjd, zans=zans
sdss_spec_block, plate, fiberid, mjd, block_flux=flux, block_ivar=ivar, $
  block_lambda=lambda, avloglam=loglam ;;, /deextinct

if(keyword_set(usevdisp)) then vdisp=zans.vdisp
k_fit_spec, flux, ivar, coeffs, vname=vname, vdisp=vdisp, templates=templates
if(arg_present(model)) then model=templates#coeffs

end
