;+
; NAME:
;   fit_sdss_fiber
; PURPOSE:
;   fit an sdss fiber spectro to a nonnegative sum of templates
; CALLING SEQUENCE:
;   fit_sdss_fiber [, plate, fiberid, mjd=, slist=, model=, coeffs=, $
;          vname=, cmodel, /usevdisp ]
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
;   cmodel - model of continuum (no lines, and with large scale
;            differences median smoothed out with 200 Angstrom box)
;   mass, age, metallicity - properties of template fit;
;                            mass is in units of solar masses; age is in yrs
; EXAMPLE:
;  To get a believable galaxy continuum, use the following, and it
;  will return the galaxy flux in "flux", the weighting in "ivar", and
;  the best fit continuum model in "cmodel."
; 
;  fit_sdss_fiber, 401, 35, mjd=51788, loglam=loglam, flux=flux, $
;     ivar=ivar, cmodel=cmodel, /usevdisp
; 
; REVISION HISTORY:
;   21-Apr-2005  MRB, NYU
;-
;------------------------------------------------------------------------------
pro fit_sdss_fiber, in_plate, in_fiberid, mjd=in_mjd, slist=slist, $
                    model=model, coeffs=coeffs, vname=vname, $
                    usevdisp=usevdisp, flux=flux, ivar=ivar, loglam=loglam, $
                    cmodel=cmodel, age=age, mass=mass, b300=b300, $
                    metallicity=metallicity
                    

if(n_tags(slist) gt 0) then begin
    plate=slist.plate
    fiberid=slist.fiberid
    mjd=slist.mjd
endif else begin
    plate=in_plate
    fiberid=in_fiberid
    if(keyword_set(in_mjd)) then $
      mjd=in_mjd
endelse

k_reconstruct_spec, dum, loglam, /init, vname=vname
readspec, plate, fiberid, mjd=mjd, zans=zans
sdss_spec_block, plate, fiberid, mjd, block_flux=flux, block_ivar=ivar, $
  block_lambda=lambda, avloglam=loglam, /deextinct
flux=flux*1.e-17
ivar=ivar*1.e+34

if(keyword_set(usevdisp)) then vdisp=zans.vdisp
k_fit_spec, flux, ivar, coeffs, vname=vname, vdisp=vdisp, templates=templates
k_reconstruct_spec, coeffs, loglam, vname=vname, age=age, $
  metallicity=metallicity, mass=mass, b300=b300
if(arg_present(mass)) then $
  mass=mass*10.^(0.4*lf_distmod(zans.z))
if(arg_present(model)) then $
  k_reconstruct_spec, coeffs, loglam, model, vname=vname
if(arg_present(cmodel)) then begin
    k_reconstruct_spec, coeffs, loglam, cmodel, vname=vname, /nolines
    cflux=flux
    ibad=where(ivar le 0., nbad)
    cflux[ibad]=cmodel[ibad]
    dloglam=loglam[1]-loglam[0]
    diff=djs_median(cflux/cmodel, width=long(0.02/dloglam))
    cmodel=cmodel*diff
endif

end
