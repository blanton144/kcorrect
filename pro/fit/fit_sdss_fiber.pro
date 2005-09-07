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
;   omega0, omegal0 - cosmological parameters for calculating distance
;                     moduli [default 0.3, 0.7]
;   range - [2] range of fit in wavelength
; OPTIONAL KEYWORDS:
;   usevdisp - use correct vdisp for fit 
;   nolines - output model without lines
;   plot - splot the results
; OUTPUTS:
;   coeffs - coefficients fit to each template 
;   model - model spectrum
;   cmodel - model of continuum (no lines, and with large scale
;            differences median smoothed out with 200 Angstrom box)
;   mass, metallicity - properties of template fit;
;                       mass is in units of solar masses and is the
;                       currently remaining stellar mass
;   b300 - star-formation within last 300Myrs relative to average
;          star-formation rate
;   b1000 - star-formation within last 1Gyrs relative to average
;           star-formation rate
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
                    cmodel=cmodel, mass=mass, b300=b300, b1000=b1000, $
                    metallicity=metallicity, nolines=nolines, plot=plot, $
                    omega0=omega0, omegal0=omegal0, range=range, $
                    vdisp=in_vdisp
                    

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

if(keyword_set(in_vdisp)) then $
  vdisp=in_vdisp

k_reconstruct_spec, dum, loglam, /init, vname=vname
readspec, plate, fiberid, mjd=mjd, zans=zans
sdss_spec_block, plate, fiberid, mjd, block_flux=flux, block_ivar=ivar, $
  block_lambda=lambda, avloglam=loglam, /deextinct, vdisp=vdisp
flux=flux*1.e-17
ivar=ivar*1.e+34

if(keyword_set(range)) then begin
    newivar=ivar*0.
    irange=where(lambda gt range[0] and $
                 lambda lt range[1], nrange)
    if(nrange eq 0) then return
    newivar[irange]=ivar[irange]
    ivar=newivar
endif

if(keyword_set(usevdisp)) then vdisp=zans.vdisp
k_fit_spec, flux, ivar, coeffs, vname=vname, vdisp=vdisp, $
  templates=templates
k_reconstruct_spec, coeffs, loglam, vname=vname, $
  metallicity=metallicity, mass=mass, b300=b300, nolines=nolines
if(arg_present(mass)) then $
  mass=mass*10.^(0.4*lf_distmod(zans.z, omega0=omega0, omegal0=omegal0))
if(arg_present(model) gt 0 or keyword_set(plot) gt 0) then $
  k_reconstruct_spec, coeffs, loglam, model, vname=vname, nolines=nolines
if(arg_present(cmodel)) then begin
    k_reconstruct_spec, coeffs, loglam, cmodel, vname=vname, /nolines
    cflux=flux
    ibad=where(ivar le 0., nbad)
    cflux[ibad]=cmodel[ibad]
    dloglam=loglam[1]-loglam[0]
    diff=djs_median(cflux/cmodel, width=long(0.02/dloglam))
    cmodel=cmodel*diff
endif

if(keyword_set(plot)) then begin
    splot, 10.^loglam, flux, xra=[1000.,10000.]
    soplot, 10.^loglam, model, color='red'
endif

end
