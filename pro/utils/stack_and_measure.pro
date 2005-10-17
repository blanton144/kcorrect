;+
; NAME:
;   stack_and_measure
; PURPOSE:
;   Stacks all SDSS imaging, and measures the object supposedly there
; CALLING SEQUENCE:
;   stack_and_measure, ra, dec [, sz=, rerun=, /fpbin]
; INPUTS:
;   ra, dec - [N] J2000 positions
;   rerun - rerun
; OPTIONAL INPUTS:
;   sz - size of stack in deg (defaults to 0.008)
; OPTIONAL KEYWORDS:
;   /fpbin - stack reconstructed frames
; COMMENTS:
;   Stacks all imaging at each RA, DEC separately
; REVISION HISTORY:
;   17-Oct-2005  Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro stack_and_measure, ra, dec, sz=sz, rerun=rerun, fpbin=fpbin, fluxes=fluxes

seed=-1L
if(NOT keyword_set(sz)) then sz=0.008

fluxes=fltarr(5, n_elements(ra))
for i=0L, n_elements(ra)-1L do begin
    name=hogg_iau_name(ra[i], dec[i], '')
    prefix=name[0]
    if(file_test(prefix+'-u.fits.gz') eq 0 OR $
       file_test(prefix+'-g.fits.gz') eq 0 OR $
       file_test(prefix+'-r.fits.gz') eq 0 OR $
       file_test(prefix+'-i.fits.gz') eq 0 OR $
       file_test(prefix+'-z.fits.gz') eq 0) then $
      smosaic_make, ra[i], dec[i], sz, sz, rerun=rerun, fpbin=fpbin, $
      prefix=prefix, pixscale=0.396/3600., seed=seed, /all

    uimage=mrdfits(prefix+'-u.fits.gz')
    gimage=mrdfits(prefix+'-g.fits.gz')
    rimage=mrdfits(prefix+'-r.fits.gz')
    iimage=mrdfits(prefix+'-i.fits.gz')
    zimage=mrdfits(prefix+'-z.fits.gz')

    nx=(size(rimage,/dim))[0]
    ny=(size(rimage,/dim))[1]

    weight=psf_gaussian(npix=nx, fwhm=2.5, centroid=[nx/2., ny/2.]-1., /norm, $
                        ndim=2)
    fluxes[0,i]=total(uimage*weight)/total(weight^2)
    fluxes[1,i]=total(gimage*weight)/total(weight^2)
    fluxes[2,i]=total(rimage*weight)/total(weight^2)
    fluxes[3,i]=total(iimage*weight)/total(weight^2)
    fluxes[4,i]=total(zimage*weight)/total(weight^2)
    
endfor

end
;------------------------------------------------------------------------------

