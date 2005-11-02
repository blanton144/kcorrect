;+
; NAME:
;   measure_deep_stack
; PURPOSE:
;   read in DEEP2 images, measure and match them to DEEP2
; REVISION HISTORY:
;   17-Oct-2005  Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro measure_deep_stack

files=file_search('J*-z.fits.gz')

fluxes=fltarr(5,n_elements(files))
ra=dblarr(n_elements(files))
dec=dblarr(n_elements(files))
for i=0L, n_elements(files)-1L do begin

    prefix=(stregex(files[i], '(.*)-z.fits.gz', /sub, /extr))[1]
    
    uimage=mrdfits(prefix+'-u.fits.gz',0,hdr)
    ra[i]=double(sxpar(hdr,'CRVAL1'))
    dec[i]=double(sxpar(hdr,'CRVAL2'))
    gimage=mrdfits(prefix+'-g.fits.gz')
    rimage=mrdfits(prefix+'-r.fits.gz')
    iimage=mrdfits(prefix+'-i.fits.gz')
    zimage=mrdfits(prefix+'-z.fits.gz')
    
    nx=(size(rimage,/dim))[0]
    ny=(size(rimage,/dim))[1]
    
    weight=psf_gaussian(npix=nx, fwhm=4.5, centroid=[nx/2., ny/2.]-1., /norm, $
                        ndim=2)
    fluxes[0,i]=total(uimage*weight)/total(weight^2)
    fluxes[1,i]=total(gimage*weight)/total(weight^2)
    fluxes[2,i]=total(rimage*weight)/total(weight^2)
    fluxes[3,i]=total(iimage*weight)/total(weight^2)
    fluxes[4,i]=total(zimage*weight)/total(weight^2)
endfor

deep=mrdfits(getenv('KCORRECT_DIR')+ $
             '/data/redshifts/deep/zcat.dr1.uniq.fits.gz',1)
spherematch, ra, dec, deep.ra, deep.dec, 3./3600., m1, m2, d12

stack=replicate(create_struct(deep[0], $
                              'modelflux', fltarr(5), $
                              'modelflux_ivar', fltarr(5), $
                              'sdss_bri', fltarr(3), $
                              'sdss_absmag', fltarr(5), $
                              'fake_deep_absmag', fltarr(3), $
                              'deep_absmag', fltarr(3), $
                              'fake_deep_absmag_ext', fltarr(3), $
                              'deep_absmag_ext', fltarr(3), $
                              'bessell_absmag', fltarr(5)), $
                n_elements(m1))
struct_assign, deep[m2], stack
stack.modelflux=fluxes[*,m1]
stack.modelflux_ivar=1./((stack.modelflux*0.1>0.01))^2
stack.sdss_bri=sdss2deep(stack.z, stack.z, calibobj=stack, flux='model')
kc=deep_kcorrect(stack.z, zcat=stack, absmag=deep_absmag, $
                 filterlist=['bessell_B.par', $
                             'bessell_R.par', $
                             'bessell_I.par'])
stack.deep_absmag=deep_absmag
kc=deep_kcorrect(stack.z, mag=stack.sdss_bri, $
                 err=replicate(0.1, 3, n_elements(stack)), $
                 absmag=fake_deep_absmag, $
                 filterlist=['bessell_B.par', $
                             'bessell_R.par', $
                             'bessell_I.par'])
stack.fake_deep_absmag=fake_deep_absmag
kc=deep_kcorrect(stack.z, zcat=stack, absmag=deep_absmag_ext)
stack.deep_absmag_ext=deep_absmag_ext
kc=deep_kcorrect(stack.z, mag=stack.sdss_bri, $
                 err=replicate(0.1, 3, n_elements(stack)), $
                 absmag=fake_deep_absmag_ext)
stack.fake_deep_absmag_ext=fake_deep_absmag_ext
kc=sdss2bessell(stack.z, calibobj=stack, absmag=bessell_absmag, flux='model')
stack.bessell_absmag=bessell_absmag
kc=sdss_kcorrect(stack.z, calibobj=stack, absmag=sdss_absmag, flux='model')
stack.sdss_absmag=sdss_absmag
                     
mwrfits, stack, 'deep_stack.fits', /create

ii=where(stack.magr lt 22)
hogg_usersym, 10, /fill
k_print, filename='deep_stack.ps'
!P.MULTi=[0,1,2]
deep_bmr=stack.magb-stack.magr
sdss_bmr=stack.sdss_bri[0]-stack.sdss_bri[1]
deep_rmi=stack.magr-stack.magi
sdss_rmi=stack.sdss_bri[1]-stack.sdss_bri[2]
djs_plot, stack[ii].z, sdss_rmi[ii]-deep_rmi[ii], psym=8, symsize=0.5, $
  xra=[0.05, 0.29], yra=[-1.,1.], xti='z', yti='\Delta (R-I)', xchar=0.0001
djs_plot, stack[ii].z, sdss_bmr[ii]-deep_bmr[ii], psym=8, symsize=0.5, $
  xra=[0.05, 0.29], yra=[-1.,1.], xti='z', yti='\Delta (B-R)'
deep_bmr0=stack.deep_absmag[0]-stack.deep_absmag[1]
fake_deep_bmr0=stack.fake_deep_absmag[0]-stack.fake_deep_absmag[1]
sdss_bmr0=stack.bessell_absmag[1]-stack.bessell_absmag[3]
deep_rmi0=stack.deep_absmag[1]-stack.deep_absmag[2]
fake_deep_rmi0=stack.fake_deep_absmag[1]-stack.fake_deep_absmag[2]
sdss_rmi0=stack.bessell_absmag[3]-stack.bessell_absmag[4]
djs_plot, stack[ii].z, sdss_rmi0[ii]-deep_rmi0[ii], psym=8, symsize=0.5, $
  xra=[0.05, 0.29], yra=[-1.,1.], xti='z', yti='\Delta (R-I)_0', xchar=0.0001
djs_plot, stack[ii].z, sdss_bmr0[ii]-deep_bmr0[ii], psym=8, symsize=0.5, $
  xra=[0.05, 0.29], yra=[-1.,1.], xti='z', yti='\Delta (B-R)_0'
k_end_print

k_print, filename=getenv('EVOLUTION_DIR')+'/tex/deep2_lowz_comp.ps'
!P.MULTI=[0,1,2]
deep_umb0=stack.deep_absmag_ext[1]-stack.deep_absmag_ext[2]
fake_deep_umb0_ext=stack.fake_deep_absmag_ext[1]-stack.fake_deep_absmag_ext[2]
fake_deep_umb0=stack.bessell_absmag[0]-stack.bessell_absmag[1]
djs_plot, stack[ii].z, deep_umb0[ii]-fake_deep_umb0_ext[ii], psym=8, $
  symsize=0.5, $
  xra=[0.05, 0.29], yra=[-1.,1.], xti='!6redshift !8z!6', $
  yti='\Delta!8(U-B)_0!6', xch=0.001
djs_xyouts, 0.07, -0.8, charsize=1.3, $
  '!8Converting to DEEP2-like data!6'
djs_plot, stack[ii].z, deep_umb0[ii]-fake_deep_umb0[ii], psym=8, symsize=0.5, $
  xra=[0.05, 0.29], yra=[-1.,1.], xti='z', yti='\Delta!8(U-B)_0!6'
djs_xyouts, 0.07, -0.8, charsize=1.3, $
  '!8Converting SDSS data directly!6'
k_end_print 
save

end
;------------------------------------------------------------------------------

