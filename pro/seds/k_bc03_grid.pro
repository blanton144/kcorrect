;+
; NAME:
;   k_bc03_grid
; PURPOSE:
;   make grid of BC03 models
; CALLING SEQUENCE:
;   k_bc03_grid
; REVISION HISTORY:
;   29-Jul-2004  Michael Blanton (NYU)
;-
;------------------------------------------------------------------------------
pro k_bc03_grid

minage=1.e-3
maxage=14.
nages=40L
ages= 10.^(alog10(minage)+(alog10(maxage)-alog10(minage))*findgen(nages)/ $
           float(nages-1L))

nmets=6
mets=[0,1,2,3,4,5]

tmp_bc03= im_read_bc03(age=1.)
nl=n_elements(tmp_bc03.flux)

wave=tmp_bc03.wave
loglam=alog10(wave)
grid=fltarr(nl, nages, nmets)

for im= 0L, nmets-1L do $
  grid[*,*,im]= (im_read_bc03(age=ages,met=mets[im])).flux


navloglam=4000L
avloglam=double(alog10(3500.)+(alog10(9500.)-alog10(3500.))* $
                (dindgen(navloglam)+0.5)/float(navloglam))
sfgrid=fltarr(navloglam, nages*nmets)

ninterloglam=20000L
interloglam=double(alog10(3500.)+(alog10(9500.)-alog10(3500.))* $
                   (dindgen(ninterloglam)+0.5)/float(ninterloglam))

for im= 0L, nmets-1L do $
  for ia= 0L, nages-1L do begin & $
  splog, string(im)+string(ia) & $
  intergrid=interpol(grid[*,ia,im], loglam, interloglam) & $
  combine1fiber, interloglam, intergrid, fltarr(ninterloglam)+1., $
  newloglam=avloglam, newflux=tmp1, maxiter=0 & $
  sfgrid[*,ia+im*nages]=tmp1 & $
  endfor

ndust=3L
dust1={dusty_str, geometry:'', dust:'', structure:'', tauv:0.}
dust=replicate(dust1,3)
dust.geometry=['dusty', $
               'dusty','dusty']
dust.dust=['MW', $
           'MW','MW']
dust.structure=['c', $
                'c','c']
dust.tauv=[0.,1.,3.]

dustygrid=fltarr(navloglam, nages*nmets, ndust)
dustfact=fltarr(navloglam, ndust)
for i=0L, ndust-1L do $
  dustfact[*,i]=exp(-witt_ext(dust[i],dust[i].tauv,10.^(avloglam)))
for i=0L, ndust-1L do $
  for j=0L, nages*nmets-1L do $
  dustygrid[*,j,i]=sfgrid[*,j]*dustfact[*,i]

outgrid=reform(dustygrid,navloglam,nages*nmets*ndust)
tauv=fltarr(nages,nmets,ndust)
for i=0L, ndust-1L do tauv[*,*,i]=dust[i].tauv
met=fltarr(nages,nmets,ndust)
for i=0L, nmets-1L do met[*,i,*]=mets[i]
age=fltarr(nages,nmets,ndust)
for i=0L, nages-1L do age[i,*,*]=ages[i]

mwrfits, outgrid, 'mmatrix.fits', /create
mwrfits, tauv, 'mmatrix.fits'
mwrfits, met, 'mmatrix.fits'
mwrfits, age, 'mmatrix.fits'


end
