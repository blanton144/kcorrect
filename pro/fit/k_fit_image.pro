;+
; NAME:
;   k_fit_image
; PURPOSE:
;   fit a lowz sample image pixel-by-pixel 
; COMMENTS:
;   currently experimental
; REVISION HISTORY:
;   06-Jul-2005  MRB, NYU
;-
;------------------------------------------------------------------------------
pro k_fit_image, inlowz

setenv, 'VAGC_REDUX=/global/data/vagc-dr2/vagc0'
name=vagc_name('lowz-atlas',ra=inlowz.ra, dec=inlowz.dec)
jpgname=vagc_name('lowz-full-irg',ra=inlowz.ra, dec=inlowz.dec)
words=strsplit((reverse(strsplit(name,'/',/extr)))[0],'.',/extr)
base=strjoin(words[0:n_elements(words)-2],'.')

filename=base+'.fits'
if(NOT file_test(filename)) then begin
    spawn, 'scp -C mb144@hola.cosmo.fas.nyu.edu:'+jpgname+' '+base+'.jpg'
    spawn, 'scp -C mb144@hola.cosmo.fas.nyu.edu:'+name+' '+filename
endif

udata=mrdfits(filename, 0)
uivar=mrdfits(filename, 1)
gdata=mrdfits(filename, 2)
givar=mrdfits(filename, 3)
rdata=mrdfits(filename, 4)
rivar=mrdfits(filename, 5)
idata=mrdfits(filename, 6)
iivar=mrdfits(filename, 7)
zdata=mrdfits(filename, 8)
zivar=mrdfits(filename, 9)

nx=(size(udata,/dim))[0]
ny=(size(udata,/dim))[1]
ndata=nx*ny

udata=udata[nx/2L-nx/4L:nx/2L+nx/4L, $
            ny/2L-ny/4L:ny/2L+ny/4L]
gdata=gdata[nx/2L-nx/4L:nx/2L+nx/4L, $
            ny/2L-ny/4L:ny/2L+ny/4L]
rdata=rdata[nx/2L-nx/4L:nx/2L+nx/4L, $
            ny/2L-ny/4L:ny/2L+ny/4L]
idata=idata[nx/2L-nx/4L:nx/2L+nx/4L, $
            ny/2L-ny/4L:ny/2L+ny/4L]
zdata=zdata[nx/2L-nx/4L:nx/2L+nx/4L, $
            ny/2L-ny/4L:ny/2L+ny/4L]
uivar=uivar[nx/2L-nx/4L:nx/2L+nx/4L, $
            ny/2L-ny/4L:ny/2L+ny/4L]
givar=givar[nx/2L-nx/4L:nx/2L+nx/4L, $
            ny/2L-ny/4L:ny/2L+ny/4L]
rivar=rivar[nx/2L-nx/4L:nx/2L+nx/4L, $
            ny/2L-ny/4L:ny/2L+ny/4L]
iivar=iivar[nx/2L-nx/4L:nx/2L+nx/4L, $
            ny/2L-ny/4L:ny/2L+ny/4L]
zivar=zivar[nx/2L-nx/4L:nx/2L+nx/4L, $
            ny/2L-ny/4L:ny/2L+ny/4L]

fake_seeing, udata, seeing=1., ampseeing=1.
fake_seeing, gdata, seeing=1., ampseeing=1.
fake_seeing, rdata, seeing=1., ampseeing=1.
fake_seeing, idata, seeing=1., ampseeing=1.
fake_seeing, zdata, seeing=1., ampseeing=1.

nx=(size(udata,/dim))[0]
ny=(size(udata,/dim))[1]
ndata=nx*ny

str=replicate({petroflux:fltarr(5), $
               petroflux_ivar:fltarr(5), $
               ra:inlowz.ra, $
               dec:inlowz.dec}, ndata)
str.petroflux[0]=reform(udata, ndata)
str.petroflux[1]=reform(gdata, ndata)
str.petroflux[2]=reform(rdata, ndata)
str.petroflux[3]=reform(idata, ndata)
str.petroflux[4]=reform(zdata, ndata)
str.petroflux_ivar[0]=reform(uivar, ndata)
str.petroflux_ivar[1]=reform(givar, ndata)
str.petroflux_ivar[2]=reform(rivar, ndata)
str.petroflux_ivar[3]=reform(iivar, ndata)
str.petroflux_ivar[4]=reform(zivar, ndata)

kc=sdss_kcorrect(replicate(inlowz.zdist, ndata), calibobj=str, mtol=mtol, $
                 b300=b300, mass=mass, coeffs=coeffs, absmag=absmag)

mass=reform(10.^(-0.4*absmag[2,*])*mtol[2,*], nx, ny)
m300=mass*reform(b300, nx,ny)

nw_rgb_make, idata, rdata, gdata, scale=mrb_sdss_scales(), $
  name=base+'_image.jpg'

nw_rgb_make, mass, mass, mass, scale=[1.e-4,1.e-4,1.e-4], /invert, $
  name=base+'_mass.jpg'

nw_rgb_make, m300, m300, m300, scale=[3.e-2,3.e-2,3.e-2], /invert, $
  name=base+'_m300.jpg'

end
