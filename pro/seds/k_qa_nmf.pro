;+
; NAME:
;   k_qa_nmf
; REVISION HISTORY:
;   30-Nov-2004  Michael Blanton (NYU)
;-
;------------------------------------------------------------------------------
pro k_qa_nmf

mmatrix=mrdfits('k_nmf_mmatrix.fits',0,hdr)
lambda=mrdfits('k_nmf_mmatrix.fits',1)
zf=mrdfits('k_nmf_mmatrix.fits',6)
nspec=long(sxpar(hdr, 'NSPEC'))
nzf=long(sxpar(hdr, 'NZ'))
nfilter=long(sxpar(hdr, 'NFILTER'))
ndusts=long(sxpar(hdr, 'NDUST'))
nmets=long(sxpar(hdr, 'NMET'))
nages=long(sxpar(hdr, 'NAGE'))

data=mrdfits('k_nmf_data.fits',1)
ivar=mrdfits('k_nmf_data.fits',2)
zhelio=mrdfits('k_nmf_data.fits',5)

ww=mrdfits('k_nmf_soln.fits')
hh=mrdfits('k_nmf_soln.fits',1)
model=mmatrix#ww#hh

iz=long(floor((nzf-1.)*(zhelio-zf[0])/(zf[nzf-1]-zf[0])+0.5))

mfluxes=fltarr(10, n_elements(zhelio))
fluxes=fltarr(10, n_elements(zhelio))
fluxes_ivar=fltarr(10, n_elements(zhelio))
for ifilter=0L, 9L do $
  mfluxes[ifilter,*]= model[iz+ifilter*nzf+nspec, lindgen(n_elements(zhelio))]
for ifilter=0L, 9L do $
  fluxes[ifilter,*]= data[iz+ifilter*nzf+nspec, lindgen(n_elements(zhelio))]
for ifilter=0L, 9L do $
  fluxes_ivar[ifilter,*]= ivar[iz+ifilter*nzf+nspec, lindgen(n_elements(zhelio))]

igood=where(fluxes[0,*] gt 0. and $
            fluxes[1,*] gt 0.)
fmn=-2.5*alog10(fluxes[0,igood]/fluxes[1,igood])
nmu=-2.5*alog10(fluxes[1,igood]/fluxes[2,igood])
umg=-2.5*alog10(fluxes[2,igood]/fluxes[3,igood])
gmr=-2.5*alog10(fluxes[3,igood]/fluxes[4,igood])
rmi=-2.5*alog10(fluxes[4,igood]/fluxes[5,igood])
imz=-2.5*alog10(fluxes[5,igood]/fluxes[6,igood])
mfmn=-2.5*alog10(mfluxes[0,igood]/mfluxes[1,igood])
mnmu=-2.5*alog10(mfluxes[1,igood]/mfluxes[2,igood])
mumg=-2.5*alog10(mfluxes[2,igood]/mfluxes[3,igood])
mgmr=-2.5*alog10(mfluxes[3,igood]/mfluxes[4,igood])
mrmi=-2.5*alog10(mfluxes[4,igood]/mfluxes[5,igood])
mimz=-2.5*alog10(mfluxes[5,igood]/mfluxes[6,igood])

splot, fmn, gmr, psym=4
soplot, mfmn, mgmr, psym=4, color='red'
splot, fmn, nmu, psym=4
soplot, mfmn, mnmu, psym=4, color='red'
splot, umg, gmr, psym=4
soplot, mumg, mgmr, psym=4, color='red'
splot, gmr, rmi, psym=4
soplot, mgmr, mrmi, psym=4, color='red'
splot, rmi, imz, psym=4
soplot, mrmi, mimz, psym=4, color='red'
;;splot, gmr, rmi, psym=4
;;soplot, mgmr, mrmi, psym=4, color='red'

end
