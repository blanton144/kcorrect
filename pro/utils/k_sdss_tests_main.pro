;+
; NAME:
;   k_sdss_tests_main
; PURPOSE:
;   runs tests on SDSS test data
; CALLING SEQUENCE:
;   k_sdss_tests_main
; REVISION HISTORY:
;   2005-04-07 MRB, NYU
;-
;------------------------------------------------------------------------------
pro k_sdss_tests_main, vname=vname

im=hogg_mrdfits(vagc_name('object_sdss_imaging'),1,nrow=28800)
sp=hogg_mrdfits(vagc_name('object_sdss_spectro'),1,nrow=28800,columns='z')

ii=where((im.vagc_select and 4) gt 0 and sp.z gt 0.01 and sp.z lt 0.4)
im=im[ii]
sp=sp[ii]
ii=shuffle_indx(n_elements(sp), num_sub=10000)
im=im[ii]
sp=sp[ii]
kc=sdss_kcorrect(sp.z, calibobj=im, band_shift=0.1, rmaggies=rmaggies, $
                 omaggies=omaggies, oivar=oivar, vname=vname, $
                 absmag=absmag, mtol=mtol)
cresid=fltarr(4, n_elements(sp))
for i=0L, 3L do $
  cresid[i,*]=(-2.5*alog10(rmaggies[i,*]/rmaggies[i+1,*]))- $
  (-2.5*alog10(omaggies[i,*]/omaggies[i+1,*]))

k_print, filename='sdss_mtol_main.ps', pold=pold, xold=xold, yold=yold, $
  axis_char_scale=1.1
!P.MULTI=[0,1,1]
hogg_usersym, 10, /fill
djs_plot, absmag[1,*]-absmag[2,*], alog10(mtol[2,*]), psym=8, symsize=0.4, $
  xra=[0.01,1.19]
k_end_print, pold=pold, xold=xold, yold=yold
stop

k_print, filename='sdss_resid_main.ps', $
  pold=pold, xold=xold, yold=yold, $
  axis_char_scale=1.1

!X.MARGIN=[0,2]
!Y.MARGIN=[0,0]
!X.OMARGIN=0
!Y.OMARGIN=0
!P.MULTI=[0,2,2]
ranges=[[-1.09, 1.09], $
        [-0.29, 0.29], $
        [-0.29, 0.29], $
        [-0.39, 0.39]]
ytitle=['\Delta!8(u-g)!6', $
        '\Delta!8(g-r)!6', $
        '\Delta!8(r-i)!6', $
        '\Delta!8(i-z)!6']
xchs=[0.001, 0.001, 1.1, 1.1]
ychs=[1.1, 0.001, 1.1, 0.001]

for i=0, 3 do begin  & $
  hogg_scatterplot, sp.z, cresid[i,*], psym=3, $
  xra=[0.041, 0.179], yra=ranges[*,i], /cond, $
  xnpix=20, ynpix=20, exp=0.5, satfrac=0.001, $
  quantiles=[0.1, 0.25, 0.5, 0.75, 0.9], ytitle=textoidl(ytitle[i]), $
  xtitle='!8z!6', xch=xchs[i], ych=ychs[i] & $
if (i mod 2) eq 1 then $
axis,!X.CRANGE[1],!Y.CRANGE[0],yaxis=1, $
  ytitle=textoidl(ytitle[i]),ycharsize=1.1 & $
endfor

k_end_print, pold=pold, xold=xold, yold=yold

im=hogg_mrdfits(vagc_name('object_sdss_imaging'),1,nrow=28800)
sp=hogg_mrdfits(vagc_name('object_sdss_spectro'),1,nrow=28800,columns='z')

dm=lf_distmod(sp.z)
absm=22.5-2.5*alog10(im.petroflux[2])-im.extinction[2]-dm
mag=22.5-2.5*alog10(im.petroflux[2])-im.extinction[2]
ii=where(sp.z gt 0.05 and sp.z lt 0.17 and $
         absm gt -22.5 and absm lt -20. and $
         (im.vagc_select and 4) gt 0 and mag lt 17.6)
sp=sp[ii]
im=im[ii]
help,im
kc=sdss_kcorrect(sp.z, calibobj=im, band_shift=0.1, vname=vname, absmag=absmag)
ii=where(absmag[2,*] gt -21.5 and absmag[2,*] lt -21.2)
sp=sp[ii]
im=im[ii]
kc=kc[*,ii]
absmag=absmag[*,ii]
dm=lf_distmod(sp.z)

k_print, filename='sdss_colors_main.ps', pold=pold, xold=xold, yold=yold, $
  axis_char_scale=1.1

!X.MARGIN=[0,2]
!Y.MARGIN=[0,0]
!X.OMARGIN=0
!Y.OMARGIN=0
!P.MULTI=[0,2,2]
ranges=[[0.11,2.9], $
        [0.21,1.19], $
        [0.11, 0.55], $
        [-0.05, 0.49]]
ytitle=['!8^{0.1}(u-g)!6', $
        '!8^{0.1}(g-r)!6', $
        '!8^{0.1}(r-i)!6', $
        '!8^{0.1}(i-z)!6']
xchs=[0.001, 0.001, 1.1, 1.1]
ychs=[1.1, 0.001, 1.1, 0.001]

for i=0, 3 do begin  & $
  hogg_scatterplot, sp.z, absmag[i,*]-absmag[i+1,*], psym=3, $
  xra=[0.041, 0.179], yra=ranges[*,i], /cond, $
  xnpix=20, ynpix=20, exp=0.5, satfrac=0.001, $
  quantiles=[0.1, 0.25, 0.5, 0.75, 0.9], ytitle=textoidl(ytitle[i]), $
  xtitle='!8z!6', xch=xchs[i], ych=ychs[i] & $
if (i mod 2) eq 1 then $
axis,!X.CRANGE[1],!Y.CRANGE[0],yaxis=1, $
  ytitle=textoidl(ytitle[i]),ycharsize=1.1 & $
endfor

k_end_print, pold=pold, xold=xold, yold=yold

end
