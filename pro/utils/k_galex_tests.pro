;+
; NAME:
;   k_galex_tests
; PURPOSE:
;   runs tests on GALEX+SDSS test data
; CALLING SEQUENCE:
;   k_galex_tests
; REVISION HISTORY:
;   2005-04-07 MRB, NYU
;-
;------------------------------------------------------------------------------
pro k_galex_tests

galex=hogg_mrdfits(getenv('VAGC_REDUX')+'/galex/galex_catalog.fits', 1, $
                   nrow=28800)
objs=hogg_mrdfits(getenv('VAGC_REDUX')+'/galex/galex_objects.fits', 1, $
                  nrow=28800)
ii=where(objs.object_position ge 0)
galex=galex[ii]
objs=objs[ii]
im=hogg_mrdfits(vagc_name('object_sdss_imaging'),1, nrow=28800)
im=im[objs.object_position]
sp=hogg_mrdfits(vagc_name('object_sdss_spectro'),1,columns='z', nrow=28800)
sp=sp[objs.object_position]

ii=where(sp.z gt 0.01 and sp.z lt 0.3)
sp=sp[ii]
im=im[ii]
galex=galex[ii]
objs=objs[ii]

kc1=galex_kcorrect(sp.z, calibobj=im, band_shift=0.1, rmaggies=rmaggies1)
kc0=galex_kcorrect(sp.z, galex=galex, calibobj=im, band_shift=0.1, $
                   rmaggies=rmaggies0, omaggies=omaggies, oivar=oivar)

cresid=fltarr(6, n_elements(sp))
for i=0L, 5L do $
  cresid[i,*]=(-2.5*alog10(rmaggies0[i,*]/rmaggies0[i+1,*]))- $
  (-2.5*alog10(omaggies[i,*]/omaggies[i+1,*]))

k_print, filename='galex_resid.ps', $
  pold=pold, xold=xold, yold=yold, $
  axis_char_scale=2.4

!X.MARGIN=[0,2]
!Y.MARGIN=[0,0]
!X.OMARGIN=0
!Y.OMARGIN=0
!P.MULTI=[0,2,3]
ranges=[[-1.09, 1.09], $
        [-1.09, 1.09], $
        [-1.09, 1.09], $
        [-0.29, 0.29], $
        [-0.29, 0.29], $
        [-0.39, 0.39]]
ytitle=['\Delta!8(F-N)!6', $
        '\Delta!8(N-u)!6', $
        '\Delta!8(u-g)!6', $
        '\Delta!8(g-r)!6', $
        '\Delta!8(r-i)!6', $
        '\Delta!8(i-z)!6']
xchs=[0.001, 0.001, 0.001, 0.001, 2.4, 2.4]
ychs=[2.4, 0.001, 2.4, 0.001, 2.4, 0.001]

for i=0, 5 do begin  & $
  hogg_scatterplot, sp.z, cresid[i,*], psym=3, $
  xra=[0.009, 0.301], yra=ranges[*,i], /cond, $
  xnpix=20, ynpix=20, exp=0.5, satfrac=0.001, $
  quantiles=[0.1, 0.25, 0.5, 0.75, 0.9], ytitle=textoidl(ytitle[i]), $
  xtitle='!8z!6', xch=xchs[i], ych=ychs[i] & $
if (i mod 2) eq 1 then $
axis,!X.CRANGE[1],!Y.CRANGE[0],yaxis=1, $
  ytitle=textoidl(ytitle[i]),ycharsize=2.4 & $
endfor

k_end_print, pold=pold, xold=xold, yold=yold

cresid=fltarr(6, n_elements(sp))
for i=0L, 5L do $
  cresid[i,*]=(-2.5*alog10(rmaggies1[i,*]/rmaggies1[i+1,*]))- $
  (-2.5*alog10(omaggies[i,*]/omaggies[i+1,*]))

k_print, filename='galex_predicted.ps', $
  pold=pold, xold=xold, yold=yold, $
  axis_char_scale=2.4

!X.MARGIN=[0,2]
!Y.MARGIN=[0,0]
!X.OMARGIN=0
!Y.OMARGIN=0
!P.MULTI=[0,2,3]
ranges=[[-1.09, 1.09], $
        [-1.09, 1.09], $
        [-1.09, 1.09], $
        [-0.29, 0.29], $
        [-0.29, 0.29], $
        [-0.39, 0.39]]
ytitle=['\Delta!8(F-N)!6', $
        '\Delta!8(N-u)!6', $
        '\Delta!8(u-g)!6', $
        '\Delta!8(g-r)!6', $
        '\Delta!8(r-i)!6', $
        '\Delta!8(i-z)!6']
xchs=[0.001, 0.001, 0.001, 0.001, 2.4, 2.4]
ychs=[2.4, 0.001, 2.4, 0.001, 2.4, 0.001]

for i=0, 5 do begin  & $
  hogg_scatterplot, sp.z, cresid[i,*], psym=3, $
  xra=[0.009, 0.301], yra=ranges[*,i], /cond, $
  xnpix=20, ynpix=20, exp=0.5, satfrac=0.001, $
  quantiles=[0.1, 0.25, 0.5, 0.75, 0.9], ytitle=textoidl(ytitle[i]), $
  xtitle='!8z!6', xch=xchs[i], ych=ychs[i] & $
if (i mod 2) eq 1 then $
axis,!X.CRANGE[1],!Y.CRANGE[0],yaxis=1, $
  ytitle=textoidl(ytitle[i]),ycharsize=2.4 & $
endfor

k_end_print, pold=pold, xold=xold, yold=yold

dm=lf_distmod(sp.z)
absm=22.5-2.5*alog10(im.petroflux[2])-im.extinction[2]-dm-kc1[4,*]
mag=22.5-2.5*alog10(im.petroflux[2])-im.extinction[2]
ii=where(sp.z gt 0.05 and sp.z lt 0.17 and $
         absm gt -21.5 and absm lt -21.2 and $
         (im.vagc_select and 4) gt 0 and mag lt 17.6)
help,ii
absmag=fltarr(7,n_elements(ii))
for i=0,6 do $
  absmag[i,*]=-2.5*alog10(omaggies[i,ii])-dm[ii]-kc0[i,ii]

k_print, filename='galex_colors_main.ps', pold=pold, xold=xold, yold=yold, $
  axis_char_scale=2.4

!X.MARGIN=[0,2]
!Y.MARGIN=[0,0]
!X.OMARGIN=0
!Y.OMARGIN=0
!P.MULTI=[0,2,3]
ranges=[[-0.49,2.49], $
        [-0.2,3.9], $
        [0.11,2.9], $
        [0.21,1.19], $
        [0.11, 0.55], $
        [-0.05, 0.49]]
ytitle=['!8^{0.1}(F-N)!6', $
        '!8^{0.1}(N-u)!6', $
        '!8^{0.1}(u-g)!6', $
        '!8^{0.1}(g-r)!6', $
        '!8^{0.1}(r-i)!6', $
        '!8^{0.1}(i-z)!6']
xchs=[0.001, 0.001, 0.001, 0.001, 2.4, 2.4]
ychs=[2.4, 0.001, 2.4, 0.001, 2.4, 0.001]

for i=0, 5 do begin  & $
  hogg_scatterplot, sp[ii].z, absmag[i,*]-absmag[i+1,*], psym=3, $
  xra=[0.041, 0.179], yra=ranges[*,i], /cond, $
  xnpix=20, ynpix=20, exp=0.5, satfrac=0.001, $
  quantiles=[0.1, 0.25, 0.5, 0.75, 0.9], ytitle=textoidl(ytitle[i]), $
  xtitle='!8z!6', xch=xchs[i], ych=ychs[i] & $
if (i mod 2) eq 1 then $
axis,!X.CRANGE[1],!Y.CRANGE[0],yaxis=1, $
  ytitle=textoidl(ytitle[i]),ycharsize=2.4 & $
endfor

k_end_print, pold=pold, xold=xold, yold=yold

save
stop

end
