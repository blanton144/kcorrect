;+
; NAME:
;   k_twomass_tests
; PURPOSE:
;   runs tests on 2MASS+SDSS test data
; CALLING SEQUENCE:
;   k_twomass_tests
; REVISION HISTORY:
;   2005-04-07 MRB, NYU
;-
;------------------------------------------------------------------------------
pro k_twomass_tests, vname=vname

twomass=hogg_mrdfits(getenv('VAGC_REDUX')+'/object_twomass.fits',1, $
                     nrow=28800)
im=hogg_mrdfits(vagc_name('object_sdss_imaging'),1, nrow=28800)
sp=hogg_mrdfits(vagc_name('object_sdss_spectro'),1,columns='z', nrow=28800)

ii=where(twomass.twomass_tag ge 0 AND sp.z gt 0.01 AND sp.z lt 0.3)
sp=sp[ii]
im=im[ii]
twomass=twomass[ii]

kc1=twomass_kcorrect(sp.z, calibobj=im, band_shift=0.1, rmaggies=rmaggies1, $
                     vname=vname)
kc0=twomass_kcorrect(sp.z, twomass=twomass, calibobj=im, band_shift=0.1, $
                     rmaggies=rmaggies0, omaggies=omaggies, oivar=oivar, $
                     vname=vname)

cresid=fltarr(7, n_elements(sp))
for i=0L, 6L do $
  cresid[i,*]=(-2.5*alog10(rmaggies0[i,*]/rmaggies0[i+1,*]))- $
  (-2.5*alog10(omaggies[i,*]/omaggies[i+1,*]))

k_print, filename='twomass_resid.ps', $
  pold=pold, xold=xold, yold=yold, $
  axis_char_scale=2.4

!X.MARGIN=[0,2]
!Y.MARGIN=[0,0]
!X.OMARGIN=0
!Y.OMARGIN=0
!P.MULTI=[0,2,4]
ranges=[[-1.09, 1.09], $
        [-0.29, 0.29], $
        [-0.29, 0.29], $
        [-0.39, 0.39], $
        [-0.39, 0.39], $
        [-0.39, 0.39], $
        [-0.39, 0.39]]
ytitle=['\Delta!8(u-g)!6', $
        '\Delta!8(g-r)!6', $
        '\Delta!8(r-i)!6', $
        '\Delta!8(i-z)!6', $
        '\Delta!8(z-J)!6', $
        '\Delta!8(J-H)!6', $
        '\Delta!8(H-K_s)!6']
xchs=[0.001, 0.001, 0.001, 0.001, 0.001, 2.4, 2.4]
ychs=[2.4, 0.001, 2.4, 0.001, 2.4, 0.001, 2.4, 0.001]

for i=0, 6 do begin  & $
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

cresid=fltarr(7, n_elements(sp))
for i=0L, 6L do $
  cresid[i,*]=(-2.5*alog10(rmaggies1[i,*]/rmaggies1[i+1,*]))- $
  (-2.5*alog10(omaggies[i,*]/omaggies[i+1,*]))

k_print, filename='twomass_predicted.ps', $
  pold=pold, xold=xold, yold=yold, $
  axis_char_scale=2.4

!X.MARGIN=[0,2]
!Y.MARGIN=[0,0]
!X.OMARGIN=0
!Y.OMARGIN=0
!P.MULTI=[0,2,3]
ranges=[[-1.09, 1.09], $
        [-0.29, 0.29], $
        [-0.29, 0.29], $
        [-0.39, 0.39], $
        [-0.39, 0.39], $
        [-0.39, 0.39], $
        [-0.39, 0.39]]
ytitle=['\Delta!8(u-g)!6', $
        '\Delta!8(g-r)!6', $
        '\Delta!8(r-i)!6', $
        '\Delta!8(i-z)!6', $
        '\Delta!8(z-J)!6', $
        '\Delta!8(J-H)!6', $
        '\Delta!8(H-K_s)!6']
xchs=[0.001, 0.001, 0.001, 0.001, 0.001, 2.4, 2.4]
ychs=[2.4, 0.001, 2.4, 0.001, 2.4, 0.001, 2.4, 0.001]

for i=0, 6 do begin  & $
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
absm=22.5-2.5*alog10(im.petroflux[2])-im.extinction[2]-dm-kc1[2,*]
mag=22.5-2.5*alog10(im.petroflux[2])-im.extinction[2]
ii=where(sp.z gt 0.05 and sp.z lt 0.17 and $
         absm gt -21.5 and absm lt -21.2 and $
         (im.vagc_select and 4) gt 0 and mag lt 17.6)
help,ii
absmag=fltarr(8,n_elements(ii))
for i=0,7 do $
  absmag[i,*]=-2.5*alog10(omaggies[i,ii])-dm[ii]-kc0[i,ii]

k_print, filename='galex_colors_main.ps', pold=pold, xold=xold, yold=yold, $
  axis_char_scale=2.4

!X.MARGIN=[0,2]
!Y.MARGIN=[0,0]
!X.OMARGIN=0
!Y.OMARGIN=0
!P.MULTI=[0,2,4]
ranges=[[0.11,2.9], $
        [0.21,1.19], $
        [0.11, 0.55], $
        [-0.05, 0.49], $
        [-1., 1.], $
        [-1., 1.], $
        [-1., 1.]]
ytitle=['!8^{0.1}(u-g)!6', $
        '!8^{0.1}(g-r)!6', $
        '!8^{0.1}(r-i)!6', $
        '!8^{0.1}(i-z)!6', $
        '!8^{0.1}(z-J)!6', $
        '!8^{0.1}(J-H)!6', $
        '!8^{0.1}(H-K_s)!6']
xchs=[0.001, 0.001, 0.001, 0.001, 0.001, 2.4, 2.4]
ychs=[2.4, 0.001, 2.4, 0.001, 2.4, 0.001, 2.4, 0.001]

for i=0, 6 do begin  & $
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
