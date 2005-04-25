;+
; NAME:
;   k_sdss_tests_lrg
; PURPOSE:
;   runs tests on SDSS test data
; CALLING SEQUENCE:
;   k_sdss_tests_main
; REVISION HISTORY:
;   2005-04-07 MRB, NYU
;-
;------------------------------------------------------------------------------
pro k_sdss_tests_lrg, vname=vname

im=hogg_mrdfits(vagc_name('object_sdss_imaging'),1,nrow=28800)
ti=hogg_mrdfits(vagc_name('object_sdss_tiling'),1,nrow=28800)
sp=hogg_mrdfits(vagc_name('object_sdss_spectro'),1,nrow=28800,columns='z')

ii=where((ti.primtarget and 32) gt 0)
im=im[ii]
sp=sp[ii]
help,im
kc=sdss_kcorrect(sp.z, calibobj=im, band_shift=0.3, rmaggies=rmaggies, $
                 omaggies=omaggies, oivar=oivar, vname=vname)
cresid=fltarr(4, n_elements(sp))
for i=0L, 3L do $
  cresid[i,*]=(-2.5*alog10(rmaggies[i,*]/rmaggies[i+1,*]))- $
  (-2.5*alog10(omaggies[i,*]/omaggies[i+1,*]))

k_print, filename='sdss_resid_lrg.ps', $
  pold=pold, xold=xold, yold=yold, $
  axis_char_scale=1.1

!X.MARGIN=[0,2]
!Y.MARGIN=[0,0]
!X.OMARGIN=0
!Y.OMARGIN=0
!P.MULTI=[0,2,2]
ranges=[[-1.49, 1.49], $
        [-0.49, 0.49], $
        [-0.49, 0.49], $
        [-0.49, 0.49]]
ytitle=['\Delta!8(u-g)!6', $
        '\Delta!8(g-r)!6', $
        '\Delta!8(r-i)!6', $
        '\Delta!8(i-z)!6']
xchs=[0.001, 0.001, 1.1, 1.1]
ychs=[1.1, 0.001, 1.1, 0.001]

for i=0, 3 do begin  & $
  hogg_scatterplot, sp.z, cresid[i,*], psym=3, $
  xra=[0.12, 0.56], yra=ranges[*,i], /cond, $
  xnpix=20, ynpix=20, exp=0.5, satfrac=0.001, $
  quantiles=[0.1, 0.25, 0.5, 0.75, 0.9], ytitle=textoidl(ytitle[i]), $
  xtitle='!8z!6', xch=xchs[i], ych=ychs[i] & $
if (i mod 2) eq 1 then $
axis,!X.CRANGE[1],!Y.CRANGE[0],yaxis=1, $
  ytitle=textoidl(ytitle[i]),ycharsize=1.1 & $
endfor

k_end_print, pold=pold, xold=xold, yold=yold

dm=lf_distmod(sp.z)
absmag=fltarr(5,n_elements(sp))
for i=0, 4 do $
  absmag[i,*]=22.5-2.5*alog10(im.petroflux[i])-im.extinction[i]-kc[i,*]-dm

k_print, filename='sdss_colors_lrg.ps', pold=pold, xold=xold, yold=yold, $
  axis_char_scale=1.1

!X.MARGIN=[0,2]
!Y.MARGIN=[0,0]
!X.OMARGIN=0
!Y.OMARGIN=0
!P.MULTI=[0,2,2]
ranges=[[-0.41,4.9], $
        [0.51,2.49], $
        [0.11, 0.95], $
        [-0.25, 0.85]]
ytitle=['!8^{0.3}(u-g)!6', $
        '!8^{0.3}(g-r)!6', $
        '!8^{0.3}(r-i)!6', $
        '!8^{0.3}(i-z)!6']
xchs=[0.001, 0.001, 1.1, 1.1]
ychs=[1.1, 0.001, 1.1, 0.001]

for i=0, 3 do begin  & $
  hogg_scatterplot, sp.z, absmag[i,*]-absmag[i+1,*], psym=3, $
  xra=[0.12, 0.56], yra=ranges[*,i], /cond, $
  xnpix=20, ynpix=20, exp=0.5, satfrac=0.001, $
  quantiles=[0.1, 0.25, 0.5, 0.75, 0.9], ytitle=textoidl(ytitle[i]), $
  xtitle='!8z!6', xch=xchs[i], ych=ychs[i] & $
if (i mod 2) eq 1 then $
axis,!X.CRANGE[1],!Y.CRANGE[0],yaxis=1, $
  ytitle=textoidl(ytitle[i]),ycharsize=1.1 & $
endfor

k_end_print, pold=pold, xold=xold, yold=yold

end
