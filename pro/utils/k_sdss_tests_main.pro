;+
; NAME:
;   k_sdss_tests_main
; PURPOSE:
;   runs tests on SDSS test data
; CALLING SEQUENCE:
;   k_sdss_tests_main
; DATA DEPENDENCIES:
;   $KCORRECT_DIR/data/test/sdss_tests_main.fits (builds if not there)
; REVISION HISTORY:
;   2005-04-07 MRB, NYU
;-
;------------------------------------------------------------------------------
pro k_sdss_tests_main, vname=vname

testfile=getenv('KCORRECT_DIR')+'/data/test/sdss_tests_main.fits'
if(NOT file_test(testfile)) then begin
    im=hogg_mrdfits(vagc_name('object_sdss_imaging'),1,nrow=28800)
    sp=hogg_mrdfits(vagc_name('object_sdss_spectro'),1,nrow=28800,columns='z')

    ii=where((im.vagc_select and 4) gt 0 and sp.z gt 0.01 and sp.z lt 0.4)
    im=im[ii]
    sp=sp[ii]
    ii=shuffle_indx(n_elements(sp), num_sub=20000)
    im=im[ii]
    sp=sp[ii]

    cat1=create_struct(im[0], sp[0])
    cat=replicate(cat1, n_elements(im))
    struct_assign, sp, cat
    struct_assign, im, cat, /nozero
    mwrfits, cat, testfile, /create
endif else begin
    cat=mrdfits(testfile,1)
endelse

vlimfile=getenv('KCORRECT_DIR')+'/data/test/sdss_tests_main_vlim.fits'
if(NOT file_test(vlimfile)) then begin
    im=hogg_mrdfits(vagc_name('object_sdss_imaging'),1,nrow=28800)
    sp=hogg_mrdfits(vagc_name('object_sdss_spectro'),1,nrow=28800,columns='z')

    kc=sdss_kcorrect(sp.z, calibobj=im, band_shift=0.1, vname=vname, $
                     absmag=absmag)
    ii=where(absmag[2,*] gt -21.5 and absmag[2,*] lt -21.2 and $
             sp.z gt 0.03 and sp.z lt 0.19)
    kc=kc[*,ii]
    absmag=absmag[*,ii]
    sp=sp[ii]
    im=im[ii]

    vlim1=create_struct(im[0], sp[0])
    vlim=replicate(vlim1, n_elements(im))
    struct_assign, sp, vlim
    struct_assign, im, vlim, /nozero
    mwrfits, vlim, vlimfile, /create
endif else begin
    vlim=mrdfits(vlimfile,1)
endelse

kcb=sdss2bessell(vlim.z, calibobj=vlim, band_shift=0.0, rmaggies=rmaggies, $
                 omaggies=omaggies, oivar=oivar, vname=vname, $
                 absmag=absmag, mtol=mtol, /vega, coeffs=coeffs)
k_print, filename='sdss_mtol_main.ps', pold=pold, xold=xold, yold=yold, $
  axis_char_scale=1.1
!P.MULTI=[0,1,1]
hogg_usersym, 10, /fill
djs_plot, absmag[1,*]-absmag[3,*], alog10(mtol[3,*]), psym=8, symsize=0.4, $
  xra=[0.01,1.79],yra=[-0.59, 0.89], xtitle='!8B-R!6', ytitle='(M/L)(R)'
djs_oplot, [0.4,1.6], -0.82+0.85*[0.4,1.6], th=4, color='red'
djs_oplot, [0.4,1.6], 0.4-0.82+0.85*[0.4,1.6], th=4, color='red'
k_end_print, pold=pold, xold=xold, yold=yold

kc=sdss_kcorrect(cat.z, calibobj=cat, band_shift=0.1, rmaggies=rmaggies, $
                 omaggies=omaggies, oivar=oivar, vname=vname, $
                 absmag=absmag, mtol=mtol)
cresid=fltarr(4, n_elements(cat))
for i=0L, 3L do $
  cresid[i,*]=(-2.5*alog10(rmaggies[i,*]/rmaggies[i+1,*]))- $
  (-2.5*alog10(omaggies[i,*]/omaggies[i+1,*]))

k_print, filename='sdss_resid_main.ps', $
  pold=pold, xold=xold, yold=yold, $
  axis_char_scale=1.1

!X.MARGIN=[0,2]
!Y.MARGIN=[0,0]
!X.OMARGIN=0
!Y.OMARGIN=0
!P.MULTI=[0,2,2]
ranges=[[-0.49, 0.49], $
        [-0.19, 0.19], $
        [-0.19, 0.19], $
        [-0.19, 0.19]]
ytitle=['\Delta!8(u-g)!6', $
        '\Delta!8(g-r)!6', $
        '\Delta!8(r-i)!6', $
        '\Delta!8(i-z)!6']
xchs=[0.001, 0.001, 1.1, 1.1]
ychs=[1.1, 0.001, 1.1, 0.001]

for i=0, 3 do begin  & $
  hogg_scatterplot, cat.z, cresid[i,*], psym=3, $
  xra=[0.01, 0.31], yra=ranges[*,i], /cond, $
  xnpix=20, ynpix=20, exp=0.5, satfrac=0.001, $
  quantiles=[0.1, 0.25, 0.5, 0.75, 0.9], ytitle=textoidl(ytitle[i]), $
  xtitle='!8z!6', xch=xchs[i], ych=ychs[i] & $
if (i mod 2) eq 1 then $
axis,!X.CRANGE[1],!Y.CRANGE[0],yaxis=1, $
  ytitle=textoidl(ytitle[i]),ycharsize=1.1 & $
endfor

k_end_print, pold=pold, xold=xold, yold=yold

kc=sdss_kcorrect(vlim.z, calibobj=vlim, band_shift=0.1, vname=vname, $
                 absmag=absmag)

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
  hogg_scatterplot, vlim.z, absmag[i,*]-absmag[i+1,*], psym=3, $
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
