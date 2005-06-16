;+
; NAME:
;   k_twomass_tests
; PURPOSE:
;   runs tests on 2MASS+SDSS test data
; CALLING SEQUENCE:
;   k_twomass_tests
; DATA DEPENDENCIES:
;   $KCORRECT_DIR/data/test/twomass_tests.fits (builds if not there)
; REVISION HISTORY:
;   2005-04-07 MRB, NYU
;-
;------------------------------------------------------------------------------
pro k_twomass_tests, vname=vname

twomassfile=getenv('KCORRECT_DIR')+'/data/test/twomass_tests.fits'
if(NOT file_test(twomassfile)) then begin
    twomass=hogg_mrdfits(getenv('VAGC_REDUX')+'/object_twomass.fits',1, $
                         nrow=28800)
    im=hogg_mrdfits(vagc_name('object_sdss_imaging'),1, nrow=28800)
    sp=hogg_mrdfits(vagc_name('object_sdss_spectro'),1,columns='z', nrow=28800)
    
    ii=where(twomass.twomass_tag ge 0 AND sp.z gt 0.01 AND sp.z lt 0.3, nii)
    sp=sp[ii]
    im=im[ii]
    twomass=twomass[ii]

    if(nii gt 10000) then begin
        ii=shuffle_indx(n_elements(twomass), num_sub=10000)
        sp=sp[ii]
        im=im[ii]
        twomass=twomass[ii]
    endif

    cat1=create_struct(im[0], sp[0], $
                       struct_trimtags(twomass[0], select_tags='*', $
                                       except_tags='ra'))
    cat=replicate(cat1, n_elements(im))
    struct_assign, sp, cat
    struct_assign, im, cat, /nozero
    struct_assign, twomass, cat, /nozero
    mwrfits, cat, twomassfile, /create
endif else begin
    cat=mrdfits(twomassfile,1)
endelse
    
kc1=twomass_kcorrect(cat.z, calibobj=cat, band_shift=0.1, rmaggies=rmaggies1, $
                     vname=vname)
kc0=twomass_kcorrect(cat.z, twomass=cat, calibobj=cat, band_shift=0.1, $
                     rmaggies=rmaggies0, omaggies=omaggies, oivar=oivar, $
                     vname=vname)

cresid=fltarr(7, n_elements(cat))
for i=0L, 6L do $
  cresid[i,*]=(-2.5*alog10(rmaggies0[i,*]/rmaggies0[i+1,*]))- $
  (-2.5*alog10(omaggies[i,*]/omaggies[i+1,*]))

k_print, filename='twomass_resid.ps', $
  pold=pold, xold=xold, yold=yold, $
  axis_char_scale=2.4

!X.MARGIN=[0,2]
!Y.MARGIN=[0,0]
!X.OMARGIN=10
!Y.OMARGIN=10
!P.MULTI=[0,2,4]
ranges=[[-1.09, 1.09], $
        [-0.29, 0.29], $
        [-0.29, 0.29], $
        [-0.39, 0.39], $
        [-1.09, 1.09], $
        [-1.09, 1.09], $
        [-1.09, 1.09]]
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
  hogg_scatterplot, cat.z, cresid[i,*], psym=3, $
  xra=[0.009, 0.301], yra=ranges[*,i], /cond, $
  xnpix=20, ynpix=20, exp=0.5, satfrac=0.001, $
  quantiles=[0.1, 0.25, 0.5, 0.75, 0.9], ytitle=textoidl(ytitle[i]), $
  xtitle='!8z!6', xch=xchs[i], ych=ychs[i] & $
if (i mod 2) eq 1 then $
axis,!X.CRANGE[1],!Y.CRANGE[0],yaxis=1, $
  ytitle=textoidl(ytitle[i]),ycharsize=2.4 & $
endfor

k_end_print, pold=pold, xold=xold, yold=yold

cresid=fltarr(7, n_elements(cat))
for i=0L, 6L do $
  cresid[i,*]=(-2.5*alog10(rmaggies1[i,*]/rmaggies1[i+1,*]))- $
  (-2.5*alog10(omaggies[i,*]/omaggies[i+1,*]))

k_print, filename='twomass_predicted.ps', $
  pold=pold, xold=xold, yold=yold, $
  axis_char_scale=2.4

!X.MARGIN=[0,2]
!Y.MARGIN=[0,0]
!X.OMARGIN=10
!Y.OMARGIN=10
!P.MULTI=[0,2,4]
ranges=[[-1.09, 1.09], $
        [-0.29, 0.29], $
        [-0.29, 0.29], $
        [-0.39, 0.39], $
        [-1.09, 1.09], $
        [-1.09, 1.09], $
        [-1.09, 1.09]]
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
  hogg_scatterplot, cat.z, cresid[i,*], psym=3, $
  xra=[0.009, 0.301], yra=ranges[*,i], /cond, $
  xnpix=20, ynpix=20, exp=0.5, satfrac=0.001, $
  quantiles=[0.1, 0.25, 0.5, 0.75, 0.9], ytitle=textoidl(ytitle[i]), $
  xtitle='!8z!6', xch=xchs[i], ych=ychs[i] & $
if (i mod 2) eq 1 then $
axis,!X.CRANGE[1],!Y.CRANGE[0],yaxis=1, $
  ytitle=textoidl(ytitle[i]),ycharsize=2.4 & $
endfor

k_end_print, pold=pold, xold=xold, yold=yold

dm=lf_distmod(cat.z)
absm=22.5-2.5*alog10(cat.petroflux[2])-cat.extinction[2]-dm-kc1[2,*]
mag=22.5-2.5*alog10(cat.petroflux[2])-cat.extinction[2]
ii=where(cat.z gt 0.05 and cat.z lt 0.17 and $
         absm gt -21.5 and absm lt -21.2 and $
         (cat.vagc_select and 4) gt 0 and mag lt 17.6)
help,ii
absmag=fltarr(8,n_elements(ii))
for i=0,7 do $
  absmag[i,*]=-2.5*alog10(omaggies[i,ii])-dm[ii]-kc0[i,ii]

k_print, filename='twomass_colors_main.ps', pold=pold, xold=xold, yold=yold, $
  axis_char_scale=2.4

!X.MARGIN=[0,2]
!Y.MARGIN=[0,0]
!X.OMARGIN=10
!Y.OMARGIN=10
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
  hogg_scatterplot, cat[ii].z, absmag[i,*]-absmag[i+1,*], psym=3, $
  xra=[0.041, 0.179], yra=ranges[*,i], /cond, $
  xnpix=20, ynpix=20, exp=0.5, satfrac=0.001, $
  quantiles=[0.1, 0.25, 0.5, 0.75, 0.9], ytitle=textoidl(ytitle[i]), $
  xtitle='!8z!6', xch=xchs[i], ych=ychs[i] & $
if (i mod 2) eq 1 then $
axis,!X.CRANGE[1],!Y.CRANGE[0],yaxis=1, $
  ytitle=textoidl(ytitle[i]),ycharsize=2.4 & $
endfor

k_end_print, pold=pold, xold=xold, yold=yold

end
