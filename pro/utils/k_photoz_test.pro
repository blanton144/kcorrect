;+
; NAME:
;   k_photoz_test
; PURPOSE:
;   test photometric redshifts on SDSS data
; CALLING SEQUENCE:
;   k_photoz_test
; REVISION HISTORY:
;   June 22, 2005 MRB
;-
;------------------------------------------------------------------------------
pro k_photoz_test, vname=vname

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

photoz=sdss_kphotoz(calibobj=cat, vname=vname)


k_print, filename='sdss_resid_main.ps', $
  pold=pold, xold=xold, yold=yold, $
  axis_char_scale=1.1

!X.MARGIN=[0,2]
!Y.MARGIN=[0,0]
!X.OMARGIN=10
!Y.OMARGIN=10
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
!X.OMARGIN=10
!Y.OMARGIN=10
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

end
