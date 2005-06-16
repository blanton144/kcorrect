;+
; NAME:
;   k_galex_tests
; PURPOSE:
;   runs tests on GALEX+SDSS test data
; CALLING SEQUENCE:
;   k_galex_tests
; DATA DEPENDENCIES:
;   $KCORRECT_DIR/data/test/galex_tests.fits (builds if not there)
; REVISION HISTORY:
;   2005-04-07 MRB, NYU
;-
;------------------------------------------------------------------------------
pro k_galex_tests, vname=vname

galexfile=getenv('KCORRECT_DIR')+'/data/test/galex_tests.fits'
if(NOT file_test(galexfile)) then begin
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

    ii=shuffle_indx(n_elements(galex), num_sub=10000)
    sp=sp[ii]
    im=im[ii]
    galex=galex[ii]
    objs=objs[ii]

    cat1=create_struct(im[0], $
                       sp[0], $
                       struct_trimtags(galex[0], select_tags='*', $
                                       except_tags='id'))
    cat=replicate(cat1, n_elements(im))
    struct_assign, sp, cat
    struct_assign, im, cat, /nozero
    struct_assign, galex, cat, /nozero
    mwrfits, cat, galexfile, /create
endif else begin
    cat=mrdfits(galexfile,1)
endelse

vlimfile=getenv('KCORRECT_DIR')+'/data/test/galex_tests_vlim.fits'
if(NOT file_test(vlimfile)) then begin
    galex=hogg_mrdfits(getenv('VAGC_REDUX')+'/galex/galex_catalog.fits', 1, $
                       nrow=28800)
    objs=hogg_mrdfits(getenv('VAGC_REDUX')+'/galex/galex_objects.fits', 1, $
                      nrow=28800)
    ii=where(objs.object_position ge 0)
    galex=galex[ii]
    objs=objs[ii]
    im=hogg_mrdfits(vagc_name('object_sdss_imaging'),1, nrow=28800)
    im=im[objs.object_position]
    kc=hogg_mrdfits(vagc_name('kcorrect',flux='petro',collision='none', $
                              band_shift='0.10'),1, nrow=28800)
    kc=kc[objs.object_position]
    
    mag=22.5-2.5*alog10(im.petroflux[2])-im.extinction[2]
    ii=where(kc.z gt 0.05 and kc.z lt 0.17 and $
             kc.absmag[2] gt -21.5 and kc.absmag[2] lt -21.2 and $
             (im.vagc_select and 4) gt 0 and mag lt 17.6)
    help,ii
    kc=kc[ii]
    im=im[ii]
    galex=galex[ii]
    objs=objs[ii]
    
    if(n_elements(galex) gt 10000) then begin
        ii=shuffle_indx(n_elements(galex), num_sub=10000)
        kc=kc[ii]
        im=im[ii]
        galex=galex[ii]
        objs=objs[ii]
    endif

    vlim1=create_struct(im[0], $
                        struct_trimtags(kc[0], select_tags='*', $
                                        except_tags=['ra', $
                                                     'dec']), $
                        struct_trimtags(galex[0], select_tags='*', $
                                        except_tags='id'))
    vlim=replicate(vlim1, n_elements(im))
    struct_assign, im, vlim
    struct_assign, kc, vlim, /nozero
    struct_assign, galex, vlim, /nozero
    mwrfits, vlim, vlimfile, /create
endif else begin
    vlim=mrdfits(vlimfile,1)
endelse

filterlist=['sdss_u0.par', $
            'sdss_g0.par', $
            'sdss_r0.par', $
            'sdss_i0.par', $
            'sdss_z0.par']
k_load_vmatrix, vmatrix, lambda, vname=vname
kc=galex_kcorrect(cat.z, calibobj=cat, galex=cat, band_shift=0.0, $
                  rmaggies=rmaggies, omaggies=omaggies, oivar=oivar, $
                  vname=vname, absmag=absmag, mtol=mtol, coeffs=coeffs)
nt=(size(coeffs,/dim))[0]
tmtol=fltarr(nt)
k_reconstruct_maggies, identity(nt), replicate(0.001,nt), rmaggies0, $
  vmatrix=vmatrix, lambda=lambda, filterlist=filterlist
for j=0, nt-1L do $
  tmtol[j]=1./rmaggies0[2,j]* $
  10.^(-0.4*k_solar_magnitudes(filterlist=filterlist[2]))
tgmr=-2.5*alog10(rmaggies0[1,*]/rmaggies0[2,*])
k_print, filename='galex_mtol.ps', pold=pold, xold=xold, yold=yold, $
  axis_char_scale=1.1
!P.MULTI=[0,1,1]
hogg_usersym, 10, /fill
djs_plot, absmag[1,*]-absmag[2,*], alog10(mtol[2,*]), psym=8, symsize=0.4, $
  xra=[-0.31,1.79],yra=[-1.59, 1.29], xtitle='!8g-r!6', ytitle='(M/L)(r)'
djs_oplot, tgmr, alog10(tmtol), psym=8, symsize=2., color='red'
k_end_print, pold=pold, xold=xold, yold=yold

kc1=galex_kcorrect(cat.z, calibobj=cat, band_shift=0.1, rmaggies=rmaggies1, $
                   vname=vname)
kc0=galex_kcorrect(cat.z, galex=cat, calibobj=cat, band_shift=0.1, $
                   rmaggies=rmaggies0, omaggies=omaggies, oivar=oivar, $
                   vname=vname)

cresid=fltarr(6, n_elements(cat))
for i=0L, 5L do $
  cresid[i,*]=(-2.5*alog10(rmaggies0[i,*]/rmaggies0[i+1,*]))- $
  (-2.5*alog10(omaggies[i,*]/omaggies[i+1,*]))

k_print, filename='galex_resid.ps', $
  pold=pold, xold=xold, yold=yold, $
  axis_char_scale=2.4

!X.MARGIN=[0,2]
!Y.MARGIN=[0,0]
!X.OMARGIN=10.
!Y.OMARGIN=10.
!P.MULTI=[0,2,3]
ranges=[[-1.99, 1.99], $
        [-1.99, 1.99], $
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

for i=0, 5 do begin  
    ii=where(abs(cresid[i,*]) lt 1.e+6)
    hogg_scatterplot, cat[ii].z, cresid[i,ii], psym=3, $
      xra=[0.009, 0.301], yra=ranges[*,i], /cond, $
      xnpix=20, ynpix=20, exp=0.5, satfrac=0.001, $
      quantiles=[0.1, 0.25, 0.5, 0.75, 0.9], ytitle=textoidl(ytitle[i]), $
      xtitle='!8z!6', xch=xchs[i], ych=ychs[i] 
    if (i mod 2) eq 1 then $
      axis,!X.CRANGE[1],!Y.CRANGE[0],yaxis=1, $
      ytitle=textoidl(ytitle[i]),ycharsize=2.4 
endfor

k_end_print, pold=pold, xold=xold, yold=yold

cresid=fltarr(6, n_elements(cat))
for i=0L, 5L do $
  cresid[i,*]=(-2.5*alog10(rmaggies1[i,*]/rmaggies1[i+1,*]))- $
  (-2.5*alog10(omaggies[i,*]/omaggies[i+1,*]))

k_print, filename='galex_predicted.ps', $
  pold=pold, xold=xold, yold=yold, $
  axis_char_scale=2.4

!X.MARGIN=[0,2]
!Y.MARGIN=[0,0]
!X.OMARGIN=10.
!Y.OMARGIN=10.
!P.MULTI=[0,2,3]
ranges=[[-1.99, 1.99], $
        [-1.99, 1.99], $
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

for i=0, 5 do begin  
    ii=where(abs(cresid[i,*]) lt 1.e+6)
    hogg_scatterplot, cat[ii].z, cresid[i,ii], psym=3, $
      xra=[0.009, 0.301], yra=ranges[*,i], /cond, $
      xnpix=20, ynpix=20, exp=0.5, satfrac=0.001, $
      quantiles=[0.1, 0.25, 0.5, 0.75, 0.9], ytitle=textoidl(ytitle[i]), $
      xtitle='!8z!6', xch=xchs[i], ych=ychs[i] 
    if (i mod 2) eq 1 then $
      axis,!X.CRANGE[1],!Y.CRANGE[0],yaxis=1, $
      ytitle=textoidl(ytitle[i]),ycharsize=2.4 
endfor

k_end_print, pold=pold, xold=xold, yold=yold

kc=galex_kcorrect(vlim.z, galex=vlim, calibobj=vlim, band_shift=0.1, $
                  rmaggies=rmaggies0, omaggies=omaggies, oivar=oivar, $
                  vname=vname, absmag=absmag)

k_print, filename='galex_colors_main.ps', pold=pold, xold=xold, yold=yold, $
  axis_char_scale=2.4

!X.MARGIN=[0,2]
!Y.MARGIN=[0,0]
!X.OMARGIN=10.
!Y.OMARGIN=10.
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

for i=0, 5 do begin  
    hogg_scatterplot, vlim.z, absmag[i,*]-absmag[i+1,*], psym=3, $
      xra=[0.041, 0.179], yra=ranges[*,i], /cond, $
      xnpix=15, ynpix=15, exp=0.5, satfrac=0.001, $
      quantiles=[0.1, 0.25, 0.5, 0.75, 0.9], ytitle=textoidl(ytitle[i]), $
      xtitle='!8z!6', xch=xchs[i], ych=ychs[i] 
    if (i mod 2) eq 1 then $
      axis,!X.CRANGE[1],!Y.CRANGE[0],yaxis=1, $
      ytitle=textoidl(ytitle[i]),ycharsize=2.
endfor

k_end_print, pold=pold, xold=xold, yold=yold


end
