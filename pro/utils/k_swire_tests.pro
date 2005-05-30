;+
; NAME:
;   k_swire_tests
; PURPOSE:
;   runs tests on SWIRE data
; CALLING SEQUENCE:
;   k_swire_tests
; DATA DEPENDENCIES:
;   $KCORRECT_DIR/data/test/swire_tests.fits (builds if not there)
; REVISION HISTORY:
;   2005-04-07 MRB, NYU
;-
;------------------------------------------------------------------------------
pro k_swire_tests, vname=vname

swirefile=getenv('KCORRECT_DIR')+'/data/test/swire_tests.fits'
if(NOT file_test(swirefile)) then begin
    swire=mrdfits(getenv('VAGC_REDUX')+'/spitzer/swire_catalog.fits',1)
    objects=mrdfits(getenv('VAGC_REDUX')+'/spitzer/swire_objects.fits',1)

    iin=where(objects.object_position ge 0)
    objects=objects[iin]
    swire=swire[iin]

    sp=mrdfits(vagc_name('object_sdss_spectro'),1, $
               row=objects.object_position)
    im=mrdfits(vagc_name('object_sdss_imaging'),1, $
               row=objects.object_position)

    igal=where(sp.z gt 0.003 and sp.z lt 0.5)
    sp=sp[igal]
    im=im[igal]
    swire=swire[igal]
    
    cat1=create_struct(im[0], sp[0], $
                       struct_trimtags(swire[0], select_tags='*', $
                                       except_tags=['ra', $
                                                    'tile', $
                                                    'dec']))
    cat=replicate(cat1, n_elements(im))
    struct_assign, sp, cat
    struct_assign, im, cat, /nozero
    struct_assign, swire, cat, /nozero
    mwrfits, cat, swirefile, /create
endif else begin
    cat=mrdfits(swirefile,1)
endelse

swire_to_maggies, cat, swire_maggies, swire_ivar
sdss_to_maggies, sdss_maggies, sdss_ivar, calibobj=cat
swire_maggies[0:4,*]=sdss_maggies
swire_ivar[0:4,*]=sdss_ivar

filterlist=['sdss_u0.par', $
            'sdss_g0.par', $
            'sdss_r0.par', $
            'sdss_i0.par', $
            'sdss_z0.par', $
            'spitzer_irac_ch1.par', $
            'spitzer_irac_ch2.par', $
            'spitzer_irac_ch3.par', $
            'spitzer_irac_ch4.par', $
            'spitzer_mips_24.par']
kcorrect, swire_maggies, swire_ivar, cat.z, kc, absmag=absmag, vname=vname, $
  filterlist=filterlist, rmaggies=rmaggies

cresid=fltarr(9, n_elements(cat))
for i=0L, 8L do $
  cresid[i,*]=(-2.5*alog10(rmaggies[i,*]/rmaggies[i+1,*]))- $
  (-2.5*alog10(swire_maggies[i,*]/swire_maggies[i+1,*]))

k_print, filename='swire_resid.ps', $
  pold=pold, xold=xold, yold=yold, $
  axis_char_scale=2.4

!X.MARGIN=[0,2]
!Y.MARGIN=[0,0]
!X.OMARGIN=0
!Y.OMARGIN=0
!P.MULTI=[0,3,3]
ranges=[[-1.09, 1.09], $
        [-1.09, 1.09], $
        [-1.09, 1.09], $
        [-1.09, 1.09], $
        [-1.09, 1.09], $
        [-1.09, 1.09], $
        [-1.09, 1.09], $
        [-1.09, 1.09], $
        [-1.09, 1.09]]
ytitle=['\Delta!8(u-g)!6', $
        '\Delta!8(g-r)!6', $
        '\Delta!8(r-i)!6', $
        '\Delta!8(i-z)!6', $
        '\Delta!8(z-[3.6])!6', $
        '\Delta!8([3.6]-[4.5])!6', $
        '\Delta!8([4.5]-[5.8])!6', $
        '\Delta!8([5.8]-[8.0])!6', $
        '\Delta!8([5.8]-[8.0])!6'  $
       ]
xchs=[0.001, 0.001, 0.001, 0.001, 2.4, 2.4]
ychs=[2.4, 0.001, 2.4, 0.001, 2.4, 0.001]

hogg_usersym, 10, /fill
for i=0, 8 do begin  & $
  djs_plot, cat.z, cresid[i,*], psym=8, symsize=0.4, $
  xra=[0.009, 0.301], yra=ranges[*,i], $
  quantiles=[0.1, 0.25, 0.5, 0.75, 0.9], ytitle=textoidl(ytitle[i]), $
  xtitle='!8z!6' ;;, xch=xchs[i], ych=ychs[i] & $
if (i mod 2) eq 1 then $
  axis,!X.CRANGE[1],!Y.CRANGE[0],yaxis=1, $
  ytitle=textoidl(ytitle[i]),ycharsize=2.4 & $
  endfor

k_end_print, pold=pold, xold=xold, yold=yold

filterlist=['sdss_u0.par', $
            'sdss_g0.par', $
            'sdss_r0.par', $
            'sdss_i0.par', $
            'sdss_z0.par', $
            'spitzer_irac_ch1.par', $
            'spitzer_irac_ch2.par', $
            'spitzer_irac_ch3.par', $
            'spitzer_irac_ch4.par', $
            'spitzer_mips_24.par']
swire_ivar[5:9,*]=0.
kcorrect, swire_maggies, swire_ivar, cat.z, kc, absmag=absmag, vname=vname, $
  filterlist=filterlist, rmaggies=rmaggies

cresid=fltarr(9, n_elements(cat))
for i=0L, 8L do $
  cresid[i,*]=(-2.5*alog10(rmaggies[i,*]/rmaggies[i+1,*]))- $
  (-2.5*alog10(swire_maggies[i,*]/swire_maggies[i+1,*]))

k_print, filename='swire_predicted.ps', $
  pold=pold, xold=xold, yold=yold, $
  axis_char_scale=2.4

!X.MARGIN=[0,2]
!Y.MARGIN=[0,0]
!X.OMARGIN=0
!Y.OMARGIN=0
!P.MULTI=[0,3,3]
ranges=[[-1.09, 1.09], $
        [-1.09, 1.09], $
        [-1.09, 1.09], $
        [-1.09, 1.09], $
        [-1.09, 1.09], $
        [-1.09, 1.09], $
        [-1.09, 1.09], $
        [-1.09, 1.09], $
        [-1.09, 1.09]]
ytitle=['\Delta!8(u-g)!6', $
        '\Delta!8(g-r)!6', $
        '\Delta!8(r-i)!6', $
        '\Delta!8(i-z)!6', $
        '\Delta!8(z-[3.6])!6', $
        '\Delta!8([3.6]-[4.5])!6', $
        '\Delta!8([4.5]-[5.8])!6', $
        '\Delta!8([5.8]-[8.0])!6', $
        '\Delta!8([5.8]-[8.0])!6'  $
       ]
xchs=[0.001, 0.001, 0.001, 0.001, 2.4, 2.4]
ychs=[2.4, 0.001, 2.4, 0.001, 2.4, 0.001]

hogg_usersym, 10, /fill
for i=0, 8 do begin  & $
  djs_plot, cat.z, cresid[i,*], psym=8, symsize=0.4, $
  xra=[0.009, 0.301], yra=ranges[*,i], $
  quantiles=[0.1, 0.25, 0.5, 0.75, 0.9], ytitle=textoidl(ytitle[i]), $
  xtitle='!8z!6' ;;, xch=xchs[i], ych=ychs[i] & $
if (i mod 2) eq 1 then $
  axis,!X.CRANGE[1],!Y.CRANGE[0],yaxis=1, $
  ytitle=textoidl(ytitle[i]),ycharsize=2.4 & $
  endfor

k_end_print, pold=pold, xold=xold, yold=yold

end
