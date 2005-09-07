;+
; NAME:
;   k_goods_tests
; PURPOSE:
;   runs tests on GOODS test data
; CALLING SEQUENCE:
;   k_goods_tests
; DATA DEPENDENCIES:
;   $KCORRECT_DIR/data/test/goods_tests.fits (builds if not there)
; REVISION HISTORY:
;   2005-04-07 MRB, NYU
;-
;------------------------------------------------------------------------------
pro k_goods_tests, vname=vname

savfile='data_k_goods_tests.sav'

if(NOT file_test(savfile)) then begin
    gphoto=rsex(getenv('KCORRECT_DIR')+ $
                '/data/redshifts/goods/mb_cdfs_isaac_ks_photz_c1.1_d3.0k.cat')
    gz=rsex(getenv('KCORRECT_DIR')+ $
            '/data/redshifts/goods/mb_z.cdfs.c1.1z.20050218.cat')
    igot=where(gz.z lt 2. and gz.z gt 0.1)
    gz=gz[igot]
    gphoto=gphoto[igot]
    goods_to_maggies, gphoto, goods_maggies, goods_ivar

    lowz=lowz_read(sample='dr4')
    ii=where(lowz.absmag[2] lt -16)
    lowz=lowz[ii]
    im=hogg_mrdfits(vagc_name('object_sdss_imaging'),1, nrow=28800)
    im=im[lowz.object_position]
    
    for i=0L, 0L do begin
        maxvmax=max(1./lowz.vmax)
        chances=(1./lowz.vmax)/maxvmax
        tmp_igoods=where(randomu(seed, n_elements(lowz)) lt chances, npick)
        tmp_goods_redshift=0.5+1.*randomu(seed, npick)
        tmp_mag=sdss2goods(lowz[tmp_igoods].zdist, tmp_goods_redshift, $
                              calibobj=im[tmp_igoods], vname=vname)
        ;; CHANGE FOR GOODS
        ipass=where(tmp_mag[4,*] lt 25.1, npass)
        if(npass gt 0) then begin
            tmp_igoods=tmp_igoods[ipass]
            tmp_goods_redshift=tmp_goods_redshift[ipass]
            if(n_elements(igoods) eq 0) then begin
                igoods=tmp_igoods
                goods_redshift=tmp_goods_redshift
            endif else begin
                igoods=[igoods, tmp_igoods]
                goods_redshift=[goods_redshift, tmp_goods_redshift]
            endelse
        endif
        help, igoods
    endfor
    sdss_redshift=lowz[igoods].zdist
    im=im[igoods]
    mag=sdss2goods(sdss_redshift, goods_redshift, calibobj=im, vname=vname)
    
    ;; CHANGE FOR GOODS
    ipass=where(mag[4,*] lt 25.1, npass)
    mag=mag[*,ipass]
    sdss_redshift=sdss_redshift[ipass]
    goods_redshift=goods_redshift[ipass]
    im=im[ipass]

    save, filename=savfile 
endif else begin
    restore, savfile
endelse 

hogg_usersym, 8, /fill

k_print, filename='k_test_sdss2goods.ps'
hogg_usersym, 10, /fill
!P.MULTI=[0,2,3]
!Y.MARGIN=0

ytitle='!8B-V!6'
djs_plot, goods_redshift, mag[0,*]-mag[1,*], psym=8, symsize=0.3, $
  xra=[0.5, 1.5], xcharsize=0.0001, ycharsize=1.4, $
  ytitle=ytitle, yra=[-0.3,4.5], quantiles=[0.05, 0.2, 0.5, 0.8, 0.95]
djs_oplot, gz.z, gphoto.bmag_magauto-gphoto.vmag_magauto, psym=8, $
  symsize=0.35, color='red'

ytitle='!8V-i!6'
djs_plot, goods_redshift, mag[1,*]-mag[2,*], psym=8, symsize=0.3, $
  xra=[0.5, 1.5], xcharsize=0.0001, ycharsize=0.0001, $
  ytitle=ytitle, yra=[-0.3,2.3], quantiles=[0.05, 0.2, 0.5, 0.8, 0.95]
djs_oplot, gz.z, gphoto.vmag_magauto-gphoto.imag_magauto, psym=8, $
  symsize=0.35, color='red'
axis,!X.CRANGE[1],!Y.CRANGE[0],yaxis=1, $
  ytitle=textoidl(ytitle),ycharsize=1.1 

ytitle='!8i-z!6'
djs_plot, goods_redshift, mag[2,*]-mag[3,*], psym=8, symsize=0.3, $
  xra=[0.5, 1.5], xcharsize=0.0001, ycharsize=1.4, $
  ytitle=ytitle, yra=[-0.3,1.5], quantiles=[0.05, 0.2, 0.5, 0.8, 0.95]
djs_oplot, gz.z, gphoto.imag_magauto-gphoto.zmag_magauto, psym=8, $
  symsize=0.35, color='red'

ytitle='!8z-J!6'
djs_plot, goods_redshift, mag[3,*]-mag[4,*], psym=8, symsize=0.3, $
  xra=[0.5, 1.5], xcharsize=0.0001, ycharsize=0.0001, $
  ytitle=ytitle, yra=[-0.3,1.9], quantiles=[0.05, 0.2, 0.5, 0.8, 0.95]
djs_oplot, gz.z, gphoto.zmag_magauto-gphoto.jmag_magauto, psym=8, $
  symsize=0.35, color='red'
axis,!X.CRANGE[1],!Y.CRANGE[0],yaxis=1, $
  ytitle=textoidl(ytitle),ycharsize=1.1 

ytitle='!8J-H!6'
djs_plot, goods_redshift, mag[4,*]-mag[5,*], psym=8, symsize=0.3, $
  xra=[0.5, 1.5], xcharsize=1.4001, ycharsize=1.4, $
  ytitle=ytitle, yra=[-0.3,1.1], quantiles=[0.05, 0.2, 0.5, 0.8, 0.95], $
  xtitle='z'
djs_oplot, gz.z, gphoto.jmag_magauto-gphoto.hmag_magauto, psym=8, $
  symsize=0.35, color='red'

ytitle='!8H-K!6'
djs_plot, goods_redshift, mag[5,*]-mag[6,*], psym=8, symsize=0.3, $
  xra=[0.5, 1.5], xcharsize=1.4001, ycharsize=0.0001, $
  ytitle=ytitle, yra=[-0.3,0.9], quantiles=[0.05, 0.2, 0.5, 0.8, 0.95], $
  xtitle='z'
djs_oplot, gz.z, gphoto.hmag_magauto-gphoto.kmag_magauto, psym=8, $
  symsize=0.35, color='red'
axis,!X.CRANGE[1],!Y.CRANGE[0],yaxis=1, $
  ytitle=textoidl(ytitle),ycharsize=1.1 

k_end_print

kc=goods_kcorrect(gz.z, goods=gphoto,absmag=absmag, vname=vname, $
                  rmaggies=rmaggies,omaggies=omaggies)
resid=fltarr(7, n_elements(gz))
for i=0L, 6L do $
  resid[i,*]=(-2.5*alog10(rmaggies[i,*]))- $
  (-2.5*alog10(omaggies[i,*]))

k_print, filename='k_test_goods_resid.ps', yold=yold, xold=xold, pold=pold
!X.MARGIN=[0,2]
!Y.MARGIN=[0,0]
!X.OMARGIN=0
!Y.OMARGIN=0
!P.MULTI=[0,2,4]
ranges=[[-0.49, 0.49], $
        [-0.49, 0.49], $
        [-0.49, 0.49], $
        [-0.49, 0.49], $
        [-0.49, 0.49], $
        [-0.49, 0.49], $
        [-0.49, 0.49]]
ytitle=['\Delta!8B!6', $
        '\Delta!8V!6', $
        '\Delta!8i!6', $
        '\Delta!8z!6', $
        '\Delta!8J!6', $
        '\Delta!8H!6', $
        '\Delta!8K!ds!n!6']
xchs=[0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 1.1, 1.1]
ychs=[1.1, 0.001, 1.1, 0.001, 1.1, 0.001, 1.1, 0.001]

for i=0, 6 do begin  
    ii=where(abs(resid[i,*]) lt 1.e+6)
    hogg_scatterplot, gz[ii].z, resid[i,ii], psym=3, $
      xra=[0.01, 2.1], yra=ranges[*,i], /cond, $
      xnpix=10, ynpix=15, exp=0.5, satfrac=0.001, $
      quantiles=[0.1, 0.25, 0.5, 0.75, 0.9], ytitle=textoidl(ytitle[i]), $
      xtitle='!8z!6', xch=xchs[i], ych=ychs[i] 
    if (i mod 2) eq 1 then $
      axis,!X.CRANGE[1],!Y.CRANGE[0],yaxis=1, $
      ytitle=textoidl(ytitle[i]),ycharsize=1.1 
endfor
k_end_print, pold=pold, xold=xold, yold=yold

end
