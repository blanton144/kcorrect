;+
; NAME:
;   paper_plots
; PURPOSE:
;   create plots for kcorrect.tex paper
; CALLING SEQUENCE:
;   paper_plots
; REVISION HISTORY:
;   15-Aug-2003  Written by M. Blanton, NYU
;-  
;------------------------------------------------------------------------------
pro paper_plots

gphoto=rsex(getenv('KCORRECT_DIR')+ $
            '/data/redshifts/goods/mb_cdfs_isaac_ks_photz_c1.1_d3.0k.cat')
gz=rsex(getenv('KCORRECT_DIR')+ $
        '/data/redshifts/goods/mb_z.cdfs.c1.1z.20050218.cat')
igot=where(gz.z lt 2. and gz.z gt 0.1)
gz=gz[igot]
gphoto=gphoto[igot]
kc=goods_kcorrect(gz.z, goods=gphoto, $
                  omaggies=maggies, rmaggies=rmaggies, band_shift=0., $
                  oivar=ivar, vname='goods')
k_print, filename=getenv('KCORRECT_DIR')+'/docs/paper/goods_special.ps', $
  xold=xold, yold=yold, pold=pold
!P.MULTI=[6,2,3]
!X.MARGIN=9.
!Y.MARGIN=3.
colors=fltarr(9, n_elements(gphoto))
mcolors=fltarr(9, n_elements(gphoto))
ioff=[0,1,3,4,5,6]
for i=0L, 5L do $
  colors[i,*]=-2.5*alog10(maggies[ioff[i],*]/maggies[2,*])
for i=0L, 5L do $
  mcolors[i,*]=-2.5*alog10(rmaggies[ioff[i],*]/rmaggies[2,*])
ytitle=['\Delta!8(B-i)!6', $
        '\Delta!8(V-i)!6', $
        '\Delta!8(z-i)!6', $
        '\Delta!8(J-i)!6', $
        '\Delta!8(H-i)!6', $
        '\Delta!8(K_s-i)!6']
hogg_usersym, 10, /fill
for i=0L, 5L do begin
    ii=where(abs(colors[i,*]-mcolors[i,*]) lt 1.e+6 )
    typerr=0.05
    if(i ne 4) then begin
        magerr0=sqrt(1./maggies[ioff[i],ii]^2/ivar[ioff[i],ii])
        magerr1=sqrt(1./maggies[2,ii]^2/ivar[2,ii])
        err=sqrt(magerr0^2+magerr1^2)
        typerr=median(err)
    endif 
    yra=[-1., 1.]*5.*typerr
;;    hogg_scatterplot, gz[ii].z, colors[i,ii]-mcolors[i,ii], psym=3, $
;;      xra=[0.02, 2.04], yra=yra, $
;;      xtitle='!8z!6', ytitle=textoidl(ytitle[i]), ycharsize=axis_char_scale, $
;;      /cond, quant=[0.1,0.25, 0.5,0.75, 0.9], exp=0.5, xnpix=8, ynpix=12, $
;;      satfrac=0.001
    djs_plot, gz[ii].z, colors[i,ii]-mcolors[i,ii], psym=8, $
      xra=[0.02, 2.04], yra=yra, symsize=0.3, $
      xtitle='!8z!6', ytitle=textoidl(ytitle[i]), ycharsize=axis_char_scale
    djs_oplot,[-1.,10.], [1.,1.]*typerr, color='red', th=6, linest=2
    djs_oplot,[-1.,10.], -[1.,1.]*typerr, color='red', th=6, linest=2
    djs_oplot,[-1.,10.], [1.,1.]*typerr, color='red', th=6, linest=2
    djs_oplot,[-1.,10.], -[1.,1.]*typerr, color='red', th=6, linest=2
endfor
k_end_print, xold=xold, yold=yold, pold=pold

gphoto=rsex(getenv('KCORRECT_DIR')+ $
            '/data/redshifts/goods/mb_cdfs_isaac_ks_photz_c1.1_d3.0k.cat')
gz=rsex(getenv('KCORRECT_DIR')+ $
        '/data/redshifts/goods/mb_z.cdfs.c1.1z.20050218.cat')
igot=where(gz.z lt 2. and gz.z gt 0.1)
gz=gz[igot]
gphoto=gphoto[igot]
kc=goods_kcorrect(gz.z, goods=gphoto, $
                  omaggies=maggies, rmaggies=rmaggies, band_shift=0., $
                  oivar=ivar)
k_print, filename=getenv('KCORRECT_DIR')+'/docs/paper/goods.ps', $
  xold=xold, yold=yold, pold=pold
!P.MULTI=[6,2,3]
!X.MARGIN=9.
!Y.MARGIN=3.
colors=fltarr(9, n_elements(gphoto))
mcolors=fltarr(9, n_elements(gphoto))
ioff=[0,1,3,4,5,6]
for i=0L, 5L do $
  colors[i,*]=-2.5*alog10(maggies[ioff[i],*]/maggies[2,*])
for i=0L, 5L do $
  mcolors[i,*]=-2.5*alog10(rmaggies[ioff[i],*]/rmaggies[2,*])
ytitle=['\Delta!8(B-i)!6', $
        '\Delta!8(V-i)!6', $
        '\Delta!8(z-i)!6', $
        '\Delta!8(J-i)!6', $
        '\Delta!8(H-i)!6', $
        '\Delta!8(K_s-i)!6']
hogg_usersym, 10, /fill
for i=0L, 5L do begin
    ii=where(abs(colors[i,*]-mcolors[i,*]) lt 1.e+6 )
    typerr=0.05
    if(i ne 4) then begin
        magerr0=sqrt(1./maggies[ioff[i],ii]^2/ivar[ioff[i],ii])
        magerr1=sqrt(1./maggies[2,ii]^2/ivar[2,ii])
        err=sqrt(magerr0^2+magerr1^2)
        typerr=median(err)
    endif 
    yra=[-1., 1.]*5.*typerr
;;    hogg_scatterplot, gz[ii].z, colors[i,ii]-mcolors[i,ii], psym=3, $
;;      xra=[0.02, 2.04], yra=yra, $
;;      xtitle='!8z!6', ytitle=textoidl(ytitle[i]), ycharsize=axis_char_scale, $
;;      /cond, quant=[0.1,0.25, 0.5,0.75, 0.9], exp=0.5, xnpix=8, ynpix=12, $
;;      satfrac=0.001
    djs_plot, gz[ii].z, colors[i,ii]-mcolors[i,ii], psym=8, $
      xra=[0.02, 2.04], yra=yra, symsize=0.3, $
      xtitle='!8z!6', ytitle=textoidl(ytitle[i]), ycharsize=axis_char_scale
    djs_oplot,[-1.,10.], [1.,1.]*typerr, color='red', th=6, linest=2
    djs_oplot,[-1.,10.], -[1.,1.]*typerr, color='red', th=6, linest=2
    djs_oplot,[-1.,10.], [1.,1.]*typerr, color='red', th=6, linest=2
    djs_oplot,[-1.,10.], -[1.,1.]*typerr, color='red', th=6, linest=2
endfor
k_end_print, xold=xold, yold=yold, pold=pold

stop

gstfile=getenv('KCORRECT_DIR')+'/data/test/gst_tests.fits'
cat=mrdfits(gstfile, 1)
kc=gst_kcorrect(cat.z, calibobj=cat, twomass=cat, galex=cat, $
                omaggies=maggies, rmaggies=rmaggies, band_shift=0., $
                oivar=ivar)
colors=fltarr(9, n_elements(cat))
mcolors=fltarr(9, n_elements(cat))
for i=0L, 8L do $
  colors[i,*]=-2.5*alog10(maggies[i,*]/maggies[i+1,*])
for i=0L, 8L do $
  mcolors[i,*]=-2.5*alog10(rmaggies[i,*]/rmaggies[i+1,*])
k_print, filename=getenv('KCORRECT_DIR')+'/docs/paper/fullfits.ps', $
  xold=xold, yold=yold, pold=pold, xsize=9., ysize=6.
!P.MULTI=[0,3,3]
!X.MARGIN=9.
!Y.MARGIN=3.
ytitle=['\Delta!8(F-N)!6', $
        '\Delta!8(N-u)!6', $
        '\Delta!8(u-g)!6', $
        '\Delta!8(g-r)!6', $
        '\Delta!8(r-i)!6', $
        '\Delta!8(i-z)!6', $
        '\Delta!8(z-J)!6', $
        '\Delta!8(J-H)!6', $
        '\Delta!8(H-K_s)!6']
yras=[[-0.99, 0.99], $
      [-0.99, 0.99], $
      [-0.29, 0.29], $
      [-0.19, 0.19], $
      [-0.19, 0.19], $
      [-0.19, 0.19], $
      [-0.29, 0.29], $
      [-0.29, 0.29], $
      [-0.29, 0.29]]
isort=sort(cat.z)
isort=isort[long(findgen(300)*n_elements(isort))/300.]
for i=0L, 8L do begin
    ii=where(abs(colors[i,*]-mcolors[i,*]) lt 1.e+6 and $
             ivar[i,*] gt 0 and ivar[i+1,*] gt 0) 
    magerr0=sqrt(1./maggies[i,ii]^2/ivar[i,ii])
    magerr1=sqrt(1./maggies[i+1,ii]^2/ivar[i+1,ii])
    err=sqrt(magerr0^2+magerr1^2)
    typerr=median(err)
    yra=[-1., 1.]*5.*typerr
    hogg_scatterplot, cat[ii].z, colors[i,ii]-mcolors[i,ii], psym=3, $
      xra=[0.02, 0.24], yra=yra, $
      xtitle='!8z!6', ytitle=textoidl(ytitle[i]), ycharsize=axis_char_scale, $
      /cond, quant=[0.1,0.25, 0.5,0.75, 0.9], exp=0.5, xnpix=20, ynpix=20, $
      satfrac=0.001
    djs_oplot,[-1.,10.], [1.,1.]*typerr, color='red', th=6, linest=2
    djs_oplot,[-1.,10.], -[1.,1.]*typerr, color='red', th=6, linest=2
endfor
k_end_print, xold=xold, yold=yold, pold=pold

cat=mrdfits(getenv('KCORRECT_DIR')+'/data/test/sdss_tests_lrg.fits',1)
kc=sdss_kcorrect(cat.z, calibobj=cat, band_shift=0.3, rmaggies=rmaggies, $
                 omaggies=maggies, oivar=oivar, /lrg)
umg=-2.5*alog10(maggies[0,*]/maggies[1,*])
gmr=-2.5*alog10(maggies[1,*]/maggies[2,*])
rmi=-2.5*alog10(maggies[2,*]/maggies[3,*])
imz=-2.5*alog10(maggies[3,*]/maggies[4,*])
mumg=-2.5*alog10(rmaggies[0,*]/rmaggies[1,*])
mgmr=-2.5*alog10(rmaggies[1,*]/rmaggies[2,*])
mrmi=-2.5*alog10(rmaggies[2,*]/rmaggies[3,*])
mimz=-2.5*alog10(rmaggies[3,*]/rmaggies[4,*])
k_print, filename=getenv('KCORRECT_DIR')+'/docs/paper/lrg_colors.ps', $
  xold=xold, yold=yold, pold=pold
!P.MULTI=[0,2,2]
!X.MARGIN=1.
!Y.MARGIN=0.
ytitle=['!8(u-g)!6', $
        '!8(g-r)!6', $
        '!8(r-i)!6', $
        '!8(i-z)!6']
isort=sort(cat.z)
isort=isort[long(findgen(300)*n_elements(isort))/300.]
hogg_scatterplot, cat.z, umg, psym=3, xra=[0.08, 0.54], yra=[-0.71, 3.2], $
  xcharsize=0.0001, ytitle=ytitle[0], ycharsize=axis_char_scale, $
  /cond, quant=[0.1,0.25, 0.5,0.75, 0.9], exp=0.5, xnpix=30, ynpix=30, $
  satfrac=0.001
djs_oplot, cat[isort].z, mumg[isort], th=6, color='red'
hogg_scatterplot, cat.z, gmr, psym=3, xra=[0.08, 0.54], yra=[0.31, 2.5], $
  xcharsize=0.0001, ycharsize=0.0001, $
  /cond, quant=[0.1,0.25, 0.5,0.75, 0.9], exp=0.5, xnpix=30, ynpix=30, $
  satfrac=0.001
djs_oplot, cat[isort].z, mgmr[isort], th=6, color='red'
axis,!X.CRANGE[1],!Y.CRANGE[0],yaxis=1, $
  ytitle=textoidl(ytitle[1]),ycharsize=axis_char_scale 
hogg_scatterplot, cat.z, rmi, psym=3, xra=[0.08, 0.54], yra=[0.11, 1.19], $
  ytitle=ytitle[2], ycharsize=axis_char_scale, xtitle='!8z!6' , $
  /cond, quant=[0.1,0.25, 0.5,0.75, 0.9], exp=0.5, xnpix=30, ynpix=30, $
  satfrac=0.001
djs_oplot, cat[isort].z, mrmi[isort], th=6, color='red'
hogg_scatterplot, cat.z, imz, psym=3, xra=[0.08, 0.54], yra=[-.11, 1.1], $
  ycharsize=0.0001, xtitle='!8z!6', $
  /cond, quant=[0.1,0.25, 0.5,0.75, 0.9], exp=0.5, xnpix=30, ynpix=30, $
  satfrac=0.001
djs_oplot, cat[isort].z, mimz[isort], th=6, color='red'
axis,!X.CRANGE[1],!Y.CRANGE[0],yaxis=1, $
  ytitle=textoidl(ytitle[3]),ycharsize=axis_char_scale
k_end_print, xold=xold, yold=yold, pold=pold

k_print, filename=getenv('KCORRECT_DIR')+'/docs/paper/sfh_templates.ps', $
  axis_char_scale=axis_char_scale, xsize=9., ysize=7.
!P.MULTI=[15,3,5]
!Y.MARGIN=0
!X.MARGIN=7
sfh_reconstruct, identity(5)*1.e+12, dage=dage, age=age, $
  sfr=sfr, met=met
for i=0L, 4L do begin
    xcharsize=axis_char_scale
    if(i lt 4) then xcharsize=0.0001
    sfr[*,i]=sfr[*,i]*age*alog(10.)
    djs_plot, age, sfr[*,i], th=4, xra=[5.e+5, 18.e+9], /xlog, /ylog, $
      xtitle='lookback time (yrs)', xcharsize=xcharsize, $
      yra=minmax(sfr[*,i])*[0.1,10.0], $
      ytitle='!8\Delta!8S/\Delta!8 !6[log_{10} !8t!6]'
    djs_plot, age, met[*,i], th=4, xra=[5.e+5, 18.e+9], /xlog, $
      yra=[-0.01, 0.068], $
      xtitle='lookback time (yrs)', xcharsize=xcharsize, $
      ytitle='metallicity'
    coeffs=fltarr(5)
    coeffs[i]=1.e+12
    k_reconstruct_spec, coeffs, loglam, flux
    ilam=where(10.^loglam gt 1500. and 10.^loglam lt 11999)
    flux=flux/max(flux[ilam])
    djs_plot, 10.^loglam, flux, xra=[1500., 11999.], $
      yra=[-0.05, 1.15]*max(flux[ilam]), xcharsize=xcharsize, $
      xtitle='!8\lambda!6 (Angstroms)', $
      ytitle='!8f(\lambda)!6 (arbitrary)'
endfor
k_end_print

k_print, filename=getenv('KCORRECT_DIR')+'/docs/paper/spec_lrg.ps'
k_reconstruct_spec, [1.], loglam, flux, vname='lrg1'
djs_plot, 10.^loglam, flux, xra=[1500., 11999.], yra=[-0.05, 1.25]*max(flux), $
  xtitle='!8\lambda!6 (Angstroms)', $
  ytitle='!8f(\lambda)!6 (ergs cm^{-2} s^{-1} \AA^{-1})'
k_end_print

k_print, filename=getenv('KCORRECT_DIR')+'/docs/paper/sfh_lrg.ps'
!P.MULTI=[2,1,2]
!Y.MARGIN=0
sfh_reconstruct, [1.e+12], vname='lrg1', dage=dage, age=age, sfr=sfr, met=met
sfr=sfr*age*alog(10.)
djs_plot, age, sfr, th=4, /xlog, /ylog, $
  xtitle='lookback time (yrs)', xcharsize=0.0001, yra=[1.e-11, 1.e+14], $
  ytitle='!8\Delta!8S/\Delta!8 !6[log_{10} !8t!6]'
djs_plot, age, met, th=4, /xlog, yra=[-0.01, 0.068], $
  xtitle='lookback time (yrs)', $
  ytitle='metallicity'
k_end_print

k_print, filename=getenv('KCORRECT_DIR')+'/docs/paper/spec_lrg.ps'
k_reconstruct_spec, [1.], loglam, flux, vname='lrg1'
djs_plot, 10.^loglam, flux, xra=[1500., 11999.], yra=[-0.05, 1.25]*max(flux), $
  xtitle='!8\lambda!6 (Angstroms)', $
  ytitle='!8f(\lambda)!6 (ergs cm^{-2} s^{-1} \AA^{-1})'
k_end_print

end
;------------------------------------------------------------------------------
