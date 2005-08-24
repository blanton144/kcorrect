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

twomassfile=getenv('KCORRECT_DIR')+'/data/test/twomass_tests.fits'
cat=mrdfits(twomassfile,1)
kc1=twomass_kcorrect(cat.z, calibobj=cat, band_shift=0.1, rmaggies=rmaggies1)

kc0=twomass_kcorrect(cat.z, twomass=cat, calibobj=cat, band_shift=0.1, $
                     rmaggies=rmaggies0, omaggies=omaggies, oivar=oivar)

cresid=fltarr(7, n_elements(cat))
for i=0L, 6L do $
  cresid[i,*]=(-2.5*alog10(rmaggies0[i,*]/rmaggies0[i+1,*]))- $
  (-2.5*alog10(omaggies[i,*]/omaggies[i+1,*]))

k_print, filename=getenv('KCORRECT_DIR')+'/docs/paper/twomass_resid.ps', $
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

k_print, filename=getenv('KCORRECT_DIR')+'/docs/paper/twomass_predicted.ps', $
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

stop

plates=[401, 401, 401, 401, $
        401, 401, 401, 401]
fibers=[101, 121, 501, 533, $
        111, 171, 45, 342]

k_print, filename=getenv('KCORRECT_DIR')+'/docs/paper/specfit.ps', $
  xold=xold, yold=yold, pold=pold, axis_char_scale=axis_char_scale
!P.MULTI=[8,2,4]
!X.MARGIN=0.
!Y.MARGIN=0.
for i=0L, n_elements(plates)-1L do begin
    xcharsize=0.0001
    ycharsize=0.0001
    if(i ge n_elements(plates)-2L) then xcharsize=1.5*axis_char_scale
    if((i mod 2) eq 0) then ycharsize=1.5*axis_char_scale
    readspec, plates[i], fibers[i], flux=flux, wave=wave, zans=zans
    alam=ext_ccm(wave)
    glactc, zans.plug_ra, zans.plug_dec, 2000., gl, gb, 1, /deg
    ebv=dust_getval(gl, gb)
    extvoebv=3.1
    ext=ebv*alam*extvoebv
    flux=flux*10.^(0.4*ext)
    flux=flux*1.e-17
    lambda=k_lambda_to_edges(wave)
    tmp_maggies=k_project_filters(lambda, flux, filterlist=['sdss_g0.par', $
                                                            'sdss_r0.par', $
                                                            'sdss_i0.par'])
    tmp_maggies=reform(tmp_maggies, n_elements(tmp_maggies))
    nmgy=[0., tmp_maggies*1.e+9, 0.]
    ivar=[0., 1./(0.05*tmp_maggies*1.e+9)^2, 0.]
    k_load_vmatrix, vm, la
    kc=sdss_kcorrect(zans.z, nmgy=nmgy, ivar=ivar, coeffs=coeffs)
                     
    mlambda=k_lambda_to_centers(la)
    mflux=vm#coeffs
    mlambda=mlambda*(1.+zans.z)
    mflux=mflux/(1.+zans.z)
    mflux=1.e+17*k_smooth(alog10(mlambda), mflux, 500.)
    flux=1.e+17*k_smooth(alog10(lambda), flux, 500.)
    fscale=max(flux)
    flux=flux/fscale
    mflux=mflux/fscale
    djs_plot, lambda, flux, th=8, xra=[3500., 9600.], $
      xtitle='\lambda (\AA)', $
      ytitle='!8f(\lambda!8)!6', $
      xcharsize=xcharsize, yra=max(flux)*[-0.05,1.15], ycharsize=ycharsize
    djs_oplot, mlambda, mflux, th=4, color='red'
endfor
k_end_print
                 
stop

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
