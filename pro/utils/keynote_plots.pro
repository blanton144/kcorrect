;+
; NAME:
;   keynote_plots
; PURPOSE:
;   make presentation plots
; CALLING SEQUENCE:
;   keynote_plots
; REVISION HISTORY:
;   2005-06-23 MRB, NYU
;-
;------------------------------------------------------------------------------
pro keynote_plots

metallicities=[0.0001, 0.0004, 0.004, 0.008, 0.02, 0.05]


swire=mrdfits(getenv('VAGC_REDUX')+'/spitzer/swire_catalog.fits',1)
objects=mrdfits(getenv('VAGC_REDUX')+'/spitzer/swire_objects.fits',1)
iin=where(objects.object_position ge 0)
sp=mrdfits(vagc_name('object_sdss_spectro'),1,row=objects[iin].object_position)
im=mrdfits(vagc_name('object_sdss_imaging'),1,row=objects[iin].object_position)
igal=where(sp.z gt 0.003 and sp.z lt 0.5)
sp=sp[igal]
im=im[igal]
swire=swire[iin[igal]]
objects=objects[iin[igal]]
swire_to_maggies, swire, swire_maggies, swire_ivar
sdss_to_maggies, sdss_maggies, sdss_ivar, calibobj=im, flux='model'
swire_maggies[0:4,*]=sdss_maggies
swire_ivar[0:4,*]=sdss_ivar

kcorrect, swire_maggies, swire_ivar, sp.z, kc, rmaggies=rmaggies, $
  vname='nmf5swire', filterlist=['sdss_u0.par', $
                                 'sdss_g0.par', $
                                 'sdss_r0.par', $
                                 'sdss_i0.par', $
                                 'sdss_z0.par', $
                                 'spitzer_irac_ch1.par', $
                                 'spitzer_irac_ch2.par', $
                                 'spitzer_irac_ch3.par', $
                                 'spitzer_irac_ch4.par', $
                                 'spitzer_mips_24.par']

gmr= -2.5*alog10(swire_maggies[1,*]/swire_maggies[2,*])
ch1mch2= -2.5*alog10(swire_maggies[5,*]/swire_maggies[6,*])
ch2mch3= -2.5*alog10(swire_maggies[6,*]/swire_maggies[7,*])
ch3mch4= -2.5*alog10(swire_maggies[7,*]/swire_maggies[8,*])
ch4mch5= -2.5*alog10(swire_maggies[8,*]/swire_maggies[9,*])

mgmr= -2.5*alog10(rmaggies[1,*]/rmaggies[2,*])
mch1mch2= -2.5*alog10(rmaggies[5,*]/rmaggies[6,*])
mch2mch3= -2.5*alog10(rmaggies[6,*]/rmaggies[7,*])
mch3mch4= -2.5*alog10(rmaggies[7,*]/rmaggies[8,*])
mch4mch5= -2.5*alog10(rmaggies[8,*]/rmaggies[9,*])

hogg_usersym, 10, /fill
k_print, filename='gmrch1ch2.ps', xold=xold, yold=yold, pold=pold
djs_plot, gmr, ch1mch2, /cond, $
  quantiles=[0.1, 0.25, 0.5, 0.75, 0.9],$
  xtitle='g-r', ytitle='[3.6]-[4.5]', $
  xra=[-0.01, 2.08], yra=[-0.71, 0.08], psym=8, symsize=0.8
djs_oplot, mgmr, mch1mch2, th=5,color='red', psym=8, symsize=0.8
k_end_print, xold=xold, yold=yold, pold=pold

spawn, 'psgif -d 255 gmrch1ch2.ps > gmrch1ch2.gif'

hogg_usersym, 10, /fill
k_print, filename='ch1ch2ch3.ps', xold=xold, yold=yold, pold=pold
djs_plot, ch1mch2, ch2mch3, /cond, $
  quantiles=[0.1, 0.25, 0.5, 0.75, 0.9],$
  ytitle='[4.5]-[5.4]', xtitle='[3.6]-[4.5]', $
  xra=[-0.91, 0.08], yra=[-0.91, 1.08], psym=8, symsize=0.8
djs_oplot, mch1mch2, mch2mch3, th=5,color='red', psym=8, symsize=0.8
k_end_print, xold=xold, yold=yold, pold=pold

spawn, 'psgif -d 255 gmrch1ch2.ps > gmrch1ch2.gif'

sample='drtwo14'
postcat=hogg_mrdfits(vagc_name('post_catalog', sample=sample, $
                               letter='safe', post='1'), 1, nrow=28800)
ii=where((postcat.letter_mask and 32) eq 0)
postcat=postcat[ii]
ii=shuffle_indx(n_elements(postcat), num_sub=10000)
postcat=postcat[ii]
im=mrdfits(vagc_name('object_sdss_imaging'),1,row=postcat.object_position)
spawn, 'date'
photoz=sdss_kphotoz(calibobj=im, vname='nmf3photoz')
spawn, 'date'
k_print, filename='main_photoz.ps', xold=xold, yold=yold, pold=pold
hogg_scatterplot, photoz, postcat.z, /cond, $
  quantiles=[0.1, 0.25, 0.5, 0.75, 0.9],$
  xtitle='photo-z', ytitle='spec-z', $
    exp=0.5, xnpix=35, ynpix=35, $
  satfrac=0.001, xra=[-0.01, 0.28], yra=[-0.01, 0.28]
djs_oplot, [0., 1.], [0.,1.], th=5,color='red'
k_end_print, xold=xold, yold=yold, pold=pold

stop
spawn, 'psgif -d 255 main_photoz.ps > main_photoz.gif'


cat=mrdfits(getenv('KCORRECT_DIR')+'/data/test/sdss_tests_lrg.fits',1)
ii=where(cat.z gt 0.01)
cat=cat[ii]
photoz=sdss_kphotoz(calibobj=cat, $
                    coeffs=coeffs, omaggies=maggies, $
                    rmaggies=rmaggies, /lrg)
k_print, filename='lrg_photoz.ps', xold=xold, yold=yold, pold=pold
hogg_scatterplot, photoz, cat.z, /cond, quantiles=[0.1, 0.25, 0.5, 0.75, 0.9],$
  xtitle='photo-z', ytitle='spec-z', $
    exp=0.2, xnpix=40, ynpix=40, $
  satfrac=0.001, xra=[-0.01, 0.5], yra=[-0.01, 0.5]
djs_oplot, [0., 1.], [0.,1.], th=5,color='red'
k_end_print, xold=xold, yold=yold, pold=pold
spawn, 'psgif -d 255 lrg_photoz.ps > lrg_photoz.gif'

stop

lowz=lowz_read(sample='drtwo14')
im=mrdfits(vagc_name('object_sdss_imaging'),1,row=lowz.object_position)
kc=sdss_kcorrect(lowz.zdist, calibobj=im, absmag=absmag, mets=mets, $
                 b300=b300)
kcb=sdss2bessell(lowz.zdist, calibobj=im, absmag=babsmag, $
                 mtol=mtol, /vega)

k_print, filename='lowz_mets.ps', xold=xold, yold=yold, pold=pold
hogg_scatterplot, absmag[2,*], alog10(mets/0.02), xra=[-12.5, -22.5], $
  yra=[-1.8, 0.7], /outliers, xnpix=50, ynpix=50, exp=0.5, satfrac=0.002, $
  ytitle='!6log!d10!n[metallicity/solar]', xtitle='!8M!dr!n!6'
k_end_print,xold=xold, yold=yold, pold=pold

k_print, filename='lowz_b300.ps', xold=xold, yold=yold, pold=pold
hogg_scatterplot, absmag[2,*], alog10(b300), xra=[-12.5, -22.5], $
  yra=[-3.8, -0.5], /outliers, xnpix=50, ynpix=50, exp=0.7, satfrac=0.004, $
  ytitle='!6log!d10!n[!8b!d300!n!6]', xtitle='!8M!dr!n!6'
k_end_print,xold=xold, yold=yold, pold=pold

k_print, filename='lowz_mtol.ps', xold=xold, yold=yold, pold=pold
hogg_scatterplot, babsmag[1,*]-babsmag[3,*], alog10(mtol[3,*]), $
  xra=[0.3, 1.7], $
  yra=[-0.6, 1.3], /outliers, xnpix=50, ynpix=50, exp=0.7, satfrac=0.004, $
  ytitle='!6log!d10!n[!8M/L!dR!n!6]', xtitle='!8B-R!6'
djs_oplot, [0.4,1.6], -0.82+0.85*[0.4,1.6], th=4, color='red'
k_end_print,xold=xold, yold=yold, pold=pold


k_print, filename='lowz_mtol_s.ps', xold=xold, yold=yold, pold=pold
hogg_scatterplot, babsmag[1,*]-babsmag[3,*], alog10(mtol[3,*]), $
  xra=[0.3, 1.7], $
  yra=[-0.6, 1.3], /outliers, xnpix=50, ynpix=50, exp=0.7, satfrac=0.004, $
  ytitle='!6log!d10!n[!8M/L!dR!n!6]', xtitle='!8B-R!6'
djs_oplot, [0.4,1.6], 0.25-0.82+0.85*[0.4,1.6], th=4, color='red'
k_end_print,xold=xold, yold=yold, pold=pold

spawn, 'psgif -d 255 lowz_mets.ps > lowz_mets.gif'
spawn, 'psgif -d 255 lowz_b300.ps > lowz_b300.gif'
spawn, 'psgif -d 255 lowz_mtol.ps > lowz_mtol.gif'
spawn, 'psgif -d 255 lowz_mtol_s.ps > lowz_mtol_s.gif'

stop


cat=mrdfits(getenv('KCORRECT_DIR')+'/data/test/sdss_tests_lrg.fits',1)
kcorrect=sdss_kcorrect(cat.z, calibobj=cat, $
                       coeffs=coeffs, omaggies=maggies, $
                       rmaggies=rmaggies, /lrg)
sdss_to_maggies, calibobj=cat, maggies
umg=-2.5*alog10(maggies[0,*]/maggies[1,*])
gmr=-2.5*alog10(maggies[1,*]/maggies[2,*])
rmi=-2.5*alog10(maggies[2,*]/maggies[3,*])
imz=-2.5*alog10(maggies[3,*]/maggies[4,*])
mumg=-2.5*alog10(rmaggies[0,*]/rmaggies[1,*])
mgmr=-2.5*alog10(rmaggies[1,*]/rmaggies[2,*])
mrmi=-2.5*alog10(rmaggies[2,*]/rmaggies[3,*])
mimz=-2.5*alog10(rmaggies[3,*]/rmaggies[4,*])
k_print, filename='lrg_colors.ps', xold=xold, yold=yold, pold=pold
!P.MULTI=[0,2,2]
!X.MARGIN=1.
!Y.MARGIN=0.
ytitle=['!8(u-g)!6', $
        '!8(g-r)!6', $
        '!8(r-i)!6', $
        '!8(i-z)!6']
isort=sort(cat.z)
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
spawn, 'psgif -d 255 lrg_colors.ps > lrg_colors.gif'
stop



hdr=headfits(getenv('KCORRECT_DIR')+ $
             '/data/templates/k_nmf_mmatrix.default.fits')
mmatrix=mrdfits(getenv('KCORRECT_DIR')+ $
                '/data/templates/k_nmf_mmatrix.default.fits',0)
lambda=mrdfits(getenv('KCORRECT_DIR')+ $
               '/data/templates/k_nmf_mmatrix.default.fits',1)
dust=mrdfits(getenv('KCORRECT_DIR')+ $
             '/data/templates/k_nmf_mmatrix.default.fits',2)
met=mrdfits(getenv('KCORRECT_DIR')+ $
             '/data/templates/k_nmf_mmatrix.default.fits',3)
age=mrdfits(getenv('KCORRECT_DIR')+ $
             '/data/templates/k_nmf_mmatrix.default.fits',4)
nspec=long(sxpar(hdr, 'NSPEC'))
ndusts=long(sxpar(hdr, 'NDUST'))
nmets=long(sxpar(hdr, 'NMET'))
nages=long(sxpar(hdr, 'NAGE'))

colors=reverse(['red', 'orange', 'magenta', 'dark green', $
                'green', 'cyan', 'blue'])

k_print, filename='mappings_mj.ps', xold=xold, yold=yold, pold=pold
djs_plot, [0], [0], /nodata, xra=[3601.,5199.], yra=[1.e-3, 1.2], $
  xtitle='!8\lambda!6 (Angstroms)', $
  ytitle='!8m_j(\lambda)!6', /ylog
ii=where(lambda[0:nspec-1L] gt 1000. and lambda[0:nspec-1L] lt 11000.)
nb=(size(mmatrix,/dim))[1]
nm=nb-nages*ndusts*nmets
ispec=nages*ndusts*nmets+lindgen(nm)
for k=0L, n_elements(ispec)-1L do begin 
    if((k mod 6) eq 2) then begin
        flux=mmatrix[ii,ispec[k]]
        scale=1./(max(flux)*1.1)
        djs_oplot, lambda[ii], scale*flux, th=3, $
          color=colors[(k/3) mod 7]
    endif
endfor
k_end_print, xold=xold, yold=yold, pold=pold
spawn, 'psgif mappings_mj.ps > mappings_mj.gif'
stop

i=0L
j=2L
k_print, filename='bc_mj_02.ps', xold=xold, yold=yold, pold=pold
djs_plot, [0], [0], /nodata, xra=[999.,11000.], yra=[1.e-3, 1.2], $
  xtitle='!8\lambda!6 (Angstroms)', $
  ytitle='!8m_j(\lambda)!6', /ylog, /xlog
ii=where(lambda[0:nspec-1L] gt 1000. and lambda[0:nspec-1L] lt 11000.)
for k=0L, nages-1L do begin 
    if((k mod 2) eq 0) then begin
        l=i*nmets*nages+j*nages+k
        flux=mmatrix[ii,l]
        scale=1./(max(flux)*1.1)
        djs_oplot, lambda[ii], scale*flux, th=3, $
          color=colors[(k/4) mod 7]
        djs_xyouts, 3300., 0.006, 'no dust', charsize=1.4
        djs_xyouts, 3300., 0.004, 'metallicity = 0.2 solar', charsize=1.4
    endif
endfor
k_end_print, xold=xold, yold=yold, pold=pold
spawn, 'psgif bc_mj_02.ps > bc_mj_02.gif'

i=0L
j=4L
k_print, filename='bc_mj_04.ps', xold=xold, yold=yold, pold=pold
djs_plot, [0], [0], /nodata, xra=[999.,11000.], yra=[1.e-3, 1.2], $
  xtitle='!8\lambda!6 (Angstroms)', $
  ytitle='!8m_j(\lambda)!6', /ylog, /xlog
ii=where(lambda[0:nspec-1L] gt 1000. and lambda[0:nspec-1L] lt 11000.)
for k=0L, nages-1L do begin 
    if((k mod 2) eq 0) then begin
        l=i*nmets*nages+j*nages+k
        flux=mmatrix[ii,l]
        scale=1./(max(flux)*1.1)
        djs_oplot, lambda[ii], scale*flux, th=3, $
          color=colors[(k/4) mod 7]
        djs_xyouts, 3300., 0.006, 'no dust', charsize=1.4
        djs_xyouts, 3300., 0.004, 'metallicity = solar', charsize=1.4
    endif
endfor
k_end_print, xold=xold, yold=yold, pold=pold
spawn, 'psgif bc_mj_04.ps > bc_mj_04.gif'

i=2L
j=4L
k_print, filename='bc_mj_24.ps', xold=xold, yold=yold, pold=pold
djs_plot, [0], [0], /nodata, xra=[999.,11000.], yra=[1.e-3, 1.2], $
  xtitle='!8\lambda!6 (Angstroms)', $
  ytitle='!8m_j(\lambda)!6', /ylog, /xlog
ii=where(lambda[0:nspec-1L] gt 1000. and lambda[0:nspec-1L] lt 11000.)
for k=0L, nages-1L do begin 
    if((k mod 2) eq 0) then begin
        l=i*nmets*nages+j*nages+k
        flux=mmatrix[ii,l]
        scale=1./(max(flux)*1.1)
        djs_oplot, lambda[ii], scale*flux, th=3, $
          color=colors[(k/4) mod 7]
        djs_xyouts, 1150., 0.8, '\tau_V = 3', charsize=1.4
        djs_xyouts, 1150., 0.6, 'metallicity = solar', charsize=1.4
    endif
endfor
k_end_print, xold=xold, yold=yold, pold=pold
spawn, 'psgif bc_mj_24.ps > bc_mj_24.gif'

i=0L
j=5L
k_print, filename='bc_mj_05.ps', xold=xold, yold=yold, pold=pold
djs_plot, [0], [0], /nodata, xra=[999.,11000.], yra=[1.e-3, 1.2], $
  xtitle='!8\lambda!6 (Angstroms)', $
  ytitle='!8m_j(\lambda!6)', /ylog, /xlog
ii=where(lambda[0:nspec-1L] gt 1000. and lambda[0:nspec-1L] lt 11000.)
for k=0L, nages-1L do begin 
    if((k mod 2) eq 0) then begin
        l=i*nmets*nages+j*nages+k
        flux=mmatrix[ii,l]
        scale=1./(max(flux)*1.1)
        djs_oplot, lambda[ii], scale*flux, th=3, $
          color=colors[(k/4) mod 7]
        djs_xyouts, 3300., 0.006, 'no dust', charsize=1.4
        djs_xyouts, 3300., 0.004, 'metallicity = 2.5 solar', charsize=1.4
    endif
endfor
k_end_print, xold=xold, yold=yold, pold=pold
spawn, 'psgif bc_mj_05.ps > bc_mj_05.gif'


stop

hdr=headfits(getenv('KCORRECT_DIR')+ $
             '/data/templates/k_nmf_mmatrix.lrg1.fits')
dust=mrdfits(getenv('KCORRECT_DIR')+ $
             '/data/templates/k_nmf_mmatrix.lrg1.fits',2)
met=mrdfits(getenv('KCORRECT_DIR')+ $
             '/data/templates/k_nmf_mmatrix.lrg1.fits',3)
age=mrdfits(getenv('KCORRECT_DIR')+ $
             '/data/templates/k_nmf_mmatrix.lrg1.fits',4)
ndusts=long(sxpar(hdr, 'NDUST'))
nmets=long(sxpar(hdr, 'NMET'))
nages=long(sxpar(hdr, 'NAGE'))
templates=mrdfits(getenv('KCORRECT_DIR')+ $
                  '/data/templates/k_nmf_soln.lrg1.fits',0)
nt=1L

dage=fltarr(nages)
dage[1:nages-2L]=0.5*(age[2:nages-1L]-age[0:nages-3L])
dage[0]=age[1]-age[0]
dage[nages-1L]=age[nages-1L]-age[nages-2L]
!P.MULTI=[0, 1, 3]
sfh_tot=fltarr(nages, nt)
sfh_met=fltarr(nages, nt)
for i=0L, nt-1L do begin
    sfh=reform(templates[0:nages*nmets*ndusts-1L,i], nages, nmets, ndusts)
    dust=reform(dust, nages, nmets, ndusts)
    for j=0L, nages-1L do $
      sfh_tot[j,i]=total(sfh[j,*,*])
    for j=0L, nages-1L do $
      sfh_met[j,i]=total(sfh[j,*,*]*metallicities[met[j,*,*]])
endfor
sfh_met=sfh_met/sfh_tot

cat=mrdfits(getenv('KCORRECT_DIR')+'/data/test/sdss_tests_lrg.fits',1)
i0=(where(cat.z gt 0.21 and cat.z lt 0.22))[0]
i1=(where(cat.z gt 0.31 and cat.z lt 0.32))[0]
i2=(where(cat.z gt 0.41 and cat.z lt 0.42))[0]
indx=[i0,i1,i2]

for k=0L, n_elements(indx)-1L do begin
    kcorrect=sdss_kcorrect(cat[indx[k]].z, calibobj=cat[indx[k]], $
                           coeffs=coeffs, omaggies=omaggies, /lrg)

    k_print, filename='pres_sfh_lrg'+strtrim(string(k),2)+'.ps', $
      xold=xold, yold=yold, pold=pold
    !P.MULTI=[0,1,2]
    !Y.MARGIN=0
    djs_plot, [0], [0], /nodata, xra=[0.9e+6, 15.e+9], yra=[1.e-9, 1.2], $
      xtitle='!6age (years)', $
      ytitle='!6SFR', /ylog, /xlog, xcharsize=0.0001
    sfh=sfh_tot#coeffs
    scale=1./(max(sfh)*1.1)
    djs_oplot, age[*,0,0], sfh*scale, th=8
    djs_oplot, age[*,0,0], (age[*,0,0]/13.e+9)^2.0, color='red', th=4
    djs_plot, [0], [0], /nodata, xra=[0.9e+6, 15.e+9], yra=[0.0001, 0.06], $
      xtitle='!6age (years)', $
      ytitle='!6metallicity', /xlog
    sfh=sfh_tot#coeffs
    met=((sfh_met*sfh_tot)#coeffs)/sfh
    djs_oplot, age[*,0,0], met, th=8
    k_end_print, xold=xold, yold=yold, pold=pold

    leff=k_lambda_eff()
    absrc=3.631*2.99792*1.e-2/leff^2
    oflux=omaggies*absrc
    leff=leff/(1.+cat[indx[k]].z)
    oflux=oflux*(1.+cat[indx[k]].z)
    k_print, filename='pres_model_lrg'+strtrim(string(k),2)+'.ps', $
      xold=xold, yold=yold, pold=pold
    djs_plot, [0], [0], /nodata, xra=[2000., 11000.], yra=[1.e-2, 1.2], $
      xtitle='!8\lambda!6 (Angstroms)', $
      ytitle='!8f_\lambda(\lambda)!6', /ylog, /xlog
    k_reconstruct_spec, coeffs, loglam, model, vname='lrg1'
    scale=1./(max(oflux)*1.6)
    djs_oplot, leff, oflux*scale, psym=6, th=6
    djs_oplot, 10.^loglam, scale*model, th=4, color='red'
    djs_oplot, leff, oflux*scale, psym=6, th=6
    k_end_print, xold=xold, yold=yold, pold=pold

    if(k eq 0) then $
      spawn, 'psgif pres_sfh_lrg'+strtrim(string(k),2)+'.ps > '+ $
      'pres_sfh_lrg'+strtrim(string(k),2)+'.gif'
    spawn, 'psgif pres_model_lrg'+strtrim(string(k),2)+'.ps > '+ $
      'pres_model_lrg'+strtrim(string(k),2)+'.gif'
endfor
stop

kc=hogg_mrdfits(vagc_name('kcorrect', flux='petro', coll='none', $
                          band_shift='0.10'),1,nrow=28800)
gmr=kc.absmag[1]-kc.absmag[2]
i0=(where(kc.z gt 0.01 and kc.z lt 0.02 and gmr gt 0.7))[0]
i1=(where(kc.z gt 0.11 and kc.z lt 0.12 and gmr lt 0.7))[0]
i2=(where(kc.z gt 0.15 and kc.z lt 0.16 and gmr lt 0.6))[0]
i3=(where(kc.z gt 0.15 and kc.z lt 0.16 and gmr lt 0.8))[0]
i4=(where(kc.z gt 0.21 and kc.z lt 0.22 and gmr lt 0.7))[0]
i5=(where(kc.z gt 0.31 and kc.z lt 0.32 and gmr lt 0.9))[0]
indx=[i0,i1,i2,i3,i4,i5]

hdr=headfits(getenv('KCORRECT_DIR')+ $
             '/data/templates/k_nmf_mmatrix.default.fits')
dust=mrdfits(getenv('KCORRECT_DIR')+ $
             '/data/templates/k_nmf_mmatrix.default.fits',2)
met=mrdfits(getenv('KCORRECT_DIR')+ $
             '/data/templates/k_nmf_mmatrix.default.fits',3)
age=mrdfits(getenv('KCORRECT_DIR')+ $
             '/data/templates/k_nmf_mmatrix.default.fits',4)
ndusts=long(sxpar(hdr, 'NDUST'))
nmets=long(sxpar(hdr, 'NMET'))
nages=long(sxpar(hdr, 'NAGE'))
templates=mrdfits(getenv('KCORRECT_DIR')+ $
                  '/data/templates/k_nmf_soln.default.fits',0)
nt=(size(templates,/dim))[1]
colors=['blue', 'magenta', 'green', 'dark green', 'dark grey']

dage=fltarr(nages)
dage[1:nages-2L]=0.5*(age[2:nages-1L]-age[0:nages-3L])
dage[0]=age[1]-age[0]
dage[nages-1L]=age[nages-1L]-age[nages-2L]
!P.MULTI=[0, 1, 3]
sfh_tot=fltarr(nages, nt)
sfh_met=fltarr(nages, nt)
for i=0L, nt-1L do begin
    sfh=reform(templates[0:nages*nmets*ndusts-1L,i], nages, nmets, ndusts)
    dust=reform(dust, nages, nmets, ndusts)
    for j=0L, nages-1L do $
      sfh_tot[j,i]=total(sfh[j,*,*])
    for j=0L, nages-1L do $
      sfh_met[j,i]=total(sfh[j,*,*]*metallicities[met[j,*,*]])
endfor
sfh_met=sfh_met/sfh_tot

for k=0L, n_elements(indx)-1L do begin
    im=mrdfits(vagc_name('object_sdss_imaging'), 1, row=indx[k])
    kcorrect=sdss_kcorrect(kc[indx[k]].z, calibobj=im, coeffs=coeffs, $
                           omaggies=omaggies)

    k_print, filename='pres_sfh_photo'+strtrim(string(k),2)+'.ps', $
      xold=xold, yold=yold, pold=pold
    !P.MULTI=[0,1,2]
    !Y.MARGIN=0
    djs_plot, [0], [0], /nodata, xra=[0.9e+6, 15.e+9], yra=[1.e-9, 1.2], $
      xtitle='!6age (years)', $
      ytitle='!6SFR', /ylog, /xlog, xcharsize=0.0001
    sfh=sfh_tot#coeffs
    scale=1./(max(sfh)*1.1)
    djs_oplot, age[*,0,0], sfh*scale, th=8
    for i=0L, n_elements(coeffs)-1L do begin 
        scoeffs=0.*coeffs 
        scoeffs[i]=coeffs[i] 
        sfh=sfh_tot#scoeffs
        djs_oplot, age[*,0,0], scale*sfh, th=4, color=colors[i] 
    endfor 
    djs_plot, [0], [0], /nodata, xra=[0.9e+6, 15.e+9], yra=[0.0001, 0.06], $
      xtitle='!6age (years)', $
      ytitle='!6metallicity', /xlog
    sfh=sfh_tot#coeffs
    met=((sfh_met*sfh_tot)#coeffs)/sfh
    djs_oplot, age[*,0,0], met, th=8
    for i=0L, n_elements(coeffs)-1L do begin 
        scoeffs=0.*coeffs 
        scoeffs[i]=coeffs[i]
        sfh=sfh_tot#scoeffs
        met=((sfh_met*sfh_tot)#scoeffs)/sfh
        djs_oplot, age[*,0,0], met, th=4, color=colors[i] 
    endfor 
    k_end_print, xold=xold, yold=yold, pold=pold

    leff=k_lambda_eff()
    absrc=3.631*2.99792*1.e-2/leff^2
    oflux=omaggies*absrc
    leff=leff/(1.+kc[indx[k]].z)
    oflux=oflux*(1.+kc[indx[k]].z)
    k_print, filename='pres_model_photo'+strtrim(string(k),2)+'.ps', $
      xold=xold, yold=yold, pold=pold
    djs_plot, [0], [0], /nodata, xra=[2000., 11000.], yra=[1.e-2, 1.2], $
      xtitle='!8\lambda!6 (Angstroms)', $
      ytitle='!8f_\lambda(\lambda)!6', /ylog, /xlog
    k_reconstruct_spec, coeffs, loglam, model 
    scale=1./(max(oflux)*1.6)
    djs_oplot, leff, oflux*scale, psym=6, th=6
    djs_oplot, 10.^loglam, scale*model, th=4, color='red'
    for i=0L, n_elements(coeffs)-1L do begin 
        scoeffs=0.*coeffs 
        scoeffs[i]=coeffs[i] 
        k_reconstruct_spec, scoeffs, loglam, smodel 
        djs_oplot, 10.^loglam, scale*smodel, th=3, color=colors[i] 
    endfor 
    djs_oplot, leff, oflux*scale, psym=6, th=6
    k_end_print, xold=xold, yold=yold, pold=pold

    spawn, 'psgif pres_sfh_photo'+strtrim(string(k),2)+'.ps > '+ $
      'pres_sfh_photo'+strtrim(string(k),2)+'.gif'
;    spawn, 'psgif pres_model_photo'+strtrim(string(k),2)+'.ps > '+ $
;      'pres_model_photo'+strtrim(string(k),2)+'.gif'
endfor
stop

fiberids=[35, 3, 4, 8, 22, 29, 56]
for k=0L, 6L do begin

    fit_sdss_fiber, 401, fiberids[k], flux=flux, loglam=loglam, ivar=ivar, $
      coeffs=coeffs, vdisp=500., model=model

    k_print, filename='pres_model_spec'+strtrim(string(k),2)+'.ps', $
      xold=xold, yold=yold, pold=pold
    djs_plot, [0], [0], /nodata, xra=[2000., 11000.], yra=[1.e-2, 1.2], $
      xtitle='!8\lambda!6 (Angstroms)', $
      ytitle='!8f_\lambda(\lambda)!6', /ylog, /xlog
    ii=where(ivar gt 0)
    scale=1./(max(flux)*1.1)
    djs_oplot, 10.^loglam[ii], scale*flux[ii], th=6
    k_end_print, xold=xold, yold=yold, pold=pold

    k_print, filename='pres_model_spec'+strtrim(string(k),2)+'_m0.ps', $
      xold=xold, yold=yold, pold=pold
    djs_plot, [0], [0], /nodata, xra=[2000., 11000.], yra=[1.e-2, 1.2], $
      xtitle='!8\lambda!6 (Angstroms)', $
      ytitle='!8f_\lambda(\lambda)!6', /ylog, /xlog
    ii=where(ivar gt 0)
    scale=1./(max(flux)*1.1)
    djs_oplot, 10.^loglam[ii], scale*flux[ii], th=6
    djs_oplot, 10.^loglam, scale*model, th=4, color='red'
    k_end_print, xold=xold, yold=yold, pold=pold

    colors=['blue', 'magenta', 'green', 'dark green', 'dark grey']
    isort=lindgen(5)
    for j=1, n_elements(coeffs) do begin 
        k_print, filename='pres_model_spec'+strtrim(string(k),2)+'_m'+$
          strtrim(string(j),2)+'.ps', $
          xold=xold, yold=yold, pold=pold 
        djs_plot, [0], [0], /nodata, xra=[2000., 11000.], yra=[1.e-2, 1.2], $
          xtitle='!8\lambda!6 (Angstroms)', $
          ytitle='!8f_\lambda(\lambda)!6', /ylog, /xlog 
        ii=where(ivar gt 0) 
        scale=1./(max(flux)*1.1) 
        djs_oplot, 10.^loglam[ii], scale*flux[ii], th=6 
        djs_oplot, 10.^loglam, scale*model, th=4, color='red' 
        for i=0L, j-1L do begin 
            scoeffs=0.*coeffs 
            scoeffs[isort[i]]=coeffs[isort[i]] 
            k_reconstruct_spec, scoeffs, loglam, smodel 
            djs_oplot, 10.^loglam, scale*smodel, th=3, color=colors[i] 
        endfor 
        k_end_print, xold=xold, yold=yold, pold=pold 
    endfor

    if(k eq 0) then begin
        spawn, 'psgif pres_model_spec'+strtrim(string(k),2)+ $
          '.ps > pres_model_spec'+strtrim(string(k),2)+ $
          '.gif'
        spawn, 'psgif pres_model_spec'+strtrim(string(k),2)+'_m'+ $
          strtrim(string(0),2)+'.ps > pres_model_spec'+ $
          strtrim(string(k),2)+'_m'+ $
          strtrim(string(0),2)+'.gif'
        for j=1, n_elements(coeffs) do $
          spawn, 'psgif pres_model_spec'+strtrim(string(k),2)+'_m'+ $
          strtrim(string(j),2)+'.ps > pres_model_spec'+ $
          strtrim(string(k),2)+'_m'+ $
          strtrim(string(j),2)+'.gif'
    endif else begin
        spawn, 'psgif pres_model_spec'+strtrim(string(k),2)+'_m'+ $
          strtrim(string(5),2)+'.ps > pres_model_spec'+ $
          strtrim(string(k),2)+'_m'+ $
          strtrim(string(5),2)+'.gif'
    endelse

endfor 
stop

sp=sdss_spectro_matched(columns=['plate', 'fiberid', 'mjd', 'z', 'class'])

i0=(where(sp.plate eq 401 and strtrim(sp.class, 2) eq 'GALAXY' and $
          sp.z gt 0.01 and sp.z lt 0.02))[0]
i1=(where(sp.plate eq 401 and strtrim(sp.class, 2) eq 'GALAXY' and $
          sp.z gt 0.11 and sp.z lt 0.12))[0]
i2=(where(sp.plate eq 401 and strtrim(sp.class, 2) eq 'GALAXY' and $
          sp.z gt 0.21 and sp.z lt 0.22))[0]
i3=(where(sp.plate eq 401 and strtrim(sp.class, 2) eq 'GALAXY' and $
          sp.z gt 0.31 and sp.z lt 0.33))[0]
i4=(where(sp.plate eq 401 and strtrim(sp.class, 2) eq 'GALAXY' and $
          sp.z gt 0.41 and sp.z lt 0.42))[0]

indx=[i0, i2, i4]

avloglam=double(alog10(2500.)+(alog10(9999.)-alog10(2500.))* $
                (dindgen(4000)+0.5)/4000.)
sdss_spec_block, sp[indx].plate, sp[indx].fiberid, sp[indx].mjd, $
  block_flux=flux, block_lambda=lambda, block_ivar=ivar, $
  avloglam=avloglam

k_print, filename='pres_rf.ps', pold=pold, xold=xold, yold=yold
djs_plot, [0], [0], /nodata, xra=[2500., 9999], yra=[-0.05, 2.65], $
  xtitle='!8\lambda!6 (Angstroms)', $
  ytitle='!8f_\lambda(\lambda)!6'
colors=['blue', 'green', 'red']
scales=[3., 2., 0.8]
for i=0L, n_elements(indx)-1L do begin 
    outflux=flux[*,i] 
    outflux=k_smooth(alog10(lambda), outflux, 600.) 
    maxflux=max(outflux) 
    outflux=outflux/maxflux*scales[i] 
    ii=where(ivar[*,i] gt 0.) 
    djs_oplot, lambda[ii], outflux[ii], th=4, color=colors[i]  
endfor
djs_xyouts, 6000., 0.5, '!8z=0.4!6', color='red', charsize=1.4
djs_xyouts, 7000., 0.8, '!8z=0.2!6', color='green', charsize=1.4
djs_xyouts, 8500., 1.3, '!8z=0.0!6', color='blue', charsize=1.4
k_end_print, pold=pold, xold=xold, yold=yold

k_print, filename='pres_rf_b0.ps', pold=pold, xold=xold, yold=yold
djs_plot, [0], [0], /nodata, xra=[2500., 9999], yra=[-0.05, 2.65], $
  xtitle='!8\lambda!6 (Angstroms)', $
  ytitle='!8f_\lambda(\lambda)!6'
colors=['blue', 'green', 'red']
scales=[3., 2., 0.8]
for i=0L, n_elements(indx)-1L do begin 
    outflux=flux[*,i] 
    outflux=k_smooth(alog10(lambda), outflux, 600.) 
    maxflux=max(outflux) 
    outflux=outflux/maxflux*scales[i] 
    ii=where(ivar[*,i] gt 0.) 
    djs_oplot, lambda[ii], outflux[ii], th=4, color=colors[i]  
endfor
djs_xyouts, 6000., 0.5, '!8z=0.4!6', color='red', charsize=1.4
djs_xyouts, 7000., 0.8, '!8z=0.2!6', color='green', charsize=1.4
djs_xyouts, 8500., 1.3, '!8z=0.0!6', color='blue', charsize=1.4
k_load_filters, ['sdss_u0.par', 'sdss_g0.par', 'sdss_r0.par', $
                 'sdss_i0.par', 'sdss_z0.par'], fn, fl, fp
for i=0L, n_elements(fn)-1L do $
  djs_oplot, fl[0:fn[i]-1,i], fp[0:fn[i]-1,i]/max(fp[*,i])*0.4, th=4, $
  color='blue'
k_end_print, pold=pold, xold=xold, yold=yold

k_print, filename='pres_rf_b1.ps', pold=pold, xold=xold, yold=yold
djs_plot, [0], [0], /nodata, xra=[2500., 9999], yra=[-0.05, 2.65], $
  xtitle='!8\lambda!6 (Angstroms)', $
  ytitle='!8f_\lambda(\lambda)!6'
colors=['blue', 'green', 'red']
scales=[3., 2., 0.8]
for i=0L, n_elements(indx)-1L do begin 
    outflux=flux[*,i] 
    outflux=k_smooth(alog10(lambda), outflux, 600.) 
    maxflux=max(outflux) 
    outflux=outflux/maxflux*scales[i] 	
    ii=where(ivar[*,i] gt 0.) 
    djs_oplot, lambda[ii], outflux[ii], th=4, color=colors[i]  
endfor
djs_xyouts, 6000., 0.5, '!8z=0.4!6', color='red', charsize=1.4
djs_xyouts, 7000., 0.8, '!8z=0.2!6', color='green', charsize=1.4
djs_xyouts, 8500., 1.3, '!8z=0.0!6', color='blue', charsize=1.4
k_load_filters, ['sdss_u0.par', 'sdss_g0.par', 'sdss_r0.par', $
                 'sdss_i0.par', 'sdss_z0.par'], fn, fl, fp
for i=0L, n_elements(fn)-1L do $
  djs_oplot, fl[0:fn[i]-1,i]/1.2, fp[0:fn[i]-1,i]/max(fp[*,i])*0.4, th=4, $
  color='green'
k_end_print, pold=pold, xold=xold, yold=yold

k_print, filename='pres_rf_b2.ps', pold=pold, xold=xold, yold=yold
djs_plot, [0], [0], /nodata, xra=[2500., 9999], yra=[-0.05, 2.65], $
  xtitle='!8\lambda!6 (Angstroms)', $
  ytitle='!8f_\lambda(\lambda)!6'
colors=['blue', 'green', 'red']
scales=[3., 2., 0.8]
for i=0L, n_elements(indx)-1L do begin 
    outflux=flux[*,i] 
    outflux=k_smooth(alog10(lambda), outflux, 600.) 
    maxflux=max(outflux) 
    outflux=outflux/maxflux*scales[i] 
    ii=where(ivar[*,i] gt 0.) 
    djs_oplot, lambda[ii], outflux[ii], th=4, color=colors[i]  
endfor
djs_xyouts, 6000., 0.5, '!8z=0.4!6', color='red', charsize=1.4
djs_xyouts, 7000., 0.8, '!8z=0.2!6', color='green', charsize=1.4
djs_xyouts, 8500., 1.3, '!8z=0.0!6', color='blue', charsize=1.4
k_load_filters, ['sdss_u0.par', 'sdss_g0.par', 'sdss_r0.par', $
                 'sdss_i0.par', 'sdss_z0.par'], fn, fl, fp
for i=0L, n_elements(fn)-1L do $
  djs_oplot, fl[0:fn[i]-1,i]/1.4, fp[0:fn[i]-1,i]/max(fp[*,i])*0.4, th=4, $
  color='red'
k_end_print, pold=pold, xold=xold, yold=yold

spawn, 'psgif pres_rf.ps > pres_rf.gif'
spawn, 'psgif pres_rf_b0.ps > pres_rf_b0.gif'
spawn, 'psgif pres_rf_b1.ps > pres_rf_b1.gif'
spawn, 'psgif pres_rf_b2.ps > pres_rf_b2.gif'

end
