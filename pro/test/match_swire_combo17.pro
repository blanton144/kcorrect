cat=k_read_tbl('SWIRE3_CDFS_cat_IRAC24_21Dec05.tbl')

c17=mrdfits(getenv('PRIMUS_DIR')+ $
            '/data/COMBO17/COMBO17_table3.fits',1) 
ii=where(c17.rmag lt 24. and c17.mc_z gt 0.1)
c17=c17[ii]

spherematch, c17.ra, c17.dec, cat.ra, cat.dec, 3./3600., m1, m2, d12

c17=c17[m1]
cat=cat[m2]

save, filename='swc17.sav'

restore, 'swc17.sav'

swire_to_maggies,cat, mm, iv

umg=((-2.5*alog10(mm[0,*]/mm[1,*])) > (-4.)) < 4.
umr=((-2.5*alog10(mm[0,*]/mm[2,*])) > (-4.)) < 4.
gmr=((-2.5*alog10(mm[1,*]/mm[2,*])) > (-1.)) < 4.
tmf=-2.5*alog10(mm[5,*]/mm[6,*])
fms=((-2.5*alog10(mm[6,*]/mm[7,*])) > (-4.)) < 4.
sme=((-2.5*alog10(mm[7,*]/mm[8,*])) > (-4.)) < 40.
emtw=((-2.5*alog10(mm[8,*]/mm[9,*])) > (-4.)) < 40.

zz=(findgen(100)+0.5)/100.
kc=sdss_kcorrect(zz, coeffs=fltarr(100)+1., /lrg, $
                 rmaggies=rmaggies)
lrggmr=-2.5*alog10(rmaggies[1,*]/rmaggies[2,*])

lowz=lowz_read(sample='dr4')
ii=where(lowz.absmag[1]-lowz.absmag[2] lt 0.55 AND $
         lowz.absmag[1]-lowz.absmag[2] gt 0.45)
im=mrdfits(vagc_name('object_sdss_imaging'),1,row=lowz[ii].object_position)
kc=sdss_kcorrect(lowz[ii].zdist, cal=im, coeffs=coeffs)

sc=total(coeffs,1)
for i=0,4 do $
  coeffs[i,*]=coeffs[i,*]/sc

mcoeffs=fltarr(5,100)
for i=0,4 do $
  mcoeffs[i,*]=median(coeffs[i,*])

kc=sdss_kcorrect(zz, coeffs=mcoeffs, rmaggies=bmaggies)
bluegmr=-2.5*alog10(bmaggies[1,*]/bmaggies[2,*])

dst=read_draine(getenv('KCORRECT_DIR')+'/data/seds/draine/'+ $
                'spec_2.2.dat')

spitzerf='spitzer_irac_ch'+['1', '2', '3', '4']+'.par'
sdssf='sdss_'+['u', 'g', 'r', 'i', 'z']+'0.par'
filters=[sdssf, spitzerf, 'spitzer_mips_24.par']

k_load_vmatrix,vmatrix,vlambda, vname='lrg1'
dstmm=fltarr(10, 100)
lrgmm=fltarr(10, 100)
for i=0L, 99L do begin & $
ll=k_lambda_to_edges(dst.lambda*(1.+zz[i])) & $
ff=dst.flux/(1.+zz[i]) & $
dstmm[*, i]= k_project_filters(ll, ff, filterlist=filters, /sil) & $
lrgmm[*, i]= k_project_filters(vlambda*(1.+zz[i]), vmatrix/(1.+zz[i]), filterlist=filters, /sil) & $
    endfor

mm=lrgmm
mm[5:7,*]=1.*lrgmm[5:7,*]+3.e+15*dstmm[5:7,*]
mm[8:9,*]=0.*lrgmm[8:9,*]+3.e+15*dstmm[8:9,*]

hogg_usersym, 10, /fill
k_print, filename='ssck1.ps', xold=xold, yold=yold, pold=pold
!P.MULTI=[0,1,3]
!Y.MARGIN=0
djs_plot, c17.mc_z, gmr, xra=[0.1, 0.95], psym=8, yra=[0.1, 2.3], $
  xch=0.0001, ytitle='!8g-r!6', charsize=1.7, syms=0.4
djs_oplot, zz, bluegmr, th=7, color='blue'
djs_oplot, zz, lrggmr, th=7, color='red'
djs_plot, c17.mc_z, tmf, xra=[0.1, 0.95], psym=8, yra=[-1.2, 0.7], $
  xch=0.0001, ytitle='!8[3.6]-[4.8]!6', charsize=1.7, syms=0.4
djs_oplot, zz, -2.5*alog10(mm[5,*]/mm[6,*]), th=7, color='red'
djs_plot, c17.mc_z, emtw, xra=[0.1, 0.95], psym=8, yra=[-0.9, 4.3], $
  xch=1.7, ytitle='!8[7.8]-[24]!6', charsize=1.7, syms=0.4, $
  xtitle='!6redshift !8z!6'
djs_oplot, zz, -2.5*alog10(mm[8,*]/mm[9,*]), th=7, color='red'
k_end_print, xold=xold, yold=yold, pold=pold

hogg_usersym, 10, /fill
k_print, filename='ssck2.ps', xold=xold, yold=yold, pold=pold, xsize=7, $
  ysize=9
!P.MULTI=[0,2,2]
!X.MARGIN=1
!Y.MARGIN=4
djs_plot, umg, fms, xra=[-1.1, 1.65], psym=8, yra=[-1.5, 1.4], $
  xtitle='!8u-g!6', ytitle='!8[4.8]-[5.6]!6', syms=0.4
djs_plot, tmf, fms, xra=[-0.8, 0.8], psym=8, yra=[-1.5, 1.4], $
  xtitle='!8[3.6]-[4.8]!6', ytitle='!8[4.8]-[5.6]!6', syms=0.4, $
  ych=0.001
axis,!X.CRANGE[1],!Y.CRANGE[0],yaxis=1, $
  ytitle='!8[4.8]-5.6]!6'
djs_plot, sme, fms, xra=[-2.1, 2.65], psym=8, yra=[-1.5, 1.4], $
  xtitle='!8[5.6]-[7.8]!6', ytitle='!8[4.8]-[5.6]!6', syms=0.4
djs_plot, emtw, fms, xra=[-0.8, 3.2], psym=8, yra=[-1.5, 1.4], $
  xtitle='!8[7.8]-[24]!6', ytitle='!8[4.8]-[5.6]!6', syms=0.4, $
  ych=0.001
axis,!X.CRANGE[1],!Y.CRANGE[0],yaxis=1, $
  ytitle='!8[4.8]-5.6]!6'
k_end_print, xold=xold, yold=yold, pold=pold

dale=k_read_dale()
dale.flux=dale.flux*1.e+30

spitzerf='spitzer_irac_ch'+['1', '2', '3', '4']+'.par'
sdssf='sdss_'+['u', 'g', 'r', 'i', 'z']+'0.par'
filters=[sdssf, spitzerf, 'spitzer_mips_24.par']

dalemm=fltarr(10, 100)
for i=0L, 99L do begin & $
ll=k_lambda_to_edges(dale.lambda*(1.+zz[i])) & $
ff=dale.flux[30]/(1.+zz[i]) & $
dalemm[*, i]= k_project_filters(ll, ff, filterlist=filters, /sil) & $
  endfor

dalegmr=-2.5*alog10(dalemm[1,*]/dalemm[2,*])
daletmf=-2.5*alog10(dalemm[5,*]/dalemm[6,*])
dalefms=-2.5*alog10(dalemm[6,*]/dalemm[7,*])
dalefmtw=-2.5*alog10(dalemm[6,*]/dalemm[9,*])

splot, c17.mc_z, fmtw,psym=4
soplot, zz, dalefmtw, th=4, color='red'

