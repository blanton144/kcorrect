kcorrect=hogg_mrdfits(vagc_name('kcorrect',collision_type='none', $
                                flux_type='petro',band_shift=0.1),1, $
                      nrowchunk=20000)
im=hogg_mrdfits(vagc_name('object_sdss_imaging'),1,nrowchunk=20000 )
sersic=hogg_mrdfits(getenv('VAGC_REDUX')+'/sersic/sersic_catalog.fits',1, $
                    nrowchunk=20000)
calibobj=retrieve_calibobj(im, columns=['ab_exp','ab_dev','fracpsf'])

restore, 'data_test_christi.sav'
spec=mrdfits('spectro_bl*',1)
sersic=hogg_mrdfits(getenv('VAGC_REDUX')+'/sersic/sersic_catalog.fits',1, $
                    nrowchunk=20000)

spherematch,im.ra,im.dec,spec.ra,spec.dec,2./3600., m1,m2,d12

ii=where(calibobj[m1].fracpsf[2] lt 0.3 and kcorrect[m1].abmaggies[7] gt 0.)

outumk= $
  -2.5*alog10(kcorrect[m1[ii]].abmaggies[0]/kcorrect[m1[ii]].abmaggies[7])- $
  kcorrect[m1[ii]].kcorrect[0]+kcorrect[m1[ii]].kcorrect[7]
outgmk= $
  -2.5*alog10(kcorrect[m1[ii]].abmaggies[1]/kcorrect[m1[ii]].abmaggies[7])- $
  kcorrect[m1[ii]].kcorrect[1]+kcorrect[m1[ii]].kcorrect[7]
outrmk= $
  -2.5*alog10(kcorrect[m1[ii]].abmaggies[2]/kcorrect[m1[ii]].abmaggies[7])- $
  kcorrect[m1[ii]].kcorrect[2]+kcorrect[m1[ii]].kcorrect[7]
outimk= $
  -2.5*alog10(kcorrect[m1[ii]].abmaggies[3]/kcorrect[m1[ii]].abmaggies[7])- $
  kcorrect[m1[ii]].kcorrect[3]+kcorrect[m1[ii]].kcorrect[7]
outzmk= $
  -2.5*alog10(kcorrect[m1[ii]].abmaggies[4]/kcorrect[m1[ii]].abmaggies[7])- $
  kcorrect[m1[ii]].kcorrect[4]+kcorrect[m1[ii]].kcorrect[7]
outjmk= $
  -2.5*alog10(kcorrect[m1[ii]].abmaggies[5]/kcorrect[m1[ii]].abmaggies[7])- $
  kcorrect[m1[ii]].kcorrect[5]+kcorrect[m1[ii]].kcorrect[7]
outhmk= $
  -2.5*alog10(kcorrect[m1[ii]].abmaggies[6]/kcorrect[m1[ii]].abmaggies[7])- $
  kcorrect[m1[ii]].kcorrect[6]+kcorrect[m1[ii]].kcorrect[7]

leff=k_lambda_eff(filterlist=['sdss_u0.par', $
                              'sdss_g0.par', $
                              'sdss_r0.par', $
                              'sdss_i0.par', $
                              'sdss_z0.par', $
                              'twomass_J.par', $
                              'twomass_H.par', $
                              'twomass_Ks.par'])

red_fac = (leff/5500.)^(-0.7)
ext_cont=fltarr(8,n_elements(spec))
for i=0, 7 do ext_cont[i,*]=spec.tauv_cont*red_fac[i]
corrmaggies=kcorrect[m1].abmaggies
for i=0, 7 do corrmaggies[i,*]=corrmaggies[i,*]*10.^(0.4*ext_cont[i,m2])

outcorrumk= $
  -2.5*alog10(corrmaggies[0,ii]/corrmaggies[7,ii])- $
  kcorrect[m1[ii]].kcorrect[0]+kcorrect[m1[ii]].kcorrect[7]
outcorrgmk= $
  -2.5*alog10(corrmaggies[1,ii]/corrmaggies[7,ii])- $
  kcorrect[m1[ii]].kcorrect[1]+kcorrect[m1[ii]].kcorrect[7]
outcorrrmk= $
  -2.5*alog10(corrmaggies[2,ii]/corrmaggies[7,ii])- $
  kcorrect[m1[ii]].kcorrect[2]+kcorrect[m1[ii]].kcorrect[7]
outcorrimk= $
  -2.5*alog10(corrmaggies[3,ii]/corrmaggies[7,ii])- $
  kcorrect[m1[ii]].kcorrect[3]+kcorrect[m1[ii]].kcorrect[7]
outcorrzmk= $
  -2.5*alog10(corrmaggies[4,ii]/corrmaggies[7,ii])- $
  kcorrect[m1[ii]].kcorrect[4]+kcorrect[m1[ii]].kcorrect[7]
outcorrjmk= $
  -2.5*alog10(corrmaggies[5,ii]/corrmaggies[7,ii])- $
  kcorrect[m1[ii]].kcorrect[5]+kcorrect[m1[ii]].kcorrect[7]
outcorrhmk= $
  -2.5*alog10(corrmaggies[6,ii]/corrmaggies[7,ii])- $
  kcorrect[m1[ii]].kcorrect[6]+kcorrect[m1[ii]].kcorrect[7]

umkstr='!8!u0.1!n(u-K!ds!n)!6'
gmkstr='!8!u0.1!n(g-K!ds!n)!6'
rmkstr='!8!u0.1!n(r-K!ds!n)!6'
imkstr='!8!u0.1!n(i-K!ds!n)!6'
zmkstr='!8!u0.1!n(z-K!ds!n)!6'
jmkstr='!8!u0.1!n(J-K!ds!n)!6'
hmkstr='!8!u0.1!n(H-K!ds!n)!6'
ba='!8b/a!6'

set_print,filename='test_christi.ps'

hogg_scatterplot,calibobj[m1[ii]].ab_exp[2], outumk, /cond, satfrac=0.01, $
  xra=[0.01,0.99], yra=[0.5,4.9], xtitle=ba, ytitle=umkstr
hogg_scatterplot,calibobj[m1[ii]].ab_exp[2], outcorrumk, /cond, satfrac=0.01, $
  xra=[0.01,0.99], yra=[0.5,4.9], xtitle=ba, ytitle=umkstr
hogg_scatterplot,calibobj[m1[ii]].ab_exp[2], outgmk, /cond, satfrac=0.01, $
  xra=[0.01,0.99], yra=[0.,3.4], xtitle=ba, ytitle=gmkstr
hogg_scatterplot,calibobj[m1[ii]].ab_exp[2], outcorrgmk, /cond, satfrac=0.01, $
  xra=[0.01,0.99], yra=[0.,3.4], xtitle=ba, ytitle=gmkstr
hogg_scatterplot,calibobj[m1[ii]].ab_exp[2], outrmk, /cond, satfrac=0.01, $
  xra=[0.01,0.99], yra=[-0.5,2.3], xtitle=ba, ytitle=rmkstr
hogg_scatterplot,calibobj[m1[ii]].ab_exp[2], outcorrrmk, /cond, satfrac=0.01, $
  xra=[0.01,0.99], yra=[-0.5,2.3], xtitle=ba, ytitle=rmkstr
hogg_scatterplot,calibobj[m1[ii]].ab_exp[2], outimk, /cond, satfrac=0.01, $
  xra=[0.01,0.99], yra=[-0.5,1.7], xtitle=ba, ytitle=imkstr
hogg_scatterplot,calibobj[m1[ii]].ab_exp[2], outcorrimk, /cond, satfrac=0.01, $
  xra=[0.01,0.99], yra=[-0.5,1.7], xtitle=ba, ytitle=imkstr
hogg_scatterplot,calibobj[m1[ii]].ab_exp[2], outzmk, /cond, satfrac=0.01, $
  xra=[0.01,0.99], yra=[-0.8,1.3], xtitle=ba, ytitle=zmkstr
hogg_scatterplot,calibobj[m1[ii]].ab_exp[2], outcorrzmk, /cond, satfrac=0.01, $
  xra=[0.01,0.99], yra=[-0.8,1.3], xtitle=ba, ytitle=zmkstr
hogg_scatterplot,calibobj[m1[ii]].ab_exp[2], outjmk, /cond, satfrac=0.01, $
  xra=[0.01,0.99], yra=[-0.5,0.9], xtitle=ba, ytitle=jmkstr
hogg_scatterplot,calibobj[m1[ii]].ab_exp[2], outcorrjmk, /cond, satfrac=0.01, $
  xra=[0.01,0.99], yra=[-0.5,0.9], xtitle=ba, ytitle=jmkstr
hogg_scatterplot,calibobj[m1[ii]].ab_exp[2], outhmk, /cond, satfrac=0.01, $
  xra=[0.01,0.99], yra=[-0.5,0.6], xtitle=ba, ytitle=hmkstr
hogg_scatterplot,calibobj[m1[ii]].ab_exp[2], outcorrhmk, /cond, satfrac=0.01, $
  xra=[0.01,0.99], yra=[-0.5,0.6], xtitle=ba, ytitle=hmkstr

end_print

jj=where((im[m1].vagc_select and 4) gt 0 AND $
         kcorrect[m1].z gt 0.01 and $
         kcorrect[m1].z lt 0.25)

gmr= $
  -2.5*alog10(kcorrect[m1[jj]].abmaggies[1]/kcorrect[m1[jj]].abmaggies[2])- $
  kcorrect[m1[jj]].kcorrect[1]+kcorrect[m1[jj]].kcorrect[2]
nn=sersic[m1[jj]].sersic_n[3]
sb=22.5-2.5*alog10(kcorrect[m1[jj]].abmaggies[3])- $
  kcorrect[m1[jj]].kcorrect[3]- $
  10.*alog10(1.+spec[m2[jj]].z)+ $
  2.5*alog10(2.*3.14159265*sersic[m1[jj]].sersic_r50[3]^2)
absm=kcorrect[m1[jj]].absmag[3]
  
data=fltarr(4,n_elements(jj))
data[0,*]=gmr
data[1,*]=sb
data[2,*]=nn
data[3,*]=absm


range=[[0.01,1.19], $
       [17.,24.], $
       [0.2,6.1], $
       [-23.5,-16.5]]

hogg_manyd_scatterplot, fltarr(n_elements(jj))+1., data, 'manyd.ps', $
  range=range, xnpix=40, ynpix=40, exp=0.5
hogg_manyd_scatterplot, fltarr(n_elements(jj))+1., data, $
  'manyd-cond.ps', $
  range=range, xnpix=40, ynpix=40, /cond, satfrac=0.005, exp=0.5
  
gmr= $
  -2.5*alog10(corrmaggies[1,jj]/corrmaggies[2,jj])- $
  kcorrect[m1[jj]].kcorrect[1]+kcorrect[m1[jj]].kcorrect[2]
nn=sersic[m1[jj]].sersic_n[3]
sb=22.5-2.5*alog10(corrmaggies[3,jj])- $
  kcorrect[m1[jj]].kcorrect[3]- $
  10.*alog10(1.+spec[m2[jj]].z)+ $
  2.5*alog10(2.*3.14159265*sersic[m1[jj]].sersic_r50[3]^2)
absm=kcorrect[m1[jj]].absmag[3]-ext_cont[3,m2[jj]]
  
data=fltarr(4,n_elements(jj))
data[0,*]=gmr
data[1,*]=sb
data[2,*]=nn
data[3,*]=absm


range=[[0.01,1.19], $
       [17.,24.], $
       [0.2,6.1], $
       [-23.5,-16.5]]

hogg_manyd_scatterplot, fltarr(n_elements(jj))+1., data, 'manyd-dustcorr.ps', $
  range=range, xnpix=40, ynpix=40, exp=0.5
hogg_manyd_scatterplot, fltarr(n_elements(jj))+1., data, $
  'manyd-dustcorr-cond.ps', $
  range=range, xnpix=40, ynpix=40, /cond, satfrac=0.005, exp=0.5
