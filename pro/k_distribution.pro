pro k_distribution

va=hogg_mrdfits('/global/data/sdss/sdssgalaxies/devel/object_catalog.fits',1, $
                nrowchunk=10000,columns=['petrocounts','reddening', $
                                         'z_princeton'])

indx=where(va.z_princeton gt 0.09 and va.z_princeton lt 0.11)
va=va[indx]

mags=va.petrocounts-va.reddening
mags[0,*]=mags[0,*]-0.042
mags[1,*]=mags[1,*]+0.036
mags[2,*]=mags[2,*]+0.015
mags[3,*]=mags[3,*]+0.013
mags[4,*]=mags[4,*]-0.002
magserr=dblarr(5,n_elements(va))+0.02
zz=va.z_princeton

kcorrect,mags,magserr,zz,kcorrect,kcorrectz=0.1,version='default', $
  coeff=coeff, ematrix=ematrix, bmatrix=bmatrix, lambda=lambda, $
  /vconstraint,/sdssfix

; then set up to project on filters
spec=(coeff##ematrix)##bmatrix
lambda=0.5*(lambda[0:n_elements(lambda)-2]+lambda[1:n_elements(lambda)-1])


bands=-2.5*alog10(k_project_filters(lambda,spec,band_shift=0.0, $
                                    filterlist=['Bfilter','Vfilter','bj']))

bands[0,*]=bands[0,*]-(k_vega2ab(filterlist=['Bfilter']))[0]
bands[1,*]=bands[1,*]-(k_vega2ab(filterlist=['Vfilter']))[0]
bands[2,*]=bands[2,*]-(k_vega2ab(filterlist=['bj']))[0]

set_print,filename='bjBV.ps'
plot,bands[0,*]-bands[1,*],bands[2,*]-bands[0,*],psym=3,xra=[0.,1.5], $
  xtitle='B-V',ytitle='bj-B'
xx=[-10.,10.]
yy=-0.28*xx
oplot,xx,yy
xx=bands[0,*]-bands[1,*]
yy=bands[2,*]-bands[0,*]-(-0.28*(bands[0,*]-bands[1,*]))
plot,xx,yy,psym=3,xra=[0.,1.5], $
  xtitle='B-V',ytitle='(bj-B)+0.28(B-V)'
end_print
stop

end
