pro test_kphotoz,failed

if(NOT file_test(getenv('KCORRECT_DIR')+'/test/sample.fits')) then begin
    failed=1
    return
endif

data=mrdfits(getenv('KCORRECT_DIR')+'/test/sample.fits',1)
print,systime()
kphotoz, data.maggies, data.maggies_ivar, redshift, $
  filterlist=['sdss_u0.par','sdss_g0.par','sdss_r0.par','sdss_i0.par', $
              'sdss_z0.par'], /noprior
; , vfile='vmatrix.test2.dat', lfile='lambda.test2.dat'
print,systime()
k_print,filename=getenv('KCORRECT_DIR')+'/test/sample_test_kphotoz.ps', $
  pold=pold,xold=xold,yold=yold
plot, data.redshift, redshift, psym=4
k_end_print,pold=pold,xold=xold,yold=yold

end
a=hogg_mrdfits(getenv('SPECTRO_DATA')+'/spAll.fits',1,nrowchunk=5000,columns=['z','modelflux','modelflux_ivar','class'],ra=[0,30000])
ii=where(a.class eq 'GALAXY' and a.modelflux[0] gt 0.)
kphotoz,a[ii].modelflux,a[ii].modelflux_ivar,photoz,vfile='vmatrix.test2.dat', $
lfile='lambda.test2.dat'
