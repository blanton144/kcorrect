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
print,systime()
k_print,filename=getenv('KCORRECT_DIR')+'/test/sample_test_kphotoz.ps', $
  pold=pold,xold=xold,yold=yold
plot, data.redshift, redshift, psym=4
k_end_print,pold=pold,xold=xold,yold=yold

end
