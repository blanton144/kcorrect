pro test_kcorrect,failed

if(NOT file_test(getenv('KCORRECT_DIR')+'/test/sample.fits')) then begin
    failed=1
    return
endif

data=mrdfits(getenv('KCORRECT_DIR')+'/test/sample.fits',1)
kcorrect, data.maggies, data.maggies_ivar, data.redshift, kcorrect, $
  filterlist=['sdss_u0.par','sdss_g0.par','sdss_r0.par','sdss_i0.par', $
              'sdss_z0.par'], coeffs=coeffs, maxiter=10000
k_print,filename=getenv('KCORRECT_DIR')+'/test/sample_test_kcorrect.ps', $
  pold=pold,xold=xold,yold=yold
for i=0, n_elements(data[0].maggies)-1L do $
  plot, data.redshift, kcorrect[i,*], psym=4
k_end_print,pold=pold,xold=xold,yold=yold

end
