; make the code test sample
pro mk_test_sample

data=mrdfits(getenv('DATA')+'/sdss/kcorrect/data/v3_1/'+ $
             'sdss_training_set.test.fits',1)
indx=shuffle_indx(n_elements(data),num_sub=1000)
data=data[indx]
out1={maggies:fltarr(5), $
      maggies_ivar:fltarr(5), $
      redshift:0.}
out=replicate(out1,n_elements(data))
out.maggies=data.maggies[0:4]
out.maggies_ivar=data.maggies_ivar[0:4]
out.redshift=data.redshift

mwrfits,out,getenv('KCORRECT_DIR')+'/test/sample.fits',/create
openw,unit,getenv('KCORRECT_DIR')+'/test/sample.dat',/get_lun
for i=0L, n_elements(out)-1L do $
  printf,unit,format='(%"%e %e %e %e %e %e %e %e %e %e %e")', $
  out[i].redshift,out[i].maggies,out[i].maggies_ivar
free_lun,unit

end
