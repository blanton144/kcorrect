pro test_kphotoz_standalone,failed

; run reconstruction of maggies at the current redshift and 
; at redshift zero
testfile=getenv('KCORRECT_DIR')+'/test/sample.dat'
cmd='cat '+testfile+' | awk '+"'"+ $
  '{printf("%e %e %e %e %e %e %e %e %e %e\n",'+ $
  '$2,$3,$4,$5,$6,$7,$8,$9,$10,$11)}'+"'"+ $
  ' | fit_photoz'+' >! '+getenv('KCORRECT_DIR')+'/test/photoz.dat'
print,cmd
spawn,cmd

; read in the results
ngals=numlines(getenv('KCORRECT_DIR')+'/test/sample.dat')
inarr=fltarr(11,ngals)
openr,unit,getenv('KCORRECT_DIR')+'/test/sample.dat',/get_lun
readf,unit,inarr
free_lun,unit
redshift=inarr[0,*]
maggies=inarr[1:5,*]
maggies_ivar=inarr[6:10,*]
inarr=0

openr,unit,getenv('KCORRECT_DIR')+'/test/photoz.dat',/get_lun
inarr=fltarr(4,ngals)
readf,unit,inarr
free_lun,unit
photoz=inarr[0,*]
coeffs=inarr[1:3,*]
inarr=0

k_print,filename=getenv('KCORRECT_DIR')+ $
  '/test/sample_test_kphotoz_standalone.ps', $
  pold=pold,xold=xold,yold=yold
plot, redshift, photoz, psym=4
k_end_print,pold=pold,xold=xold,yold=yold

end
