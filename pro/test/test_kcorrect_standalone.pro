pro test_kcorrect_standalone,failed

; run reconstruction of maggies at the current redshift and 
; at redshift zero
testfile=getenv('KCORRECT_DIR')+'/test/sample.dat'
cmd='cat '+testfile+' | fit_coeffs '+ $
  ' >! '+getenv('KCORRECT_DIR')+'/test/coeffs.dat'
print,cmd
spawn,cmd
cmd='cat '+testfile+' | fit_coeffs | reconstruct_maggies'+ $
  ' >! '+getenv('KCORRECT_DIR')+'/test/maggies.dat'
print,cmd
spawn,cmd
cmd='cat '+testfile+' | fit_coeffs | reconstruct_maggies --redshift 0.'+ $
  ' >! '+getenv('KCORRECT_DIR')+'/test/maggies.z0.dat'
print,cmd
spawn,cmd

; read in the results
ngals=numlines(getenv('KCORRECT_DIR')+'/test/maggies.dat')
openr,unit,getenv('KCORRECT_DIR')+'/test/maggies.dat',/get_lun
inarr=fltarr(6,ngals)
readf,unit,inarr
free_lun,unit
redshift=inarr[0,*]
maggies=inarr[1:5,*]
inarr=0

openr,unit,getenv('KCORRECT_DIR')+'/test/maggies.z0.dat',/get_lun
inarr=fltarr(6,ngals)
readf,unit,inarr
free_lun,unit
maggies_z0=inarr[1:5,*]
inarr=0

; construct kcorrections
kcorrect=-2.5*alog10(maggies/maggies_z0)

k_print,filename=getenv('KCORRECT_DIR')+ $
  '/test/sample_test_kcorrect_standalone.ps', $
  pold=pold,xold=xold,yold=yold
for i=0, n_elements(maggies[*,0])-1L do $
  plot, redshift, kcorrect[i,*], psym=4
k_end_print,pold=pold,xold=xold,yold=yold

end
