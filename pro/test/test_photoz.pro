pro test_photoz

data=mrdfits('sdss_training_set.test.fits',1)
print,systime()
kphotoz,data.maggies[0:4],data.maggies_ivar[0:4],photoz, $
  vfile='vmatrix.smooth3.dat', lfile='lambda.smooth3.dat', chi2=chi2, $
  /verbose
print,systime()
splot,data.redshift,photoz,psym=4

i=1000
i=ii[420]
kphotoz,data[i].maggies[0:4],data[i].maggies_ivar[0:4],photoz, $
  vfile='vmatrix.smooth3.dat', lfile='lambda.smooth3.dat', chi2=chi2, $
  /verbose
test_redshift=0.001*(dindgen(1000)+0.5) 
test_maggies=data[i].maggies[0:4]#replicate(1.,1000) 
test_maggies_ivar=data[i].maggies_ivar[0:4]#replicate(1.,1000) 
kcorrect,test_maggies,test_maggies_ivar,test_redshift,kcorrect, $ 
  chi2=chi2v, vfile='vmatrix.smooth3.dat', lfile='lambda.smooth3.dat' 
splot,test_redshift,chi2v

data=mrdfits('sdss_training_set.test.fits',1)
indx=shuffle_indx(n_elements(data),num_sub=1)
indx=lindgen(n_elements(data))
kphotoz, data[indx].maggies[0:4], data[indx].maggies_ivar[0:4], photoz_3, $
  vfile='vmatrix.test3.dat', lfile='lambda.test3.dat',chi2=chi2
kphotoz, data[indx].maggies[0:4], data[indx].maggies_ivar[0:4], photoz_3_a, $
  vfile='vmatrix.test3.dat', lfile='lambda.test3.dat',chi2=chi2, $
  lprior=1.*alog((lindgen(1000)+0.5)/1000.)
indx2=where(photoz_3 lt 0.03 and data[indx].redshift gt 0.1)
indx2=lindgen(50)
set_print,filename='blah.ps'
!P.MULTI=[0,1,2]
for i=0, n_elements(indx2)-1L do begin & $
test_redshift=0.001*(dindgen(1000)+0.5) & $
test_maggies=data[indx[indx2[i]]].maggies[0:4]#replicate(1.,1000) & $
test_maggies_ivar=data[indx[indx2[i]]].maggies_ivar[0:4]#replicate(1.,1000) & $
kcorrect,test_maggies,test_maggies_ivar,test_redshift,kcorrect, $ 
  chi2=chi2v, vfile='vmatrix.test3.dat', lfile='lambda.test3.dat' & $
plot,test_redshift,chi2v,xra=[0.,0.5] & $
oplot,[data[indx[indx2[i]]].redshift]#replicate(1.,2),[0.,1000.] & $
plot,test_redshift,chi2v-2.*alog(test_redshift),xra=[0.,0.5] & $
oplot,[data[indx[indx2[i]]].redshift]#replicate(1.,2),[0.,1000.] & $
print,data[indx[indx2[i]]].redshift, chi2[indx2[i]] & $
  endfor
end_print

kphotoz, data[indx].maggies[0:4], data[indx].maggies_ivar[0:4], photoz_std
curr_photoz=photoz_std

plot,data[indx].redshift,curr_photoz,psym=4
oplot,[0.,1.],[0.,1]
diff=curr_photoz-data[indx].redshift
plot,-2.5*alog10(data[indx].maggies[2]),diff,psym=4

kphotoz, data[indx].maggies[0:4], data[indx].maggies_ivar[0:4], photoz_2, $
  vfile='vmatrix.test2.dat', lfile='lambda.test2.dat'
curr_photoz=photoz_2

plot,data[indx].redshift,curr_photoz,psym=4
oplot,[0.,1.],[0.,1]
diff=curr_photoz-data[indx].redshift
plot,-2.5*alog10(data[indx].maggies[2]),diff,psym=4

kphotoz, data[indx].maggies[0:4], data[indx].maggies_ivar[0:4], photoz_3, $
  vfile='vmatrix.test3.dat', lfile='lambda.test3.dat'
curr_photoz=photoz_3

plot,data[indx].redshift,curr_photoz,psym=4
oplot,[0.,1.],[0.,1]
diff=curr_photoz-data[indx].redshift
plot,-2.5*alog10(data[indx].maggies[2]),diff,psym=4

kphotoz, data[indx].maggies[0:4], data[indx].maggies_ivar[0:4], photoz_4, $
  vfile='vmatrix.test4.dat', lfile='lambda.test4.dat'
curr_photoz=photoz_4

plot,data[indx].redshift,curr_photoz,psym=4
oplot,[0.,1.],[0.,1]
diff=curr_photoz-data[indx].redshift
plot,-2.5*alog10(data[indx].maggies[2]),diff,psym=4

kphotoz, data[indx].maggies[0:4], data[indx].maggies_ivar[0:4], photoz_5, $
  vfile='vmatrix.test5.dat', lfile='lambda.test5.dat'
curr_photoz=photoz_5

plot,data[indx].redshift,curr_photoz,psym=4
oplot,[0.,1.],[0.,1]
diff=curr_photoz-data[indx].redshift
plot,-2.5*alog10(data[indx].maggies[2]),diff,psym=4

diff=photoz_2-data[indx].redshift & print, djsig(diff,sigrej=100)
diff=photoz_3-data[indx].redshift & print, djsig(diff,sigrej=100)
diff=photoz_4-data[indx].redshift & print, djsig(diff,sigrej=100)
diff=photoz_5-data[indx].redshift & print, djsig(diff,sigrej=100)

end
