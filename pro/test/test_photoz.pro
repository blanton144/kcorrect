pro test_photoz

data=mrdfits('sdss_training_set.test.fits',1)
indx=shuffle_indx(n_elements(data),num_sub=300)
indx=lindgen(n_elements(data))

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
