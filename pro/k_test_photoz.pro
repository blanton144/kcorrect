pro k_test_photoz

a=hogg_mrdfits('/data/sdss/spectro/spAll.fits',1,nrowchunk=3000, $
               columns=['counts_model','counts_modelerr','reddening','z', $
                        'class','primtarget','objc_flags','objc_flags2', $
                        'parent','nchild'])
indx=where(a.class eq 'GALAXY' and a.primtarget gt 0)
help,indx
model=a[indx].counts_model-a[indx].reddening
modelerr=a[indx].counts_modelerr
kphotoz,model,modelerr,photoz,coeff=coeff,/sdssfix
window,retain=2
plot,a[indx].z,photoz,psym=3,xra=[0.,0.5]
save,a,indx,photoz,model,modelerr,coeff,filename='k_test_photoz.sav'

end
