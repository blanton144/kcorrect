;+
; NAME:
;   k_fit_training_set
; PURPOSE:
;   Use the training set to develop a set of templates
; CALLING SEQUENCE:
;   k_fit_training_set
; INPUTS:
; OPTIONAL INPUTS:
; OUTPUTS:
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   18-Jan-2003  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_fit_training_set

gals=mrdfits('training_set.test.fits',1)

maggies=10.^(-0.4*gals.mag)
maggies_err=maggies*0.4*alog(10.)*gals.mag_err

k_load_ascii_table,vmatrix,'vmatrix.test.dat'
k_load_ascii_table,lambda,'lambda.test.dat'
nv=n_elements(vmatrix)/(n_elements(lambda)-1L)

indx=where(gals.mag[0,*]-gals.mag[1,*] lt 1. and $
           gals.mag[1,*]-gals.mag[2,*] gt 0.5 and $
           gals.redshift gt 0.25)
indx=lindgen(n_elements(gals))
coeffs=k_fit_nonneg(maggies[*,indx],maggies_err[*,indx],vmatrix, $
                    lambda,redshift=gals[indx].redshift,version='test', $
                    vpath='.',filterlist=['sdss_u0','sdss_g0','sdss_r0', $
                                          'sdss_i0','sdss_z0']+'.dat', $
                    chi2=chi2,rmatrix=rmatrix,zvals=zvals,maxiter=5000)

k_reconstruct_maggies,coeffs,gals[indx].redshift,rec_maggies, $
  zvals=zvals,rmatrix=rmatrix,ematrix=identity(nv)

umg_model=-2.5*alog10(rec_maggies[0,*]/rec_maggies[1,*])
umg_obs=gals[indx].mag[0]-gals[indx].mag[1]
gmr_model=-2.5*alog10(rec_maggies[1,*]/rec_maggies[2,*])
gmr_obs=gals[indx].mag[1]-gals[indx].mag[2]
rmi_model=-2.5*alog10(rec_maggies[2,*]/rec_maggies[3,*])
rmi_obs=gals[indx].mag[2]-gals[indx].mag[3]
imz_model=-2.5*alog10(rec_maggies[3,*]/rec_maggies[4,*])
imz_obs=gals[indx].mag[3]-gals[indx].mag[4]

set_print,filename='k_fit_training_set_restricted.ps'
plot,umg_obs,gmr_obs,psym=3,xra=[0.,3.2],yra=[-0.3,2.2]
oplot,umg_model,gmr_model,psym=4,color=djs_icolor('red')

plot,gmr_obs,rmi_obs,psym=3,xra=[-0.3,2.2],yra=[-0.1,1.0]
oplot,gmr_model,rmi_model,psym=4,color=djs_icolor('red')

plot,rmi_obs,imz_obs,psym=3,xra=[-0.1,1.0],yra=[-0.2,0.6]
oplot,rmi_model,imz_model,psym=4,color=djs_icolor('red')

plot,gals[indx].redshift,umg_model-umg_obs,psym=3,yra=[-0.4,0.4]
plot,gals[indx].redshift,gmr_model-gmr_obs,psym=3,yra=[-0.4,0.4]
plot,gals[indx].redshift,rmi_model-rmi_obs,psym=3,yra=[-0.4,0.4]
plot,gals[indx].redshift,imz_model-imz_obs,psym=3,yra=[-0.4,0.4]

end_print

; get median
indx=where(total(coeffs[*,0:2000],1,/double) lt 1.D-30)
k_reconstruct_maggies,coeffs[*,indx],replicate(0.,n_elements(indx)), $
  rec_maggies0,zvals=zvals,rmatrix=rmatrix,ematrix=identity(nv), $
  band_shift=replicate(0.1,n_elements(indx))
;ranscale=(0.995+0.1*randomu(seed,n_elements(indx)))
for i=0L, nv-1L do $
  coeffs[i,indx]=coeffs[i,indx]/rec_maggies0[2,*]

coeffs=coeffs*1.D+45
center=dblarr(nv)
for i=0L, nv-1L do $
  center[i]=djs_avsigclip(coeffs[i,indx],sigrej=3.)
centered_coeffs=dblarr(nv,n_elements(indx))
for i=0L, nv-1L do $
  centered_coeffs[i,indx]=(coeffs[i,indx])-center[i]

nv=73
em_pca,centered_coeffs,nv,eigenvec,newcoeffs
ematrix=eigenvec
scale=sqrt(total(ematrix^2,1,/double))
for i=0,nv-1 do ematrix[*,i]=ematrix[*,i]/scale[i]
;newcoeffs=(centered_coeffs)##invert(ematrix)
newcenter=center#(ematrix)
for i=0,nv-1 do newcoeffs[i,*]=newcoeffs[i,*]/scale[i]+newcenter[i]
k_reconstruct_maggies,newcoeffs,gals[indx].redshift,rec_maggies, $
  zvals=zvals,rmatrix=rmatrix,ematrix=ematrix

umg_model=-2.5*alog10(rec_maggies[0,*]/rec_maggies[1,*])
umg_obs=gals[indx].mag[0]-gals[indx].mag[1]
gmr_model=-2.5*alog10(rec_maggies[1,*]/rec_maggies[2,*])
gmr_obs=gals[indx].mag[1]-gals[indx].mag[2]
rmi_model=-2.5*alog10(rec_maggies[2,*]/rec_maggies[3,*])
rmi_obs=gals[indx].mag[2]-gals[indx].mag[3]
imz_model=-2.5*alog10(rec_maggies[3,*]/rec_maggies[4,*])
imz_obs=gals[indx].mag[3]-gals[indx].mag[4]

set_print,filename='k_fit_training_set_restricted.ps'
plot,umg_obs,gmr_obs,psym=3,xra=[0.,3.2],yra=[-0.3,2.2]
oplot,umg_model,gmr_model,psym=4,color=djs_icolor('red')

plot,gmr_obs,rmi_obs,psym=3,xra=[-0.3,2.2],yra=[-0.1,1.0]
oplot,gmr_model,rmi_model,psym=4,color=djs_icolor('red')

plot,rmi_obs,imz_obs,psym=3,xra=[-0.1,1.0],yra=[-0.2,0.6]
oplot,rmi_model,imz_model,psym=4,color=djs_icolor('red')

plot,gals[indx].redshift,umg_model-umg_obs,psym=3,yra=[-0.4,0.4]
plot,gals[indx].redshift,gmr_model-gmr_obs,psym=3,yra=[-0.4,0.4]
plot,gals[indx].redshift,rmi_model-rmi_obs,psym=3,yra=[-0.4,0.4]
plot,gals[indx].redshift,imz_model-imz_obs,psym=3,yra=[-0.4,0.4]

end_print


;median_coeffs=dblarr(nv)
;for i=0L, nv-1L do median_coeffs[i]=median(coeffs[i,*])

; refit?

; PCA 
indx=lindgen(2000)
coeffs=coeffs*1.D+45
coeffs=(coeffs > 1.D-1) < 1.D+1
num=2000 & nv=70
indx=lindgen(num) & coeffs=double(exp(1.*randomu(seed,nv,num)) )
for i=0L, nv-1L do $
  centered_coeffs[i,indx]=alog10(coeffs[i,indx])-center[i]
;newcoeffs=pcomp(centered_coeffs,coefficients=ematrix,/double, $
;                eigenvalues=eigenvalues, /covariance, $
;                variances=var)
pca,transpose(centered_coeffs),eigenvalues,ematrix,percent,newcoeffs, $
  blah,matrix=matrix, /silent, /covariance
scale=sqrt(total(ematrix^2,1,/double))
for i=0,nv-1 do ematrix[*,i]=ematrix[*,i]/scale[i]
;newcoeffs=(centered_coeffs)##invert(ematrix)
newcenter=center##(ematrix)
for i=0,nv-1 do newcoeffs[i,*]=newcoeffs[i,*]/scale[i]+newcenter[i]
print,coeffs[*,0]-10.^(newcoeffs[*,0]##invert(ematrix)),minmax(coeffs)
blah=newcoeffs##invert(ematrix)
plot,blah[0,*],blah[1,*],psym=4
oplot,alog10(coeffs[0,*]),alog10(coeffs[1,*]),psym=4,color=255
for i=1, 9 do begin & window,i-1 & plot,newcoeffs[0,*],newcoeffs[i,*],psym=4 & wait,0.2 & endfor
nuse=nv-1
k_reconstruct_maggies,newcoeffs[0:nuse-1,*],gals[indx].redshift,rec_maggies, $
  zvals=zvals,rmatrix=rmatrix,ematrix=(ematrix[*,0:nuse-1])

umg_model=-2.5*alog10(rec_maggies[0,*]/rec_maggies[1,*])
umg_obs=gals[indx].mag[0]-gals[indx].mag[1]
gmr_model=-2.5*alog10(rec_maggies[1,*]/rec_maggies[2,*])
gmr_obs=gals[indx].mag[1]-gals[indx].mag[2]
rmi_model=-2.5*alog10(rec_maggies[2,*]/rec_maggies[3,*])
rmi_obs=gals[indx].mag[2]-gals[indx].mag[3]
imz_model=-2.5*alog10(rec_maggies[3,*]/rec_maggies[4,*])
imz_obs=gals[indx].mag[3]-gals[indx].mag[4]

set_print,filename='k_fit_training_set_restricted.ps'
plot,umg_obs,gmr_obs,psym=3,xra=[0.,3.2],yra=[-0.3,2.2]
oplot,umg_model,gmr_model,psym=4,color=djs_icolor('red')

plot,gmr_obs,rmi_obs,psym=3,xra=[-0.3,2.2],yra=[-0.1,1.0]
oplot,gmr_model,rmi_model,psym=4,color=djs_icolor('red')

plot,rmi_obs,imz_obs,psym=3,xra=[-0.1,1.0],yra=[-0.2,0.6]
oplot,rmi_model,imz_model,psym=4,color=djs_icolor('red')

plot,gals[indx].redshift,umg_model-umg_obs,psym=3,yra=[-0.4,0.4]
plot,gals[indx].redshift,gmr_model-gmr_obs,psym=3,yra=[-0.4,0.4]
plot,gals[indx].redshift,rmi_model-rmi_obs,psym=3,yra=[-0.4,0.4]
plot,gals[indx].redshift,imz_model-imz_obs,psym=3,yra=[-0.4,0.4]

end_print

stop

end
;------------------------------------------------------------------------------
