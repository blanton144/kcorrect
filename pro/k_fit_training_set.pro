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

flux_lmin=3500.
flux_lmax=7500.
savfile='data_k_fit_training_set.sav'

if(not file_test(savfile)) then begin
    gals=(mrdfits('training_set.test.fits',1))
;   HACK to rid of bad u-band (bad uband!)
    indx=where(gals.redshift gt 0.25,count)
    if(count gt 0) then gals[indx].mag_err[0]=5.

    maggies=10.^(-0.4*gals.mag)
    maggies_err=maggies*0.4*alog(10.)*gals.mag_err

    k_load_ascii_table,vmatrix,'vmatrix.test.dat'
    k_load_ascii_table,lambda,'lambda.test.dat'
    nv=n_elements(vmatrix)/(n_elements(lambda)-1L)

;   HACK for the moment to normalize vmatrix appropriately
    flux_indx=where(lambda ge flux_lmin and lambda lt flux_lmax)
    flux_indxp1=flux_indx+1l
    flux=dblarr(nv)
    for i=0L, nv-1L do $
      flux[i]=total((lambda[flux_indxp1]-lambda[flux_indx]) $
                    *vmatrix[flux_indx,i],/double)
    for i=0L, nv-1L do $
        vmatrix[*,i]=vmatrix[*,i]/flux[i]

;   Fit nonnegative model
    coeffs=k_fit_nonneg(maggies[*,*],maggies_err[*,*],vmatrix, $
                        lambda,redshift=gals[*].redshift, $
                        filterlist=['sdss_u0','sdss_g0','sdss_r0', $
                                    'sdss_i0','sdss_z0']+'.dat', $
                        chi2=chi2,rmatrix=rmatrix,zvals=zvals,maxiter=5000)

; rescale coeffs and vmatrix based on fit
    for i=0L,nv-1L do begin
        sigma=djsig(coeffs[i,*],sigrej=1.e+12)
        coeffs[i,*]=coeffs[i,*]/sigma
        vmatrix[*,i]=vmatrix[*,i]*sigma
    endfor
    
; Scale distribution to plane of constant flux
    flux_indx=where(lambda ge flux_lmin and lambda lt flux_lmax)
    flux_indxp1=flux_indx+1l
    flux=dblarr(nv)
    for i=0L, nv-1L do $
      flux[i]=total((lambda[flux_indxp1]-lambda[flux_indx]) $
                    *vmatrix[flux_indx,i],/double)
    flux_total=flux#coeffs
    for i=0L, nv-1L do $
      coeffs[i,*]=coeffs[i,*]/flux_total
    
    fluxed_center=dblarr(nv)
    for i=0L, nv-1L do $
      fluxed_center[i]=median(coeffs[i,*])
    
; now get the conditions for each galaxy
    condition_mean=dblarr(nv,n_elements(gals))
    condition_var=dblarr(nv,n_elements(gals))
    for i=0L, n_elements(gals)-1L do begin
        condition_mean[*,i]=fluxed_center*flux_total[i]
        condition_var[*,i]=condition_mean[10,i]*6.
    endfor
    
; now refit
    coeffs=k_fit_nonneg(maggies[*,*],maggies_err[*,*],vmatrix, $
                        lambda,redshift=gals[*].redshift, $
                        filterlist=['sdss_u0','sdss_g0','sdss_r0', $
                                    'sdss_i0','sdss_z0']+'.dat', $
                        chi2=chi2,rmatrix=rmatrix,zvals=zvals,maxiter=5000, $
                        condition_mean=condition_mean, $
                        condition_var=condition_var)

    save,filename=savfile
endif else begin
    restore,savfile
endelse

indx=lindgen(3000)
gals=gals[indx]
coeffs=coeffs[*,indx]
maggies=maggies[*,indx]
maggies_err=maggies_err[*,indx]

; Output QA 
set_print,filename='qa_k_fit_training_set_full_coeffs.ps'
k_reconstruct_maggies,coeffs,gals[*].redshift,rec_maggies, $
  zvals=zvals,rmatrix=rmatrix,ematrix=identity(nv)
umg_model=-2.5*alog10(rec_maggies[0,*]/rec_maggies[1,*])
umg_obs=gals[*].mag[0]-gals[*].mag[1]
gmr_model=-2.5*alog10(rec_maggies[1,*]/rec_maggies[2,*])
gmr_obs=gals[*].mag[1]-gals[*].mag[2]
rmi_model=-2.5*alog10(rec_maggies[2,*]/rec_maggies[3,*])
rmi_obs=gals[*].mag[2]-gals[*].mag[3]
imz_model=-2.5*alog10(rec_maggies[3,*]/rec_maggies[4,*])
imz_obs=gals[*].mag[3]-gals[*].mag[4]
plot,umg_obs,gmr_obs,psym=3,xra=[0.,3.2],yra=[-0.3,2.2]
oplot,umg_model,gmr_model,psym=4,color=djs_icolor('red')
plot,gmr_obs,rmi_obs,psym=3,xra=[-0.3,2.2],yra=[-0.1,1.0]
oplot,gmr_model,rmi_model,psym=4,color=djs_icolor('red')
plot,rmi_obs,imz_obs,psym=3,xra=[-0.1,1.0],yra=[-0.2,0.6]
oplot,rmi_model,imz_model,psym=4,color=djs_icolor('red')
plot,gals[*].redshift,umg_model-umg_obs,psym=3,yra=[-0.8,0.8]
plot,gals[*].redshift,gmr_model-gmr_obs,psym=3,yra=[-0.8,0.8]
plot,gals[*].redshift,rmi_model-rmi_obs,psym=3,yra=[-0.8,0.8]
plot,gals[*].redshift,imz_model-imz_obs,psym=3,yra=[-0.8,0.8]
end_print

; Scale coeffs to reasonable values
;median_coeffs=median(coeffs)
;coeffs=coeffs/median_coeffs

order=[0,5,10,1+lindgen(4),6+lindgen(4),11+lindgen(62)]
vmatrix=vmatrix[*,order]
coeffs=coeffs[order,*]
k_ortho_templates,vmatrix,lambda,bmatrix,bflux,bdotv=bdotv, $
  sublmin=2000.,sublmax=12000.
bcoeffs=bdotv#coeffs

; find flux direction and rotation to those coordinates
ematrix=identity(nv)
k_ortho_etemplates,ematrix,bflux

; Scale distribution to plane of constant flux
newcoeffs=transpose(ematrix)#bcoeffs
flux_total=newcoeffs[0,*]
pcacoeffs=dblarr(nv-1,n_elements(indx))
for i=0L, nv-2L do $
  pcacoeffs[i,*]=newcoeffs[i+1,*]/flux_total

; Find center of the distribution and get "centered" coeffs
center=dblarr(nv-1L)
for i=0L, nv-2L do begin & $
  tmp_mean=median(pcacoeffs[i,*]) & $
  djs_iterstat,pcacoeffs[i,*],sigrej=5.,mean=tmp_mean & $
  center[i]=tmp_mean & $
endfor
;flux_center=flux##transpose(center)
;center=center/flux_center[0]
centered_coeffs=dblarr(nv-1L,n_elements(gals))
for i=0L, nv-2L do $
  centered_coeffs[i,*]=(pcacoeffs[i,*])-center[i]

; Now trim the outliers
dist2=total(centered_coeffs^2,1,/double)
var=mean(dist2)
use_indx=where(dist2 lt 25.*var)

; Pull out the eigencomponents
npca=3
em_pca,centered_coeffs[*,use_indx],npca,tmp_eigenvec,eigencoeffs_indx, $
  maxiter=100
eigenvec=dblarr(nv,npca+1L)
eigenvec[1:nv-1,1:npca]=tmp_eigenvec
eigenvec[0,0]=1.
newematrix=dblarr(nv,npca+2L)
newematrix[*,1L:npca+1L]=ematrix#eigenvec
newematrix[*,0]=ematrix#[0.,center]
eigencoeffs=dblarr(npca+2L,n_elements(indx))
eigencoeffs[0,*]=1.
eigencoeffs[1,*]=1.
eigencoeffs[2:npca+1,*]=transpose(tmp_eigenvec)#centered_coeffs
for i=0,npca+1 do $
  eigencoeffs[i,*]=eigencoeffs[i,*]*flux_total

; Output QA 
nuse=npca-2L
k_reconstruct_maggies,eigencoeffs[0:nuse+1,*],gals[*].redshift,rec_maggies, $
  zvals=zvals,lambda=lambda,bmatrix=bmatrix,ematrix=newematrix[*,0:nuse+1], $
  filterlist=['sdss_u0','sdss_g0','sdss_r0','sdss_i0','sdss_z0']+'.dat'
plot,maggies[2,*],maggies[2,*]/rec_maggies[2,*],psym=3
set_print,filename='qa_k_fit_training_set_pca_coeffs.ps'
umg_model=-2.5*alog10(rec_maggies[0,*]/rec_maggies[1,*])
umg_obs=gals[*].mag[0]-gals[*].mag[1]
gmr_model=-2.5*alog10(rec_maggies[1,*]/rec_maggies[2,*])
gmr_obs=gals[*].mag[1]-gals[*].mag[2]
rmi_model=-2.5*alog10(rec_maggies[2,*]/rec_maggies[3,*])
rmi_obs=gals[*].mag[2]-gals[*].mag[3]
imz_model=-2.5*alog10(rec_maggies[3,*]/rec_maggies[4,*])
imz_obs=gals[*].mag[3]-gals[*].mag[4]
plot,umg_obs,gmr_obs,psym=3,xra=[0.,3.2],yra=[-0.3,2.2]
oplot,umg_model,gmr_model,psym=4,color=djs_icolor('red')
plot,gmr_obs,rmi_obs,psym=3,xra=[-0.3,2.2],yra=[-0.1,1.0]
oplot,gmr_model,rmi_model,psym=4,color=djs_icolor('red')
plot,rmi_obs,imz_obs,psym=3,xra=[-0.1,1.0],yra=[-0.2,0.6]
oplot,rmi_model,imz_model,psym=4,color=djs_icolor('red')
plot,gals[*].redshift,umg_model-umg_obs,psym=3,yra=[-0.8,0.8]
plot,gals[*].redshift,gmr_model-gmr_obs,psym=3,yra=[-0.8,0.8]
plot,gals[*].redshift,rmi_model-rmi_obs,psym=3,yra=[-0.8,0.8]
plot,gals[*].redshift,imz_model-imz_obs,psym=3,yra=[-0.8,0.8]
!P.MULTI=[0,2,2]
!X.MARGIN=2.
!Y.MARGIN=2.
plot,eigencoeffs[2,*],eigencoeffs[3,*],psym=3,xtitle='e2',ytitle='e3'
plot,eigencoeffs[2,*],eigencoeffs[4,*],psym=3,xtitle='e2',ytitle='e4'
plot,eigencoeffs[3,*],eigencoeffs[4,*],psym=3,xtitle='e3',ytitle='e4'
ecoeff=[1.,1.,-0.5,0.,0.]
spec=bmatrix#(newematrix#ecoeff)
plot,lambda,spec,xra=[1000.,13000.],yra=max(spec)*[-0.1,1.1], $
  color=djs_icolor('black')
for i=0, 9 do begin & $
    cc=-0.5+0.64*(i+1.)/10. & $
    ecoeff=[1.,1.,cc,0.,0.] & $
    spec=bmatrix#(newematrix#ecoeff) & $
    oplot,lambda,spec,color=i & $
endfor
!P.MULTI=[1,1,1]
end_print

vmatrix=bmatrix#newematrix[*,0:2]
coeffs=eigencoeffs[0:2,*]
maggies_ivar=1./maggies_err^2
k_tweak_templates, maggies, maggies_ivar, gals.redshift, $
  coeffs, vmatrix, lambda, filterlist=filterlist, $
  maggies_factor=maggies_factor, vmatrix_tweaked=vmatrix_tweaked

k_reconstruct_maggies,coeffs,gals[*].redshift,rec_maggies, $
  zvals=zvals,lambda=lambda,bmatrix=vmatrix_tweaked,ematrix=identity(3), $
  filterlist=['sdss_u0','sdss_g0','sdss_r0','sdss_i0','sdss_z0']+'.dat'
for i=0L, 4L do $
  rec_maggies[i,*]=rec_maggies[i,*]/maggies_factor[i]
set_print,filename='qa_k_fit_training_set_tweak_coeffs.ps'
umg_model=-2.5*alog10(rec_maggies[0,*]/rec_maggies[1,*])
umg_obs=gals[*].mag[0]-gals[*].mag[1]
gmr_model=-2.5*alog10(rec_maggies[1,*]/rec_maggies[2,*])
gmr_obs=gals[*].mag[1]-gals[*].mag[2]
rmi_model=-2.5*alog10(rec_maggies[2,*]/rec_maggies[3,*])
rmi_obs=gals[*].mag[2]-gals[*].mag[3]
imz_model=-2.5*alog10(rec_maggies[3,*]/rec_maggies[4,*])
imz_obs=gals[*].mag[3]-gals[*].mag[4]
plot,umg_obs,gmr_obs,psym=3,xra=[0.,3.2],yra=[-0.3,2.2]
oplot,umg_model,gmr_model,psym=4,color=djs_icolor('red')
plot,gmr_obs,rmi_obs,psym=3,xra=[-0.3,2.2],yra=[-0.1,1.0]
oplot,gmr_model,rmi_model,psym=4,color=djs_icolor('red')
plot,rmi_obs,imz_obs,psym=3,xra=[-0.1,1.0],yra=[-0.2,0.6]
oplot,rmi_model,imz_model,psym=4,color=djs_icolor('red')
plot,gals[*].redshift,umg_model-umg_obs,psym=3,yra=[-0.8,0.8]
plot,gals[*].redshift,gmr_model-gmr_obs,psym=3,yra=[-0.8,0.8]
plot,gals[*].redshift,rmi_model-rmi_obs,psym=3,yra=[-0.8,0.8]
plot,gals[*].redshift,imz_model-imz_obs,psym=3,yra=[-0.8,0.8]

!P.MULTI=[0,1,2]
ecoeff=[1.,1.,-0.5]
spec=vmatrix#ecoeff
plot,lambda,spec,xra=[1000.,13000.],yra=max(spec)*[-0.1,1.1], $
  color=djs_icolor('black')
for i=0, 9 do begin & $
    cc=-0.5+0.64*(i+1.)/10. & $
    ecoeff=[1.,1.,cc] & $
    spec=vmatrix#ecoeff & $
    oplot,lambda,spec,color=i & $
endfor

ecoeff=[1.,1.,-0.5]
spec=vmatrix_tweaked#ecoeff
plot,lambda,spec,xra=[1000.,13000.],yra=max(spec)*[-0.1,1.1], $
  color=djs_icolor('black')
for i=0, 9 do begin & $
  cc=-0.5+0.64*(i+1.)/10. & $
  ecoeff=[1.,1.,cc] & $
  spec=vmatrix_tweaked#ecoeff & $
  oplot,lambda,spec,color=i & $
  endfor

plot,lambda,vmatrix_tweaked[*,0]/vmatrix[*,0],xra=[2000.,12000.],/xlog,/ylog
end_print
stop

chi2=dblarr(100,n_elements(indx))
for i=0L, 99L do begin & $
    redshift=0.+1.*(double(i)+0.5)/100. & $
    k_photoz2, maggies, 1./maggies_err^2, redshift, lambda, vmatrix_tweaked, $
      tmp_fit, tmp_chi2 & $
    chi2[i,*]=tmp_chi2 & $
  endfor

photoz=dblarr(n_elements(indx))
minchi2=dblarr(n_elements(indx))
for i=0L, n_elements(indx)-1L do begin & $
  minchi2[i]=min(chi2[*,i],iminchi2) & $
  photoz[i]=0.+1.*(double(iminchi2)+0.5)/100. & $
  endfor    
    
stop

end
;------------------------------------------------------------------------------
