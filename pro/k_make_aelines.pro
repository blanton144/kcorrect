;+
; NAME:
;   k_make_aelines
; PURPOSE:
;   Make the 3-template A+E+lines
; CALLING SEQUENCE:
; INPUTS:
; OPTIONAL INPUTS:
; KEYWORDS:
; OUTPUTS:
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   24-jan-2002  WRitten by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_make_aelines, vmatrix, lambda

if(NOT keyword_set(pegasepath)) then $
  pegasepath=getenv('DATA')+'/specmodels/PEGASE.2'
if(NOT keyword_set(nl)) then nl=5000L
if(NOT keyword_set(lmin)) then lmin=2500.
if(NOT keyword_set(lmax)) then lmax=12000.
if(NOT keyword_set(veganame)) then veganame='lcbvega.ori'

vmatrix=dblarr(nl,2)

lambda=lmin+(lmax-lmin)*dindgen(nl+1)/double(nl)
dl=lambda[1]-lambda[0]
interp_lambda=0.5*(lambda[0:nl-1]+lambda[1:nl])

pegfile=pegasepath+'/mrb_spectra.0.008.dat'  
read_peg,pegfile,peg=peg
vmatrix[*,2]=k_interp_pegase(peg,2.,nl=nl,lmin=lmin,lmax=lmax,/nocontinuum)
vmatrix[*,0]=k_interp_pegase(peg,12000.,nl=nl,lmin=lmin,lmax=lmax)

;k_read_basel,vega_lambda,vega_flux, $
;  getenv('KCORRECT_DIR')+'/data/basel/'+veganame
;vega_lambda=vega_lambda*10.d
;cspeed=2.99792d+18              ; ang per sec
;vega_flux=3.14159*4.*vega_flux*cspeed/(vega_lambda)^2 ; from "hnu" to flambda
;vega_flux=vega_flux*6.5043898d-17         ; normalize (to fit Hayes)
;vmatrix[*,1]= $
;  interpol(vega_flux,vega_lambda,interp_lambda)
vmatrix[*,1]=k_interp_pegase(peg,100.,nl=nl,lmin=lmin,lmax=lmax)

nv=n_elements(vmatrix)/(n_elements(lambda)-1L)
for i=0,nv-1 do $
  vmatrix[*,i]=vmatrix[*,i]/total(vmatrix[*,i],/double)

dust1={dusty_str, geometry:'', dust:'', structure:'', tauv:0.}
dust=replicate(dust1,1)
dust.geometry=['dusty']
dust.dust=['MW']
dust.structure=['c']
dust.tauv=[1.5]
vmatrix[*,0]=vmatrix[*,0]*exp(-witt_ext(dust[0],dust[0].tauv,lambda[0:nl-1]))

nv=2
vmatrix=vmatrix[*,0:nv-1]
k_write_ascii_table,vmatrix,'vmatrix.aelines.dat'
k_write_ascii_table,lambda,'lambda.aelines.dat'

gals=mrdfits('training_set.test.fits',1)

maggies=10.^(-0.4*gals.mag)
maggies_err=maggies*0.4*alog(10.)*gals.mag_err

indx=where(gals.redshift lt 0.05)
;coeffs=k_fit_nonneg(maggies[1:3,indx],maggies_err[1:3,indx],vmatrix, $
                    ;lambda,redshift=gals[indx].redshift, $
                    ;filterlist=['sdss_g0','sdss_r0','sdss_i0']+'.dat', $
                    ;chi2=chi2,rmatrix=rmatrix,zvals=zvals,maxiter=5000)
k_fit_coeffs,maggies[1:3,indx],1./maggies_err[1:3,indx]^2, $
  gals[indx].redshift,coeffs2,bmatrix=vmatrix,lambda=lambda, $
  filterlist=['sdss_g0','sdss_r0','sdss_i0']+'.dat', $
  ematrix=identity(nv)

;k_reconstruct_maggies,coeffs2,gals[indx].redshift,rec_maggies, $
;  zvals=zvals,rmatrix=rmatrix,ematrix=identity(nv)

k_reconstruct_maggies,coeffs2,gals[indx].redshift,rec_maggies, $
  bmatrix=vmatrix,lambda=lambda, $
  filterlist=['sdss_g0','sdss_r0','sdss_i0']+'.dat', $
  ematrix=identity(nv)

;umg_model=-2.5*alog10(rec_maggies[0,*]/rec_maggies[1,*])
gmr_model=-2.5*alog10(rec_maggies[0,*]/rec_maggies[1,*])
rmi_model=-2.5*alog10(rec_maggies[1,*]/rec_maggies[2,*])
;imz_model=-2.5*alog10(rec_maggies[3,*]/rec_maggies[4,*])
;umg_obs=gals[indx].mag[0]-gals[indx].mag[1]
gmr_obs=gals[indx].mag[1]-gals[indx].mag[2]
rmi_obs=gals[indx].mag[2]-gals[indx].mag[3]
;imz_obs=gals[indx].mag[3]-gals[indx].mag[4]

;window,0
;plot,umg_obs,gmr_obs,psym=4,xra=[0.,3.],yra=[0.,1.0]
;oplot,umg_model,gmr_model,psym=4,color=255

window,0
plot,gmr_obs,rmi_obs,psym=4,xra=[0.,1.],yra=[0.,0.6]
oplot,gmr_model,rmi_model,psym=4,color=255

;window,2
;plot,rmi_obs,imz_obs,psym=4,xra=[0.,0.6],yra=[0.,0.6]
;oplot,rmi_model,imz_model,psym=4,color=255

nv=2
k_load_ascii_table,vmatrix,'vmatrix.aelines.dat'
k_load_ascii_table,lambda,'lambda.aelines.dat'

im_g=mrdfits('image9400070101970152.fit',1)
im_r=mrdfits('image9400070101970152.fit',2)
im_i=mrdfits('image9400070101970152.fit',3)

tsfield=mrdfits('tsField.fit.gz',1)

maggies=fltarr(3,n_elements(im_g))
maggies_err=fltarr(3,n_elements(im_g))+1.
maggies[0,*]=im_g/55.247* $
  10.^(0.4*(tsfield.aa[1]+tsfield.airmass[1]*tsfield.kk[1]))- $
  10.^(-0.4*tsfield.sky[1])
maggies[1,*]=im_r/55.247* $
  10.^(0.4*(tsfield.aa[2]+tsfield.airmass[2]*tsfield.kk[1]))- $
  10.^(-0.4*tsfield.sky[2])
maggies[2,*]=im_i/55.247* $
  10.^(0.4*(tsfield.aa[3]+tsfield.airmass[3]*tsfield.kk[1]))- $
  10.^(-0.4*tsfield.sky[3])

k_fit_coeffs,maggies[0:2,*],1./maggies_err[0:2,*]^2, $
  replicate(0.024,n_elements(im_g)),coeffs2,bmatrix=vmatrix,lambda=lambda, $
  filterlist=['sdss_g0','sdss_r0','sdss_i0']+'.dat', $
  ematrix=identity(nv)

im0=reform(coeffs2[0,*],200,200)
im1=reform(coeffs2[1,*],200,200)

stop

end
;------------------------------------------------------------------------------
