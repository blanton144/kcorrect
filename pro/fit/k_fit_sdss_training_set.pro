;+
; NAME:
;   k_fit_sdss_training_set
; PURPOSE:
;   Use the SDSS training set to develop a set of templates
; CALLING SEQUENCE:
;   k_fit_sdss_training_set
; INPUTS:
; OPTIONAL INPUTS:
; OUTPUTS:
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
;   ToDo -- track units better
;           map coefficients BACK to SFH. 
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   18-Jan-2003  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_fit_sdss_training_set, name=name

if(NOT keyword_set(name)) then name='test'
if(NOT keyword_set(npca)) then npca=5
if(NOT keyword_set(sublmin)) then sublmin=2000.
if(NOT keyword_set(sublmax)) then sublmax=12000.

savfile='data_k_fit_sdss_training_set.sav'

if(not file_test(savfile)) then begin
    gals=(mrdfits('sdss_training_set.'+name+'.fits',1))

;   HACK to improve precision
    err_band=[0.034,0.017,0.017,0.017,0.027]
    for i=0, 4 do begin
        newvar=1./gals.mag_ivar[i]-err_band[i]^2
        gals.mag_ivar[i]=1./newvar
    endfor

;   HACK to rid of bad u-band (bad uband!)
    indx=where(gals.mag[0]-gals.mag[1] lt 1.6 and $
               gals.mag[1]-gals.mag[2] gt 1.4,count)
    if(count gt 0) then gals[indx].mag_ivar[0]=0.

    maggies=10.^(-0.4*gals.mag)
    maggies_ivar=gals.mag_ivar/(maggies*0.4*alog(10.))^2

    k_read_ascii_table,vmatrix,'vmatrix.'+name+'.dat'
    k_read_ascii_table,lambda,'lambda.'+name+'.dat'
    nv=n_elements(vmatrix)/(n_elements(lambda)-1L)

;   HACK for reasonable numbers 
    vmatrix=vmatrix/1.e+38
    maggies=maggies*1.e+9
    maggies_ivar=maggies_ivar/1.e+18

;   fit nonnegative model
    coeffs=k_fit_nonneg(maggies[*,0:4999],maggies_ivar[*,0:4999],vmatrix, $
                        lambda,redshift=gals[0:4999].redshift, $
                        filterlist=['sdss_u0.par','sdss_g0.par', $
                                    'sdss_r0.par','sdss_i0.par', $
                                    'sdss_z0.par'], $
                        chi2=chi2,rmatrix=rmatrix,zvals=zvals,maxiter=10000)

    save,filename=savfile
endif else begin
    splog,'restoring'
    restore,savfile
    splog,'done'
endelse

if (0) then begin
gals_orig=gals
indx=lindgen(5000)
gals=gals[indx]
coeffs=coeffs[*,indx]
maggies_orig=maggies
maggies=maggies[*,indx]
maggies_ivar_orig=maggies_ivar
maggies_ivar=maggies_ivar[*,indx]
nv=n_elements(coeffs)/n_elements(gals)

;specs=vmatrix#coeffs
;filterlist=['sdss_u0.par','sdss_g0.par','sdss_r0.par','sdss_i0.par', $
;            'sdss_z0.par']
;k_projection_table, rmatrix, specs, lambda, zvals, filterlist
;
;igal=long(5000.*randomu(seed,1))
;chi2=dblarr(1000)
;for i=0, 999 do begin & $
;    k_reconstruct_maggies,1.,gals_orig[igal].redshift, $
;      recmaggies,rmatrix=rmatrix[*,i,*],zvals=zvals & $
;    chi2[i]=total((recmaggies-maggies_orig[*,igal])^2*maggies_ivar[*,igal]) & $
;endfor
;plothist,alog10(chi2),bin=0.1
;
;stop

;   for kicks, reconstruct maggies
compare_maggies,coeffs,gals.redshift,maggies,maggies_ivar, $
  rmatrix=rmatrix,zvals=zvals,recmaggies=recmaggies, $
  filename='compare_maggies.ps'

; find flux direction, and scale
k_ortho_templates,coeffs,vmatrix,lambda,bcoeffs,bmatrix,bflux,bdotv=bdotv, $
  bdotb=bdotb,sublmin=sublmin,sublmax=sublmax

; Scale distribution to plane of constant flux
flux_total=bcoeffs[0,*]
pcacoeffs=fltarr(nv-1,n_elements(indx))
for i=0L, nv-2L do $
  pcacoeffs[i,*]=bcoeffs[i+1,*]/flux_total

; Find center of the distribution and get "centered" coeffs
center=fltarr(nv-1L)
for i=0L, nv-2L do begin 
    djs_iterstat,pcacoeffs[i,*],sigrej=5.,mean=tmp_mean 
    center[i]=tmp_mean 
endfor
centcoeffs=fltarr(nv-1L,n_elements(gals))
for i=0L, nv-2L do $
  centcoeffs[i,*]=(pcacoeffs[i,*])-center[i]

; Now trim the outliers
dist2=total(centcoeffs^2,1,/double)
var=mean(dist2)
use_indx=where(dist2 lt 25.*var)

; PCA the coefficients
em_pca,centcoeffs[*,use_indx],npca,tmp_eigenvec,eigencoeffs_indx,maxiter=100

; Interpret the results
eigenvec=fltarr(nv,npca+1L)
eigenvec[1:nv-1,1:npca]=tmp_eigenvec
eigenvec[0,0]=1.
eigenvec[1:nv-1,0]=eigenvec[1:nv-1,0]+center
eigenmatrix=bmatrix#eigenvec
eigencoeffs=fltarr(npca+1L,n_elements(indx))
eigencoeffs[0,*]=1.
eigencoeffs[1:npca,*]=transpose(tmp_eigenvec)#centcoeffs
for i=0,npca do $
  eigencoeffs[i,*]=eigencoeffs[i,*]*flux_total

data=fltarr(npca,n_elements(gals))
for i=0,npca-1 do $
  data[i,*]=eigencoeffs[i+1,*]/eigencoeffs[0,*]
kmeans_streams, 10, data, group, group_mean=group_mean

loadct,0
hogg_manyd_scatterplot,fltarr(n_elements(gals))+1., $
  data,'eigencoeffs_scaled.ps',exponent=0.25,xnpix=30,ynpix=30, $
  range=range
for i=0, 9 do $
  hogg_manyd_scatterplot,fltarr(n_elements(where(group eq i)))+1., $
  data[*,where(group eq i)],'eigencoeffs_scaled'+strtrim(string(i),2)+'.ps', $
  exponent=0.25,xnpix=30,ynpix=30,range=range

endif

if(0) then begin
l=0
chi2=dblarr(10,10,10)+1.d+30
mean=dblarr(10,10,10)+1.d+30
median=dblarr(10,10,10)+1.d+30
for i=0, 9 do begin
    for j=i+1, 9 do begin
        for k=j+1, 9 do begin
            indx0=where(group eq i)
            indx1=where(group eq j)
            indx2=where(group eq k)
            avgeigencoeffs=fltarr(npca+1,3)
            avgeigencoeffs[0,*]=1.
            avgeigencoeffs[1:npca,0]= $
              total(eigencoeffs[1:npca,indx0]/ $
                    (replicate(1.,npca)#eigencoeffs[0,indx0]),2)$
              /double(n_elements(indx0))
            avgeigencoeffs[1:npca,1]= $
              total(eigencoeffs[1:npca,indx1]/ $
                    (replicate(1.,npca)#eigencoeffs[0,indx1]),2)$
              /double(n_elements(indx1))
            avgeigencoeffs[1:npca,2]= $
              total(eigencoeffs[1:npca,indx2]/ $
                    (replicate(1.,npca)#eigencoeffs[0,indx2]),2)$
              /double(n_elements(indx2))
            avgspec=fltarr(n_elements(lambda)-1L,3)
            avgspec[*,0]=eigenmatrix#avgeigencoeffs[*,0]
            avgspec[*,1]=eigenmatrix#avgeigencoeffs[*,1]
            avgspec[*,2]=eigenmatrix#avgeigencoeffs[*,2]
            print,i,j,k
            set_plot,'x'
            !P.MULTI=[3,1,3]
            plot,lambda,avgspec[*,0],xra=[2000.,12000.]
            plot,lambda,avgspec[*,1],xra=[2000.,12000.]
            plot,lambda,avgspec[*,2],xra=[2000.,12000.]
            chi2gals=0
            rmatrix=0
            avgcoeffs=k_fit_nonneg(maggies,maggies_ivar,avgspec, $
                                   lambda,redshift=gals.redshift, $
                                   filterlist=['sdss_u0.par','sdss_g0.par', $
                                               'sdss_r0.par','sdss_i0.par', $
                                               'sdss_z0.par'], $
                                   chi2=chi2gals,rmatrix=rmatrix,zvals=zvals, $
                                   maxiter=10000)
            djs_iterstat,chi2gals,sigrej=3,mean=tmp_mean,median=tmp_median
            chi2[i,j,k]=total(chi2gals,/double)
            mean[i,j,k]=tmp_mean
            median[i,j,k]=tmp_median
            print,i,j,k,chi2[i,j,k],mean[i,j,k],median[i,j,k]
        endfor
    endfor
endfor

minmean=min(mean,minl)

minl=820
k=minl / 100L
j=(minl / 10L) mod 10L
i=minl mod 10L

indx0=where(group eq i)
indx1=where(group eq j)
indx2=where(group eq k)
avgeigencoeffs=fltarr(npca+1,3)
avgeigencoeffs[0,*]=1.
avgeigencoeffs[1:npca,0]= $
  total(eigencoeffs[1:npca,indx0]/ $
        (replicate(1.,npca)#eigencoeffs[0,indx0]),2)$
  /double(n_elements(indx0))
avgeigencoeffs[1:npca,1]= $
  total(eigencoeffs[1:npca,indx1]/ $
        (replicate(1.,npca)#eigencoeffs[0,indx1]),2)$
  /double(n_elements(indx1))
avgeigencoeffs[1:npca,2]= $
  total(eigencoeffs[1:npca,indx2]/ $
        (replicate(1.,npca)#eigencoeffs[0,indx2]),2)$
  /double(n_elements(indx2))
avgspec=fltarr(n_elements(lambda)-1L,3)
avgspec[*,0]=eigenmatrix#avgeigencoeffs[*,0]
avgspec[*,1]=eigenmatrix#avgeigencoeffs[*,1]
avgspec[*,2]=eigenmatrix#avgeigencoeffs[*,2]
print,i,j,k
set_plot,'x'
!P.MULTI=[3,1,3]
plot,lambda,avgspec[*,0],xra=[2000.,12000.]
plot,lambda,avgspec[*,1],xra=[2000.,12000.]
plot,lambda,avgspec[*,2],xra=[2000.,12000.]
chi2gals=0
rmatrix=0
avgcoeffs=k_fit_nonneg(maggies,maggies_ivar,avgspec, $
                       lambda,redshift=gals.redshift, $
                       filterlist=['sdss_u0.par','sdss_g0.par', $
                                   'sdss_r0.par','sdss_i0.par', $
                                   'sdss_z0.par'], $
                       chi2=chi2gals,rmatrix=rmatrix,zvals=zvals, $
                       maxiter=10000)

compare_maggies,avgcoeffs,gals.redshift,maggies,maggies_ivar, $
  lambda=lambda,vmatrix=avgspec, $
  filterlist=['sdss_u0.par','sdss_g0.par','sdss_r0.par','sdss_i0.par', $
              'sdss_z0.par'], $
  recmaggies=recmaggies, filename='compare_maggies_avgspec.ps'

goodindx=where(chi2gals lt 20.)
k_tweak_templates, maggies[*,goodindx], maggies_ivar[*,goodindx], $
  gals[goodindx].redshift, avgcoeffs[*,goodindx], $
  avgspec, lambda, filterlist=['sdss_u0.par','sdss_g0.par','sdss_r0.par', $
                               'sdss_i0.par','sdss_z0.par'], $
  maggies_factor=maggies_factor, tweakpars=tweakpars, $
  vmatrix_tweaked=avgspec_tweaked

compare_maggies,avgcoeffs,gals.redshift,maggies,maggies_ivar, $
  lambda=lambda,vmatrix=avgspec_tweaked, $
  filterlist=['sdss_u0.par','sdss_g0.par','sdss_r0.par','sdss_i0.par', $
              'sdss_z0.par'], $
  recmaggies=recmaggies, filename='compare_maggies_avgspec_tweaked.ps'

endif
restore,'blah.sav'

stop 

ngauss=3
gauss_mixtures,ngauss,data,gauss_amp,gauss_mean,gauss_covar

fakedata=fake_catalog(-100L,gauss_amp,gauss_mean,gauss_covar,ndata=20000)
hogg_manyd_scatterplot,fltarr(20000)+1., $
  fakedata,'fake_eigencoeffs_scaled.ps',exponent=0.25,xnpix=30,ynpix=30, $
  range=range


; Now pick out some special cases:
indx0=where(eigencoeffs[2,*]/eigencoeffs[0,*] gt 0.12 and $
            eigencoeffs[2,*]/eigencoeffs[0,*] lt 0.18 and $
            eigencoeffs[1,*]/eigencoeffs[0,*] gt -0.03 and $
            eigencoeffs[1,*]/eigencoeffs[0,*] lt 0.03, count0)
indx1=where(eigencoeffs[2,*]/eigencoeffs[0,*] lt -0.12 and $
            eigencoeffs[2,*]/eigencoeffs[0,*] gt -0.20 and $
            eigencoeffs[1,*]/eigencoeffs[0,*] gt 0.230 and $
            eigencoeffs[1,*]/eigencoeffs[0,*] lt 0.270, count1)
indx2=where(eigencoeffs[1,*]/eigencoeffs[0,*] gt -0.60 and $
            eigencoeffs[1,*]/eigencoeffs[0,*] lt -0.50, count2)

;specs=vmatrix#coeffs
;avgspec=fltarr(n_elements(lambda)-1L,3)
;avgspec[*,0]=total(specs[*,indx0],2)/double(n_elements(indx0))
;avgspec[*,1]=total(specs[*,indx1],2)/double(n_elements(indx1))
;avgspec[*,2]=total(specs[*,indx2],2)/double(n_elements(indx2))
avgspec=fltarr(n_elements(lambda)-1L,3)
avgspec[*,0]=eigenmatrix#avgeigencoeffs[*,0]
avgspec[*,1]=eigenmatrix#avgeigencoeffs[*,1]
avgspec[*,2]=eigenmatrix#avgeigencoeffs[*,2]
avgcoeffs=k_fit_nonneg(maggies,maggies_ivar,avgspec, $
                       lambda,redshift=gals.redshift, $
                       filterlist=['sdss_u0.par','sdss_g0.par', $
                                   'sdss_r0.par','sdss_i0.par', $
                                   'sdss_z0.par'], $
                       chi2=chi2,rmatrix=rmatrix,zvals=zvals,maxiter=10000)

compare_maggies,avgcoeffs,gals.redshift,maggies,maggies_ivar, $
  lambda=lambda,vmatrix=avgspec, $
  filterlist=['sdss_u0.par','sdss_g0.par','sdss_r0.par','sdss_i0.par', $
              'sdss_z0.par'], $
  recmaggies=recmaggies, filename='compare_maggies_avgspec.ps'

stop

; Output the results 
compare_maggies,eigencoeffs,gals.redshift,maggies,maggies_ivar, $
  lambda=lambda,vmatrix=eigenmatrix, $
  filterlist=['sdss_u0.par','sdss_g0.par','sdss_r0.par','sdss_i0.par', $
              'sdss_z0.par'], $recmaggies=recmaggies, $
  filename='compare_maggies_eigen.ps'

hogg_manyd_scatterplot,fltarr(n_elements(gals))+1., $
  eigencoeffs,'eigencoeffs.ps',exponent=0.25,xnpix=40,ynpix=40

data=fltarr(npca+2,n_elements(gals))
for i=1,npca do $
  data[i+1,*]=eigencoeffs[i,*]/eigencoeffs[0,*]
data[1,*]=gals.redshift
ld=lumdis(gals.redshift,0.3,0.7)
data[0,*]=alog(eigencoeffs[0,*]*ld^2)
range=[[-9.5,-2.5], $
       [0.,0.5], $
       [-0.7,0.32], $
       [-0.24,0.24], $
       [-0.2,0.2], $
       [-0.09,0.09], $
       [-0.03,0.03]]
hogg_manyd_scatterplot,fltarr(n_elements(gals))+1., $
  data,'eigencoeffs_scaled.ps',exponent=0.25,xnpix=30,ynpix=30, $
  range=range
show_indx=where(data[2,*] gt 0.1)
hogg_manyd_scatterplot,fltarr(n_elements(show_indx))+1., $
  data[*,show_indx],'eigencoeffs_scaled_red.ps',exponent=0.25,xnpix=30, $
  ynpix=30,range=range

pts=fltarr(3,n_elements(gals))
for i=1,3 do $
  pts[i-1,*]=eigencoeffs[i,*]/eigencoeffs[0,*]

openw,unit,'train.pts',/get_lun
writeu,unit,pts
close,unit

tvlct,180.*(dindgen(256)+0.5)/256.,0.95+fltarr(256),0.95+fltarr(256),/hsv
tvlct,r,g,b,/get
r=double(r)/255.
g=double(g)/255.
b=double(b)/255.
scalevals=(dindgen(256)+0.5)/256.

scale=eigencoeffs[5,*]/eigencoeffs[0,*]
minscale=scale[(sort(scale))[long(n_elements(scale)*0.1)]]
maxscale=scale[(sort(scale))[long(n_elements(scale)*0.9)]]
scale=scale>(minscale)
scale=scale<(maxscale)
scale=1.-(scale-minscale)/(maxscale-minscale)

col=fltarr(4,n_elements(gals))
col[3,*]=1.
col[0,*]=interpol(r,scalevals,scale)
col[1,*]=interpol(g,scalevals,scale)
col[2,*]=interpol(b,scalevals,scale)

openw,unit,'train.col',/get_lun
writeu,unit,col
close,unit

stop

end
