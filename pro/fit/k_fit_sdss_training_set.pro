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
;           check fit to SDSS ugriz of CNOC2 (right now CNOC2 doesn't
;           get g band right)
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   18-Jan-2003  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_fit_click, avgindx, nkmeans
avgindx[0]=avgindx[0]+1L
if(avgindx[0] ge nkmeans and n_elements(avgindx) gt 1) then begin
    tmp_avgindx=avgindx[1:n_elements(avgindx)-1L]
    k_fit_click, tmp_avgindx, nkmeans-1L
    avgindx[1:n_elements(avgindx)-1L]=tmp_avgindx
    avgindx[0]=avgindx[1]+1L
endif
end
;
pro k_fit_avgspec, avgindx, gals, group, lambda, eigenmatrix, eigencoeffs, $
                   navg, npca, avgcoeffs, avgspec, chi2gals,chi2, mean, $
                   median, filterlist=filterlist

;   create average spectra
    print,avgindx
    avgspec=fltarr(n_elements(lambda)-1L,navg)
    for i=0L, navg-1L do begin
        indx=where(group eq avgindx[i])
        avgeigencoeffs=fltarr(npca+1)
        avgeigencoeffs[0]=1.
        avgeigencoeffs[1:npca]= $
          total(eigencoeffs[1:npca,indx]/ $
                (replicate(1.,npca)#eigencoeffs[0,indx]),2)$
          /double(n_elements(indx))
        avgspec[*,i]=eigenmatrix#avgeigencoeffs
    endfor

;   qa plot
    set_plot,'x'
    !P.MULTI=[navg,1,navg]
    for i=0L, navg-1L do $
      plot,lambda,avgspec[*,i]

;   fit spectra to gals
    chi2gals=0
    rmatrix=0
    avgcoeffs=k_fit_nonneg(gals.maggies,gals.maggies_ivar, $
                           avgspec,lambda,redshift=gals.redshift, $
                           filterlist=filterlist, $
                           chi2=chi2gals,rmatrix=rmatrix, $
                           zvals=zvals, maxiter=10000)

;   record stats
    djs_iterstat,chi2gals,sigrej=3,mean=mean,median=median
    chi2=total(chi2gals,/double)
end
;
pro k_fit_sdss_training_set, name=name, navg=navg

if(NOT keyword_set(name)) then name='test'
if(NOT keyword_set(navg)) then navg=3
dontsave_name=name
dontsave_navg=navg

if(NOT keyword_set(npca)) then npca=5
if(NOT keyword_set(nkmeans)) then nkmeans=10
if(NOT keyword_set(sublmin)) then sublmin=2000.
if(NOT keyword_set(sublmax)) then sublmax=12000.
if(NOT keyword_set(filterlist)) then $
  filterlist=['sdss_u0.par','sdss_g0.par', $
              'sdss_r0.par','sdss_i0.par', $
              'sdss_z0.par','twomass_J.par', $
              'twomass_H.par','twomass_Ks.par']
  
savfile='data_k_fit_sdss_training_set.sav'

if(not file_test(savfile)) then begin
;   read in sdss stuff
    gals=(mrdfits('sdss_training_set.'+name+'.fits',1))

    k_read_ascii_table,vmatrix,'vmatrix.'+name+'.dat'
    k_read_ascii_table,lambda,'lambda.'+name+'.dat'
    nv=n_elements(vmatrix)/(n_elements(lambda)-1L)

;   HACK for reasonable numbers 
    vmatrix=vmatrix/1.e+38

;   fit nonnegative model
    use_indx=shuffle_indx(n_elements(gals),num_sub=2500)
    ;add_indx=where(gals.redshift gt 0.5 OR $
                   ;gals.sdss_spectro_tag eq -1)
    ;use_indx=[use_indx,add_indx]
    ;use_indx=lindgen(n_elements(gals))
    coeffs=k_fit_nonneg(gals[use_indx].maggies, $
                        gals[use_indx].maggies_ivar,vmatrix, $
                        lambda,redshift=gals[use_indx].redshift, $
                        filterlist=filterlist, $
                        chi2=chi2,rmatrix=rmatrix,zvals=zvals,maxiter=10000, $
                        /verbose)

    dontsave_navg=0
    dontsave_name=0
    save,filename=savfile
endif else begin
    splog,'restoring'
    restore,savfile
    splog,'done'
    if(keyword_set(dontsave_name)) then name=dontsave_name
    if(keyword_set(dontsave_navg)) then navg=dontsave_navg
endelse

orig_chi2=chi2

fixindx=where(gals[use_indx].maggies_ivar[0] gt 0. or $
              gals[use_indx].maggies_ivar[1] gt 0. or $
              gals[use_indx].maggies_ivar[2] gt 0. or $
              gals[use_indx].maggies_ivar[3] gt 0. or $
              gals[use_indx].maggies_ivar[4] gt 0. or $
              gals[use_indx].maggies_ivar[5] gt 0. or $
              gals[use_indx].maggies_ivar[6] gt 0. or $
              gals[use_indx].maggies_ivar[7] gt 0.)
use_indx=use_indx[fixindx]
coeffs=coeffs[*,fixindx]
orig_chi2=orig_chi2[fixindx]
gals=gals[use_indx]

;   for kicks, reconstruct maggies
compare_maggies,coeffs,gals.redshift,gals.maggies,gals.maggies_ivar, $
  rmatrix=rmatrix,zvals=zvals,recmaggies=recmaggies, $
  filename='compare_maggies.ps', psym=3, xra=[0.,0.8], $
  filterlist=filterlist

; find flux direction, and scale
k_ortho_templates,coeffs,vmatrix,lambda,bcoeffs,bmatrix,bflux,bdotv=bdotv, $
  bdotb=bdotb,sublmin=sublmin,sublmax=sublmax

; Scale distribution to plane of constant flux
flux_total=bcoeffs[0,*]
pcacoeffs=fltarr(nv-1,n_elements(gals))
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
eigencoeffs=fltarr(npca+1L,n_elements(gals))
eigencoeffs[0,*]=1.
eigencoeffs[1:npca,*]=transpose(tmp_eigenvec)#centcoeffs
for i=0,npca do $
  eigencoeffs[i,*]=eigencoeffs[i,*]*flux_total

; qa plots
loadct,0
compare_maggies,eigencoeffs,gals.redshift,gals.maggies, $
  gals.maggies_ivar, lambda=lambda,vmatrix=eigenmatrix, $
  filterlist=filterlist, xra=[0.,0.8], $ 
  recmaggies=recmaggies, filename='compare_maggies_eigen.ps'
data=fltarr(npca+1,n_elements(gals))
for i=0,npca-1 do $
  data[i,*]=eigencoeffs[i+1,*]/eigencoeffs[0,*]
data[npca,*]=gals.redshift
range=[[-0.7,0.32], $
       [-0.24,0.24], $
       [-0.2,0.2], $
       [-0.09,0.09], $
       [-0.03,0.03], $
       [0.,0.7]]
hogg_manyd_scatterplot,fltarr(n_elements(gals))+1., $
  data,'eigencoeffs_scaled.ps',exponent=0.25,xnpix=30,ynpix=30, $
  range=range

; k-means data
data=fltarr(npca,n_elements(gals))
for i=0,npca-1 do $
  data[i,*]=eigencoeffs[i+1,*]/eigencoeffs[0,*]
kmeans_streams, nkmeans, data, group, group_mean=group_mean

ndistinct=factorial(nkmeans)/(factorial(nkmeans-navg)*factorial(navg))
chi2=dblarr(ndistinct)+1.d+30
mean=dblarr(ndistinct)+1.d+30
median=dblarr(ndistinct)+1.d+30
avgindx=fltarr(navg,ndistinct)
tmp_avgindx=navg-1L-lindgen(navg)
for m=0L, ndistinct-1L do begin
    
    k_fit_avgspec, tmp_avgindx, gals, group, lambda, eigenmatrix, $
      eigencoeffs, navg, npca, avgcoeffs, avgspec, chi2gals,tmp_chi2, $
      tmp_mean, tmp_median, filterlist=filterlist

;   record results
    chi2[m]=tmp_chi2
    mean[m]=tmp_mean
    median[m]=tmp_median
    print,m,chi2[m],mean[m],median[m]
    
;   update avgindx
    avgindx[*,m]=tmp_avgindx
    k_fit_click, tmp_avgindx, nkmeans
endfor 

minmean=min(mean,minm)
tmp_avgindx=avgindx[*,minm]
k_fit_avgspec, tmp_avgindx, gals, group, lambda, eigenmatrix, eigencoeffs, $
  navg, npca, avgcoeffs, avgspec, chi2gals,tmp_chi2, tmp_mean, tmp_median, $
  filterlist=filterlist

compare_maggies,avgcoeffs,gals.redshift,gals.maggies,gals.maggies_ivar, $
  lambda=lambda,vmatrix=avgspec, filterlist=filterlist, xra=[0.,0.8], $
  recmaggies=recmaggies, filename='compare_maggies_avgspec_'+ $
  strtrim(string(navg),2)+'.ps'

goodindx=where(orig_chi2 lt 200. and gals.redshift lt 0.8)
k_tweak_templates, gals[goodindx].maggies, gals[goodindx].maggies_ivar, $
  gals[goodindx].redshift, avgcoeffs[*,goodindx], $
  avgspec, lambda, filterlist=filterlist, $
  maggies_factor=maggies_factor, tweakpars=tweakpars, $
  vmatrix_tweaked=avgspec_tweaked

compare_maggies,avgcoeffs,gals.redshift,gals.maggies,gals.maggies_ivar, $
  lambda=lambda,vmatrix=avgspec_tweaked, filterlist=filterlist, xra=[0.,0.8], $
  recmaggies=recmaggies, filename='compare_maggies_avgspec_tweaked_'+ $
  strtrim(string(navg),2)+'.ps'

set_print,filename='tweaked_specs_'+strtrim(string(navg),2)+'.ps'
!P.MULTI=[0,1,navg]
for i=0,navg-1L do begin & $
  plot,lambda,avgspec[*,i],/xlog,/ylog & $
  oplot,lambda,avgspec_tweaked[*,i],color=djs_icolor('red') & $
endfor
for i=0,navg-1L do $
  plot,lambda,avgspec_tweaked[*,i]/avgspec[*,i],/xlog
end_print

k_write_ascii_table, avgspec_tweaked, 'vmatrix.'+name+ $
  strtrim(string(navg),2)+'.dat'
k_write_ascii_table, lambda, 'lambda.'+name+strtrim(string(navg),2)+'.dat'

save,filename='done_k_fit_sdss_training_set_'+strtrim(string(navg),2)+'.sav'

end
