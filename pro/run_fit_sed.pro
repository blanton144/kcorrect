;+
; NAME:
;   run_fit_sed
; PURPOSE:
;   Fit SED to broad-band colors and redshifts using k_fit_sed, given 
;   data in the VAGC.
; CALLING SEQUENCE:
;   run_fit_sed, outname
; INPUTS:
;   outname  -  name to attach to output files
; OPTIONAL INPUTS:
;   nophotozplates - don't use the photoz plates
;   maxiter - maximum number of iterations
;   outpath - directory for results
;   zlimits - limits for redshifts to use
;   nz - number of redshift bins to use for redshift histogram
;        equalization
;   templatelist - list of templates to begin with
;   nk - number of bands
;   nl - number of wavelength bins
;   lambdalim - wavelength range to consider
;   smoothtemplate - how much to smooth the initial templates
;   subsmoothtemplate - for a subrange, use a different smoothing
;   subsmoothlimits - limits of that subrange
;   nt - number of templates
;   shiftband - addition to mags to get AB mags
;   errband - "calibration" errors to associate with mags
;   errlimit - largest error to consider for fitting objects
;   maglimit - largest magnitude to consider for fitting objects
;   scale - fraction of objects to use of maximum possible
;   modelzlim - at what redshift to start using model mags
;   nozlim - don't consider objects in this range
;   etemplatepath - look here for the etemplates
; OUTPUTS:
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
;   k_fit_sed
;   k_write_ascii_table
; REVISION HISTORY:
;   05-Jan-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro run_fit_sed,outname,spfile=spfile,nophotozplates=nophotozplates,zlimits=zlimits,nz=nz,templatelist=templatelist,filtfile=filtfile,nl=nl,lambdalim=lambdalim,smoothtemplate=smoothtemplate,nt=nt,shiftband=shiftband,errband=errband,errlimit=errlimit,maglimit=maglimit,outpath=outpath, savfile=savfile, nk=nk,scale=scale, nsp=nsp,maxiter=maxiter, nozlim=nozlim,subsmoothlimits=subsmoothlimits,subsmoothtemplate=subsmoothtemplate,etemplatepath=etemplatepath,plotmaggies=plotmaggies,useconstraint=useconstraint, chi2=chi2, insavfile=insavfile, preset_ematrix=preset_ematrix

if(NOT keyword_set(nophotozplates)) then mustdo=[669,670,671,672] 

if(NOT keyword_set(maxiter)) then maxiter=30l
if(NOT keyword_set(outpath)) then outpath='.'
if(NOT keyword_set(zlimits)) then zlimits=[0.,0.5]
if(NOT keyword_set(nz)) then nz=8l
if(NOT keyword_set(templatelist)) then $
  templatelist=[ $
                 'ssp_salp_z02.flux.0220.dat', $
                 'ssp_salp_z02.flux.0104.dat', $
                 'ssp_salp_z02.flux.0123.dat', $
                 'ssp_salp_z02.flux.0154.dat', $
                 'ssp_salp_z02.flux.0093.dat', $
                 'ssp_salp_z008.flux.0220.dat', $
                 'ssp_salp_z008.flux.0104.dat', $
                 'ssp_salp_z008.flux.0123.dat', $
                 'ssp_salp_z008.flux.0154.dat', $
                 'ssp_salp_z008.flux.0093.dat' $
               ]
if(NOT keyword_set(nk)) then nk=5L
if(NOT keyword_set(nl)) then nl=500L
if(NOT keyword_set(lambdalim)) then lambdalim=[1000.,12000.]
if(NOT keyword_set(smoothtemplate)) then smoothtemplate=300.d
if(NOT keyword_set(subsmoothtemplate)) then subsmoothtemplate=150.d
if(NOT keyword_set(subsmoothlimits)) then subsmoothlimits=[3000.,5000.]
if(NOT keyword_set(nt)) then nt=4L
if(NOT keyword_set(shiftband)) then shiftband=[-0.042,0.036,0.015,0.013,-0.002]
if(NOT keyword_set(errband)) then errband=[0.05,0.02,0.02,0.02,0.03]
if(NOT keyword_set(errlimit)) then errlimit=dblarr(nk)+2.0d
if(NOT keyword_set(maglimit)) then maglimit=dblarr(nk)+24.0d
if(NOT keyword_set(scale)) then scale=1.d
if(NOT keyword_set(modelzlim)) then modelzlim=0.25
if(NOT keyword_set(nozlim)) then nozlim=[0.28,0.32]
if(NOT keyword_set(etemplatepath)) then $
  etemplatepath=getenv('KCORRECT_DIR')+'/data/etemplates' 

if(NOT keyword_set(filtfile)) then begin
    filtfile=etemplatepath+'/filterlist.'+outname+'.dat'
    spawn,'cat '+filtfile+' | wc -l',nfilters
    nk=nfilters[0]-1
    filterlist=strarr(nk)
    openr,unit,filtfile,/get_lun
    readf,unit,nk
    readf,unit,filterlist
    close,unit
    free_lun,unit
endif

; get data to search through
savfile='data_run_fit_sed.sav'
if(NOT file_test(savfile)) then begin

; get redshifts for all spectra
    sp=hogg_mrdfits(getenv('SDSS_VAGCDIR')+'/sdss_spectro_catalog.fits',1, $
                    columns=['sdss_spectro_z', 'sdss_spectro_zerr', $
                             'sdss_spectro_zwarning','sdss_spectro_class'],nrowchunk=40000)
    
; find imaging for all spectra objects
    im_matches=find_matches('sdss_spectro','sdss_imaging')
    
; get photometry
    im=hogg_mrdfits(getenv('SDSS_VAGCDIR')+'/sdss_imaging_catalog.fits',1, $
                    columns=['sdss_imaging_petrocounts', $
                             'sdss_imaging_petrocountserr', $
                             'sdss_imaging_counts_model', $
                             'sdss_imaging_counts_modelerr', $
                             'sdss_imaging_psfcounts', $
                             'sdss_imaging_psfcountserr', $
                             'sdss_imaging_fibercounts', $
                             'sdss_imaging_fibercountserr', $
                             'sdss_imaging_reddening'], nrowchunk=40000)
    im=im[im_matches]

    help,im,sp
    
    save,im,sp,filename=savfile
endif else begin
    restore,savfile
endelse

; choose objects
klog,'choose galaxies'
indx=where(sp.sdss_spectro_z gt zlimits[0] and $
           sp.sdss_spectro_z lt zlimits[1] and $
           strtrim(sp.sdss_spectro_class,2) eq 'GALAXY' and $
           im.sdss_imaging_petrocounts[2] gt 0.)
sp=sp[indx]
im=im[indx]
help,sp
help,im

; Cut down the sample in redshift
klog,'equalize redshift histogram'
num=lonarr(nz)
usesp=lonarr(n_elements(sp))
for i = 0l, nz-1l do begin
    zlo=zlimits[0]+double(i)*(zlimits[1]-zlimits[0])/double(nz)
    zhi=zlimits[0]+double(i+1l)*(zlimits[1]-zlimits[0])/double(nz)
    indx=where(sp.z gt zlo and sp.z lt zhi,count)
    num[i]=count
    klog,num[i]
endfor
nuse=long(double(num[nz-1l])*scale)
for i = 0l, nz-1l do begin
    zlo=zlimits[0]+double(i)*(zlimits[1]-zlimits[0])/double(nz)
    zhi=zlimits[0]+double(i+1l)*(zlimits[1]-zlimits[0])/double(nz)
    indx=where(sp.z gt zlo and sp.z lt zhi,count)
    if(i lt nz-1l) then begin
        indxuse=long(double(n_elements(indx))*randomu(seed,nuse))
        sortindxuse=indxuse[sort(indxuse)]
        uniqindxuse=sortindxuse[uniq(sortindxuse)]
    endif else begin
        uniqindxuse=lindgen(nuse)
    endelse
    usesp[indx[uniqindxuse]]=1
    klog,total(usesp)
endfor
if(keyword_set(mustdo)) then begin
    for i = 0l, n_elements(mustdo)-1l do begin
        mustindx=where(sp.sdss_spectro_plate eq mustdo[i],count)
        if(count gt 0) then usesp[mustindx]=1
    endfor
endif
indx=where(usesp gt 0)
sp=sp[indx]
im=im[indx]
help,sp
help,im

; trim out bad redshifts 
klog,'trim undesired redshifts'
indx=where(sp.sdss_spectro_z lt nozlim[0] or sp.sdss_spectro_z gt nozlim[1],count)
if(count gt 0) then begin
  sp=sp[indx]
  im=im[indx]
endif
help,sp
help,im

; use model where desired
klog,'use model '
mag=im.sdss_imaging_petrocounts
mag_err=im.sdss_imaging_petrocountserr
indx=where(sp.sdss_spectro_z ge modelzlim)
if(count gt 0) then begin
    klog,'using model'
    mag[*,indx]=im[indx].sdss_imaging_counts_model
    mag_err[*,indx]=im[indx].sdss_imaging_counts_modelerr
endif

; Trim off *anything* with bad errors, magnitudes
klog,'cutting bad errors'
goodindx=where(abs(mag_err[*,0]) lt errlimit[0] and $
               abs(mag_err[*,1]) lt errlimit[1] and $
               abs(mag_err[*,2]) lt errlimit[2] and $
               abs(mag_err[*,3]) lt errlimit[3] and $
               abs(mag_err[*,4]) lt errlimit[4] and $
               abs(mag[*,0]) lt maglimit[0] and $
               abs(mag[*,1]) lt maglimit[1] and $
               abs(mag[*,2]) lt maglimit[2] and $
               abs(mag[*,3]) lt maglimit[3] and $
               abs(mag[*,4]) lt maglimit[4])
sp=sp[goodindx]
im=im[goodindx]
mag=mag[*,goodindx]
mag_err=mag_err[*,goodindx]
help,sp
help,im
help,mag
help,mag_err

    
; Cut out weird colors (any 3-sigma points from the mean colors)
klog,'cutting weird colors'
for i = 0l, nk-2l do begin
    color=mag[i]-mag_err[i+1l]
    result=moment(color)
    klog,result[0],sqrt(result[1])
    goodindx=where((color-result[0])^2/result[1] lt 9.)
    sp=sp[goodindx]
    im=im[goodindx]
    mag=mag[*,goodindx]
    mag_err=mag_err[*,goodindx]
    help,sp
    help,im
    help,mag
    help,mag_err
endfor

; now make maggies to pass in 
galaxy_maggies=dblarr(nk,n_elements(sp))
galaxy_invvar=dblarr(nk,n_elements(sp))
for k=0l,nk-1l do begin
    galaxy_maggies[k,*]=10.d^(-0.4d*(mag[k]-im.sdss_imaging_reddening[k] $
                                     +shiftband[k]))
    galaxy_invvar[k,*]=galaxy_maggies[k,*]*0.4d*alog(10.d)* $
      sqrt(mag_err[k]^2+errband[k]^2)
    galaxy_invvar[k,*]=1.d/(galaxy_invvar[k,*]^2)
endfor

lambda=lambdalim[0]+dindgen(nl+1l)*(lambdalim[1]-lambdalim[0])/double(nl)
k_fit_sed,galaxy_maggies,galaxy_invvar,sp.sdss_spectro_z,templatelist, $
  filterlist, coeff, ematrix, bmatrix, bflux, lambda, nt=nt, $
  reconstruct_maggies=reconstruct_maggies, plotmaggies=plotmaggies, $
  smoothtemplate=smoothtemplate, subsmoothtemplate=subsmoothtemplate, $
  subsmoothlimits=subsmoothlimits, maxiter=maxiter, chi2=chi2, $
  useconstraint=useconstraint, preset_ematrix=preset_ematrix
z=sp.sdss_spectro_z

; get the ellipse of coeffs
ngals=n_elements(z)
scaled=coeff[1:nt-1L,*]/(replicate(1.,nt-1)#coeff[0,*])
meanscaled=total(scaled,2,/double)/double(ngals)
varscaled= 0d
for i=0L, n_elements(z)-1L do begin & delta=scaled[*,i]-meanscaled & varscaled= varscaled+delta#delta & endfor
varscaled=varscaled/double(ngals)

k_write_ascii_table,ematrix,outpath+'/ematrix.'+outname+'.dat'
k_write_ascii_table,bmatrix,outpath+'/bmatrix.'+outname+'.dat'
k_write_ascii_table,bflux,outpath+'/bflux.'+outname+'.dat'
k_write_ascii_table,lambda,outpath+'/lambda.'+outname+'.dat'
k_write_ascii_table,coeff,outpath+'/coeff.'+outname+'.dat'
k_write_ascii_table,z,outpath+'/z.'+outname+'.dat'
k_write_ascii_table,varscaled,outpath+'/scaledvar.'+outname+'.dat'
k_write_ascii_table,meanscaled,outpath+'/scaledmean.'+outname+'.dat'

end
;------------------------------------------------------------------------------
