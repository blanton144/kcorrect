;+
; NAME:
;   k_make_sdss_training_set
; PURPOSE:
;   Extract training set from VAGC and create a save file
; CALLING SEQUENCE:
;   k_make_sdss_training_set [, modelzlim=, zlimits=, nzchunks=, $
;       shiftband=, errband=, errlimit=, maglimit=, /nophotozplates ]
; INPUTS:
; OPTIONAL INPUTS:
;   modelzlim - at greater than this redshift, use model mags (default 0.15)
;   zlimits - redshift limits to impose on sample (default 1.e-3, 0.5)
;   nzchunks - number of chunks used to equalize redshift dist (default 8)
;   shiftband - mag=mag+shiftband (to get to ABy system)
;               (default -0.042,0.036,0.015,0.013,-0.002)
;   errband - error to add in quadrature (default 0.05,0.02,0.02, 0.02, 0.03)
;   maglimit - magnitude limit to impose in each band
;              (default [30,30,30,30,30])
;   errlimit - error limit (default [1,1,1,1,1])
;   name - name of sample (default 'test')
; OUTPUTS:
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
; EXAMPLES:
; BUGS:
;   Needs to have VAGC extraction updated
; PROCEDURES CALLED:
; REVISION HISTORY:
;   18-Jan-2003  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_make_sdss_training_set,modelzlim=modelzlim,zlimits=zlimits, $
                             nzchunks=nzchunks, shiftband=shiftband, $
                             errband=errband, errlimit=errlimit, $
                             maglimit=maglimit, name=name

if(NOT keyword_set(nophotozplates)) then mustdoplates=[669,670,671,672] 
if(NOT keyword_set(zlimits)) then zlimits=[1.e-3,0.5]
if(NOT keyword_set(nzchunks)) then nzchunks=8
if(NOT keyword_set(shiftband)) then shiftband=[-0.042,0.036,0.015,0.013,-0.002]
if(NOT keyword_set(errband)) then errband=[0.05,0.02,0.02,0.02,0.03]
if(NOT keyword_set(errlimit)) then errlimit=fltarr(5)+1.0d
if(NOT keyword_set(maglimit)) then maglimit=fltarr(5)+30.0d
if(NOT keyword_set(modelzlim)) then modelzlim=0.15
if(NOT keyword_set(name)) then name='test'

; get data to search through
savfile='sdss_training_set.'+name+'.sav'
if(NOT file_test(savfile)) then begin
    
;   get redshifts for all spectra
    sp=hogg_mrdfits(getenv('SDSS_VAGCDIR')+'/sdss_spectro_catalog.fits',1, $
                    columns=['sdss_spectro_z', 'sdss_spectro_plate', $
                             'sdss_spectro_zerr','sdss_spectro_zwarning', $
                             'sdss_spectro_class','sdss_spectro_tag', $
                             'sdss_spectro_tag_primary'], nrowchunk=40000)
    
;   find imaging for all spectra objects
    im_matches=find_matches('sdss_spectro','sdss_imaging')

;   now select all the possible spectra we want
    sp_indx=where(sp.sdss_spectro_tag eq sp.sdss_spectro_tag_primary and $
                  sp.sdss_spectro_class eq 'GALAXY' and $
                  sp.sdss_spectro_z gt zlimits[0] and $
                  sp.sdss_spectro_z le zlimits[1] and $
                  sp.sdss_spectro_zwarning eq 0, sp_count) 
    
; get photometry
    im=hogg_mrdfits(getenv('SDSS_VAGCDIR')+'/sdss_imaging_catalog.fits',1, $
                    columns=['sdss_imaging_tag', $
                             'sdss_imaging_ra2000', $
                             'sdss_imaging_dec2000', $
                             'sdss_imaging_petrocounts', $
                             'sdss_imaging_petrocountserr', $
                             'sdss_imaging_counts_model', $
                             'sdss_imaging_counts_modelerr', $
                             'sdss_imaging_psfcounts', $
                             'sdss_imaging_psfcountserr', $
                             'sdss_imaging_fibercounts', $
                             'sdss_imaging_fibercountserr', $
                             'sdss_imaging_reddening'], nrowchunk=40000)

    help,im,sp
    matched_indx=where(im_matches[sp_indx] ge 0, matched_count)
    if(matched_count gt 0) then begin
        sp_indx=sp_indx[matched_indx]
        
        im=im[im_matches[sp_indx]]
        sp=sp[sp_indx]
        
        help,im,sp
        save,im,sp,filename=savfile 
    endif else begin
        klog,'no matches to imaging.'
        return
    endelse
endif else begin
    restore,savfile
endelse

; now equalize the redshift bins (but include all mustdo plates)
mustdo=lonarr(n_elements(im))
for i=0, n_elements(mustdoplates)-1L do begin
    inplate_indx=where(sp.sdss_spectro_plate eq mustdoplates[i],inplate_count)
    if(inplate_count gt 0) then mustdo[inplate_indx]=1
endfor
includegal=lonarr(n_elements(im))
ninchunk=fltarr(nzchunks)
for i=nzchunks-1L,0,-1 do begin
    zchunk_lo=zlimits[0]+(zlimits[1]-zlimits[0])*float(i)/float(nzchunks)
    zchunk_hi=zlimits[0]+(zlimits[1]-zlimits[0])*float(i+1)/float(nzchunks)
    chunk_indx=where(sp.sdss_spectro_z gt zchunk_lo and $
                     sp.sdss_spectro_z lt zchunk_hi, chunk_count)
    ninchunk[i]=chunk_count
    if(i ne nzchunks-1L) then begin
        dice=randomu(seed,ninchunk[i])
        subsample=float(ninchunk[nzchunks-1L])/float(ninchunk[i])
        help,i,subsample[0]
        include_indx=where(dice lt subsample[0], include_count)
        if(include_count gt 0) then includegal[chunk_indx[include_indx]]=1
    endif else begin
        includegal[chunk_indx]=1L
    endelse
endfor
include_indx=where(includegal gt 0 or mustdo gt 0,include_count)
if(include_count eq 0) then begin
    klog,'no remaining galaxies after equalizing redshifts.'
    return
endif
im=im[include_indx]
sp=sp[include_indx]
help,im

; set up which magnitudes are used
klog,'use model past z='+string(modelzlim)
mag=im.sdss_imaging_petrocounts-im.sdss_imaging_reddening
mag_err=im.sdss_imaging_petrocountserr
model_indx=where(sp.sdss_spectro_z ge modelzlim,model_count)
help,model_count
if(model_count gt 0) then begin
    mag[*,model_indx]=im[model_indx].sdss_imaging_counts_model- $
      im[model_indx].sdss_imaging_reddening
    mag_err[*,model_indx]=im[model_indx].sdss_imaging_counts_modelerr
endif

; cut out bad errors
klog,'cutting bad errors'
goodindx=where(abs(mag_err[0,*]) lt errlimit[0] and $
               abs(mag_err[1,*]) lt errlimit[1] and $
               abs(mag_err[2,*]) lt errlimit[2] and $
               abs(mag_err[3,*]) lt errlimit[3] and $
               abs(mag_err[4,*]) lt errlimit[4] and $
               abs(mag[0,*]) lt maglimit[0] and $
               abs(mag[1,*]) lt maglimit[1] and $
               abs(mag[2,*]) lt maglimit[2] and $
               abs(mag[3,*]) lt maglimit[3] and $
               abs(mag[4,*]) lt maglimit[4])
sp=sp[goodindx]
im=im[goodindx]
mag=mag[*,goodindx]
mag_err=mag_err[*,goodindx]
help,im

; finally, shift the bands, add errors in quadrature
for i=0L, n_elements(shiftband)-1L do $
  mag[i,*]=mag[i,*]+shiftband[i]
for i=0L, n_elements(errband)-1L do $
  mag_err[i,*]=sqrt(mag_err[i,*]^2+errband[i]^2)

; now output 
outfile='sdss_training_set.'+name+'.fits'
hdr=strarr(1)
sxaddpar,hdr,'DATE',systime(),'date of creation'
sxaddpar,hdr,'KVERSION',k_version(),'version of kcorrect'
outstr1={mag:fltarr(5), mag_ivar:fltarr(5), redshift:0.D, ra:0.D, dec:0.D, $ 
         objc_rowc:0.D, objc_colc:0.D, sdss_imaging_tag:long64(0L), $
         sdss_spectro_tag:long64(0L)}
outstr=replicate(outstr1,n_elements(sp))
outstr.ra=im.sdss_imaging_ra2000
outstr.dec=im.sdss_imaging_dec2000
outstr.redshift=sp.sdss_spectro_z
outstr.mag=mag
outstr.mag_ivar=1./mag_err^2
outstr.sdss_imaging_tag=im.sdss_imaging_tag
outstr.sdss_spectro_tag=sp.sdss_spectro_tag
mwrfits,dummy,outfile,hdr,/create
mwrfits,outstr,outfile

end
;------------------------------------------------------------------------------
