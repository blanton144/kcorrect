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
    sp=get_catalog('sdss/sdss_spectro', $
                   columns=['z', 'plate', 'zerr','zwarning', $
                            'class','sdss_spectro_tag', $
                            'sdss_spectro_tag_primary'], nrowchunk=40000)
    
;   find imaging for all spectra objects
    im_matches=find_matches('sdss_spectro','sdss_imaging')
    twomass_matches=find_matches('sdss_spectro','twomass')

;   now select all the possible spectra we want
    sp_indx=where(sp.sdss_spectro_tag eq sp.sdss_spectro_tag_primary and $
                  sp.class eq 'GALAXY' and $
                  sp.z gt zlimits[0] and $
                  sp.z le zlimits[1] and $
                  sp.zwarning eq 0, sp_count) 
    
; get photometry
    im=get_catalog('sdss/sdss_imaging', $
                   columns=['sdss_imaging_tag', $
                            'ra', $
                            'dec', $
                            'petrocounts', $
                            'petrocountserr', $
                            'counts_model', $
                            'counts_modelerr', $
                            'psfcounts', $
                            'psfcountserr', $
                            'fibercounts', $
                            'fibercountserr', $
                            'reddening'], nrowchunk=40000)

; get twomass photometry
    twomass=get_catalog('twomass/twomass', $
                        columns=['twomass_tag', $
                                 'k_m_k20fe', $
                                 'k_msig_k20fe', $
                                 'h_m_k20fe', $
                                 'h_msig_k20fe', $
                                 'j_m_k20fe', $
                                 'j_msig_k20fe'],nrowchunk=40000)
    
    help,im,sp,twomass
    matched_indx=where(im_matches[sp_indx] ge 0, matched_count)
    if(matched_count gt 0) then begin
        sp_indx=sp_indx[matched_indx]

        twomass=twomass[twomass_matches[sp_indx]]
        im=im[im_matches[sp_indx]]
        sp=sp[sp_indx]
        
        help,im,sp,twomass
        save,im,sp,twomass,filename=savfile 
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
    inplate_indx=where(sp.plate eq mustdoplates[i],inplate_count)
    if(inplate_count gt 0) then mustdo[inplate_indx]=1
endfor
includegal=lonarr(n_elements(im))
ninchunk=fltarr(nzchunks)
for i=nzchunks-1L,0,-1 do begin
    zchunk_lo=zlimits[0]+(zlimits[1]-zlimits[0])*float(i)/float(nzchunks)
    zchunk_hi=zlimits[0]+(zlimits[1]-zlimits[0])*float(i+1)/float(nzchunks)
    chunk_indx=where(sp.z gt zchunk_lo and $
                     sp.z lt zchunk_hi, chunk_count)
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
twomass=twomass[include_indx]
sp=sp[include_indx]
help,im

; set up which magnitudes are used
klog,'use model past z='+string(modelzlim)
mag=im.petrocounts-im.reddening
mag_err=im.petrocountserr
model_indx=where(sp.z ge modelzlim,model_count)
help,model_count
if(model_count gt 0) then begin
    mag[*,model_indx]=im[model_indx].counts_model- $
      im[model_indx].reddening
    mag_err[*,model_indx]=im[model_indx].counts_modelerr
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
twomass=twomass[goodindx]
mag=mag[*,goodindx]
mag_err=mag_err[*,goodindx]
help,im

; finally, shift the bands, add errors in quadrature
for i=0L, n_elements(shiftband)-1L do $
  mag[i,*]=mag[i,*]+shiftband[i]
for i=0L, n_elements(errband)-1L do $
  mag_err[i,*]=sqrt(mag_err[i,*]^2+errband[i]^2)

; add CNOC2 data (NEED TO APPLY REDDENING AND AB SHIFTS -- and to
;                                                          2MASS)
cnoc2=mrdfits(getenv('KCORRECT_DIR')+'/data/redshifts/cnoc2/cnoc2.fits',1);
cnoc2_obj=mrdfits(getenv('KCORRECT_DIR')+ $
                  '/data/redshifts/cnoc2/cnoc2-obj.fits',1)
cnoc2_childobj=mrdfits(getenv('KCORRECT_DIR')+ $
                       '/data/redshifts/cnoc2/cnoc2-childobj.fits',1)
cnoc2_indx=where(cnoc2_obj.matchdist*3600. lt 10.)
cnoc2=cnoc2[cnoc2_indx]
cnoc2_obj=cnoc2_obj[cnoc2_indx]
cnoc2_childobj=cnoc2_childobj[cnoc2_indx]

; now output 
; make into maggies here ....
outfile='sdss_training_set.'+name+'.fits'
hdr=strarr(1)
sxaddpar,hdr,'DATE',systime(),'date of creation'
sxaddpar,hdr,'KVERSION',k_version(),'version of kcorrect'
outstr1={maggies:fltarr(8), maggies_ivar:fltarr(8), redshift:0.D, ra:0.D, $
         dec:0.D, objc_rowc:0.D, objc_colc:0.D, sdss_imaging_tag:long64(0L), $
         sdss_spectro_tag:long64(0L)}
outstr=replicate(outstr1,n_elements(sp)+n_elements(cnoc2))
isdss=0L+lindgen(n_elements(sp))
icnoc2=n_elements(sp)+lindgen(n_elements(cnoc2))
outstr[isdss].ra=im.ra
outstr[isdss].dec=im.dec
outstr[isdss].redshift=sp.z
outstr[isdss].maggies[0:4]=10.^(-0.4*mag[0:4])
outstr[isdss].maggies_ivar[0:4]=1./ $
  (mag_err[0:4]*outstr[isdss].maggies[0:4]*0.4*alog(10.))^2
twomass_indx=where(twomass.j_m_k20fe gt 0. and $
                   twomass.h_m_k20fe gt 0. and $
                   twomass.k_m_k20fe gt 0., twomass_count)
if(twomass_count gt 0) then begin
    outstr[isdss[twomass_indx]].maggies[5]= $
      10.^(-0.4*(twomass[twomass_indx].j_m_k20fe+ $
                 (k_vega2ab(filterlist='twomass_J.par',/kurucz))[0]))
    outstr[isdss[twomass_indx]].maggies_ivar[5]= $
      1./(0.4*alog(10.)*outstr[isdss[twomass_indx]].maggies[5]* $
          twomass[twomass_indx].j_msig_k20fe)^2
    outstr[isdss[twomass_indx]].maggies[6]= $
      10.^(-0.4*(twomass[twomass_indx].h_m_k20fe+ $
                 (k_vega2ab(filterlist='twomass_H.par',/kurucz))[0]))
    outstr[isdss[twomass_indx]].maggies_ivar[6]= $
      1./(0.4*alog(10.)*outstr[isdss[twomass_indx]].maggies[6]* $
          twomass[twomass_indx].h_msig_k20fe)^2
    outstr[isdss[twomass_indx]].maggies[7]= $
      10.^(-0.4*(twomass[twomass_indx].k_m_k20fe+ $
                 (k_vega2ab(filterlist='twomass_Ks.par',/kurucz))[0]))
    outstr[isdss[twomass_indx]].maggies_ivar[7]= $
      1./(0.4*alog(10.)*outstr[isdss[twomass_indx]].maggies[7]* $
          twomass[twomass_indx].k_msig_k20fe)^2
endif
outstr[isdss].sdss_imaging_tag=im.sdss_imaging_tag
outstr[isdss].sdss_spectro_tag=sp.sdss_spectro_tag
outstr[icnoc2].ra=cnoc2_childobj.ra
outstr[icnoc2].dec=cnoc2_childobj.dec
outstr[icnoc2].redshift=cnoc2.z
outstr[icnoc2].maggies[0:4]=cnoc2_childobj.modelflux
outstr[icnoc2].maggies_ivar[0:4]=cnoc2_childobj.modelflux_ivar
outstr[icnoc2].sdss_imaging_tag=-1
outstr[icnoc2].sdss_spectro_tag=-1
mwrfits,dummy,outfile,hdr,/create
mwrfits,outstr,outfile

end
;------------------------------------------------------------------------------
