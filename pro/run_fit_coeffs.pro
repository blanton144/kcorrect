;+
; NAME:
;   run_fit_coeffs
;
; PURPOSE:
;   Fit SED to the spAll file using k_fit_coeffs, given eigentemplates as 
;   determined by k_fit_sed. 
;
; CALLING SEQUENCE:
;   k_fit_coeffs, outname
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL INPUT/OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   k_fit_coeffs
;
; REVISION HISTORY:
;   05-Jan-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro run_fit_coeffs,version,spfile=spfile,chunksize=chunksize,zlimits=zlimits,shiftband=shiftband,errband=errband,errlimit=errlimit,maglimit=maglimit,vpath=vpath,savfile=savfile,nsp=nsp, evenz=evenz, primtargetmask=primtargetmask,modelzlim=modelzlim, nprimtargetmask=nprimtargetmask

nk=5l
;if(NOT keyword_set(version)) then version='default'
;if(NOT keyword_set(vpath)) then vpath='.'
if(NOT keyword_set(savfile)) then savfile='.'+'/all.'+version+'.sav'
if(NOT keyword_set(zstep)) then zstep=0.08
if(NOT keyword_set(chunksize)) then chunksize=10000l
if(NOT keyword_set(spfile)) then spfile='/data/sdss/spectro/spAll.fits'
if(NOT keyword_set(shiftband)) then shiftband=dblarr(nk)
if(NOT keyword_set(errband)) then errband=[0.05,0.02,0.02,0.02,0.03]
if(NOT keyword_set(errlimit)) then errlimit=dblarr(nk)+0.8d
if(NOT keyword_set(maglimit)) then maglimit=dblarr(nk)+22.5d
if(NOT keyword_set(zlimits)) then zlimits=[0.,0.5]
if(NOT keyword_set(nz)) then nz=8l
if(NOT keyword_set(scale)) then scale=1.d
if(NOT keyword_set(modelzlim)) then modelzlim=0.25

columns=['z','petrocounts','petrocountserr','counts_model','counts_modelerr', $
         'reddening','class', 'ra', 'dec', 'plate','tile','fiberid', 'run', $
         'rerun', 'camcol','field','id','primtarget','fibercounts', $
         'fibercountserr']

; Read the necessary columns from spAll.fits
openr,unit,spfile,/get_lun
mrd_hread,unit,hdrstr
mrd_hread,unit,hdrstr
close,unit
free_lun,unit
hdr=hdr2struct(hdrstr)
if(NOT keyword_set(nsp)) then nsp=hdr.naxis2
nchunks=nsp/chunksize
for i = 0l, nchunks do begin
    nlo=i*chunksize
    nhi=(i+1l)*chunksize-1l
    if(nlo lt nsp) then begin
        if(nhi ge nsp) then nhi=nsp-1
        klog,nlo,nhi
        sptmp=mrdfits(spfile, 1, range=[nlo,nhi], columns=columns)
        indx=where(sptmp.class eq 'GALAXY' and $
                   sptmp.z gt zlimits[0] and $
                   sptmp.z lt zlimits[1] and $
                   (abs(sptmp.petrocounts[0,*]) lt 50. or $
                    abs(sptmp.petrocounts[1,*]) lt 50. or $
                    abs(sptmp.petrocounts[2,*]) lt 50. or $
                    abs(sptmp.petrocounts[3,*]) lt 50. or $
                    abs(sptmp.petrocounts[4,*]) lt 50.) and $
                   (abs(sptmp.counts_model[0,*]) lt 50. or $
                    abs(sptmp.counts_model[1,*]) lt 50. or $
                    abs(sptmp.counts_model[2,*]) lt 50. or $
                    abs(sptmp.counts_model[3,*]) lt 50. or $
                    abs(sptmp.counts_model[4,*]) lt 50.) and $
                   sptmp.petrocounts[2] gt 0.,count)
        if(count gt 0) then begin
            if(keyword_set(sp)) then begin
                sp=[sp,sptmp[indx]]
                sptmp=0l
            endif else begin 
                sp=sptmp[indx]
                help,/struct,sp
            endelse 
        endif
    endif
endfor

if(keyword_set(primtargetmask)) then begin
    indx=where((sp.primtarget and primtargetmask) gt 0,count)
    if(count eq 0) then return
    sp=sp[indx]
endif

if(keyword_set(nprimtargetmask)) then begin
    indx=where((sp.primtarget and nprimtargetmask) eq 0,count)
    if(count eq 0) then return
    sp=sp[indx]
endif

if(keyword_set(evenz)) then begin
; Cut down the sample in redshift
    seed=100
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
            mustindx=where(sp.plate eq mustdo[i],count)
            if(count gt 0) then usesp[mustindx]=1
        endfor
    endif
    indx=where(usesp gt 0)
    sp=sp[indx]
    help,sp
endif

; Find the good magnitudes, and determine the average colors in 
; bins of redshift; use these to reset bad determinations of the 
; flux to something more sensible
mags=sp.petrocounts
magserr=sp.petrocountserr
k_fix_mags,sp.z,mags,magserr,maglimit,errlimit,zstep
sp.petrocounts=mags
sp.petrocountserr=magserr
mags=sp.counts_model
magserr=sp.counts_modelerr
k_fix_mags,sp.z,mags,magserr,maglimit,errlimit,zstep
sp.counts_model=mags
sp.counts_modelerr=magserr

mags=sp.counts_model
magserr=sp.counts_modelerr
indx=where(sp.z lt modelzlim,count)
if(count gt 0) then begin
    mags[*,indx]=sp[indx].petrocounts
    magserr[*,indx]=sp[indx].petrocountserr
endif

galaxy_maggies=dblarr(nk,n_elements(sp))
galaxy_invvar=dblarr(nk,n_elements(sp))
for k=0l,nk-1l do begin
    galaxy_maggies[k,*]=10.d^(-0.4d*(mags[k,*]-sp.reddening[k] $
                                  +shiftband[k]))
    galaxy_invvar[k,*]=galaxy_maggies[k,*]*0.4d*alog(10.d)* $
      sqrt(magserr[k,*]^2+errband[k]^2)
    galaxy_invvar[k,*]=1.d/(galaxy_invvar[k,*]^2)
endfor

kcorrect,galaxy_maggies,galaxy_invvar,sp.z,recmags,coeff=coeff, $
  version=version, vpath=vpath, /maggies, /invvar, /addgrgap
galaxy_z=sp.z
galaxy_petrocounts=sp.petrocounts
galaxy_petrocountserr=sp.petrocountserr
galaxy_counts_model=sp.counts_model
galaxy_counts_modelerr=sp.counts_modelerr

save,galaxy_maggies,galaxy_invvar,sp,coeff,ematrix,bmatrix,lambda, $
  filterlist,filename=savfile

if(keyword_set(outpts)) then begin
    out=fltarr(nt-1l,n_elements(sp.z))
    for i=0, nt-2l do $
      out[i,*]=coeff[i+1l,*]/coeff[0,*]
    openw,11,outpath+'/'+outpts
    writeu,11,out
    close,11
    out=0d
endif

end
;------------------------------------------------------------------------------
