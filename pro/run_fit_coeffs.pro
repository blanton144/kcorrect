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
pro run_fit_coeffs,version,spfile=spfile,chunksize=chunksize,zlimits=zlimits,shiftband=shiftband,errband=errband,errlimit=errlimit,maglimit=maglimit,vpath=vpath,savfile=savfile,nsp=nsp

nk=5l
if(NOT keyword_set(version)) then version='default'
if(NOT keyword_set(vpath)) then vpath='.'
if(NOT keyword_set(savfile)) then savfile=vpath+'/all.'+version+'.sav'
if(NOT keyword_set(zstep)) then zstep=0.08
if(NOT keyword_set(chunksize)) then chunksize=10000l
if(NOT keyword_set(spfile)) then spfile='/data/sdss/spectro/spAll.fits'
if(NOT keyword_set(shiftband)) then shiftband=dblarr(nk)
if(NOT keyword_set(errband)) then errband=[0.05,0.02,0.02,0.02,0.03]
if(NOT keyword_set(errlimit)) then errlimit=dblarr(nk)+0.8d
if(NOT keyword_set(maglimit)) then maglimit=dblarr(nk)+22.5d
if(NOT keyword_set(zlimits)) then zlimits=[0.,0.5]

columns=['z','petrocounts','petrocountserr','counts_model','counts_modelerr', $
         'reddening','class', 'ra', 'dec', 'plate','tile','fiberid', 'run', $
         'rerun', 'camcol','field','id','primtarget']

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
mags=0.d
magserr=0.d

galaxy_maggies=dblarr(nk,n_elements(sp))
galaxy_invvar=dblarr(nk,n_elements(sp))
for k=0l,nk-1l do begin
    galaxy_maggies[k,*]=10.d^(-0.4d*(sp.petrocounts[k]-sp.reddening[k] $
                                  +shiftband[k]))
    galaxy_invvar[k,*]=galaxy_maggies[k,*]*0.4d*alog(10.d)* $
      sp.petrocountserr[k]
    galaxy_invvar[k,*]=1.d/(galaxy_invvar[k,*]^2+errband[k]^2)
endfor

k_fit_coeffs,galaxy_maggies,galaxy_invvar,sp.z,coeff,version=version, $
  vpath=vpath,ematrix=ematrix,bmatrix=bmatrix,lambda=lambda, $
  filterlist=filterlist
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
