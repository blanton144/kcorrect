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
pro run_fit_coeffs,outname,spfile=spfile,chunksize=chunksize,filterlist=filterlist,shiftband=shiftband,errband=errband,errlimit=errlimit,maglimit=maglimit,outpath=outpath,savfile=savfile

if(NOT keyword_set(outpath)) then outpath='.'
if(NOT keyword_set(zstep)) then zstep=0.08
if(NOT keyword_set(chunksize)) then chunksize=10000l
if(NOT keyword_set(spfile)) then spfile='/data/sdss/spectro/spAll.fits'
if(NOT keyword_set(shiftband) then shiftband=dblarr(nk)
if(NOT keyword_set(errband) then errband=dblarr(nk)
if(NOT keyword_set(errlimit) then errlimit=dblarr(nk)+0.8d
if(NOT keyword_set(maglimit) then errlimit=dblarr(nk)+30.5d

columns=['z','petrocounts','petrocountserr','reddening','class', $
         'ra', 'dec', 'plate']

k_load_ascii_table,ematrix,outpath+'/ematrix.'+outname+'.dat'
k_load_ascii_table,bmatrix,outpath+'/bmatrix.'+outname+'.dat'
k_load_ascii_table,lambda,outpath+'/lambda.'+outname+'.dat'

filtfile=getenv('KCORRECT_DIR')+'/data/etemplates/filterlist.'+outname+'.dat'
spawn,'cat '+filtfile+' | wc -l',nfilters
nk=nfilters[0]
filterlist=strarr(nk)
openr,unit,getenv('KCORRECT_DIR')+'/data/etemplates/filterlist.'+outname+ $
  '.dat',/get_lun
readf,unit,filterlist
close,unit
free_lun,unit
filterlist=getenv('KCORRECT_DIR')+'/data/filters/'+filterlist

; Read the necessary columns from spAll.fits
openr,unit,spfile,/get_lun
mrd_hread,unit,hdrstr
mrd_hread,unit,hdrstr
close,unit
free_lun,unit
hdr=hdr2struct(hdrstr)
nsp=hdr.naxis2
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
;sp.petrocounts[1]=sp.petrocounts[1]+0.08

; Find the good magnitudes, and determine the average colors in 
; bins of redshift; use these to reset bad determinations of the 
; flux to something more sensible
goodindx=where(abs(sp.petrocountserr[0]) lt errlimit[0] and $
               abs(sp.petrocountserr[1]) lt errlimit[1] and $
               abs(sp.petrocountserr[2]) lt errlimit[2] and $
               abs(sp.petrocountserr[3]) lt errlimit[3] and $
               abs(sp.petrocountserr[4]) lt errlimit[4] and $
               abs(sp.petrocounts[0]) lt maglimit[0] and $
               abs(sp.petrocounts[1]) lt maglimit[1] and $
               abs(sp.petrocounts[2]) lt maglimit[2] and $
               abs(sp.petrocounts[3]) lt maglimit[3] and $
               abs(sp.petrocounts[4]) lt maglimit[4])
good=lonarr(n_elements(sp))
good[goodindx]=1l
nz=long((max(sp.z)-min(sp.z))/zstep)
zbounds=min(sp.z)+(max(sp.z)-min(sp.z))*dindgen(nz+1l)/double(nz)
avgcolors=dblarr(n_elements(zbounds),nk-1l)
for i=0l, nz-1l do begin
    indxz=where(sp[goodindx].z ge zbounds[i] and $
                sp[goodindx].z le zbounds[i+1],count)
    if(count gt 0) then begin 
        for k=0l, nk-2l do $
          avgcolors[i,k]= $
          djs_avsigclip(sp[goodindx[indxz]].petrocounts[k] $
                        -sp[goodindx[indxz]].petrocounts[k+1])
    endif else begin
        for k=0l, nl-2l do $
          avgcolors[i,k]=0.5
    endelse
endfor
badindx=where(good ne 1l,count)
if(count gt 0) then begin
    for i=0, count-1 do begin
        zindx=long(double(nz)*(sp[i].z-min(sp.z))/(max(sp.z)-min(sp.z)))
        tmpindx=where(sp[badindx[i]].petrocountserr-errlimit lt 0.d and $
                      sp[badindx[i]].petrocounts-maglimit lt 0.d,ngood)
        isbad=lonarr(nk)+1l
        isbad[tmpindx]=0l
        while (ngood lt nk) do begin
; go backwards through the bands (for ugriz case, bases mags on most
; stable bands
            for k=ngood-1l, 0l, -1l do begin
                if(tmpindx[k] gt 0) then begin
                    if(isbad[tmpindx[k]-1]) then begin
                        sp[badindx[i]].petrocounts[tmpindx[k]-1l]= $
                          sp[badindx[i]].petrocounts[tmpindx[k]]+ $
                          avgcolor[zindx,tmpindx[k]-1]
                        sp[badindx[i]].petrocountserr[tmpindx[k]-1l]= $
                          errlimit[tmpindx[k]-1l]
                    endif
                endif
                if(tmpindx[k] lt nk-1) then begin
                    if(isbad[tmpindx[k]+1]) then begin
                        sp[badindx[i]].petrocounts[tmpindx[k]+1l]= $
                          sp[badindx[i]].petrocounts[tmpindx[k]]- $
                          avgcolor[zindx,tmpindx[k]]
                        sp[badindx[i]].petrocountserr[tmpindx[k]-1l]= $
                          errlimit[tmpindx[k]+1l]
                    endif
                endif
            endfor
            tmpindx=where(sp[badindx[i]].petrocountserr-errlimit lt 0.d and $
                          sp[badindx[i]].petrocounts-maglimit lt 0.d,ngood)
        end
    endfor
endif
sp=sp[goodindx]
help,sp

galaxy_flux=dblarr(nk,n_elements(sp))
galaxy_invvar=dblarr(nk,n_elements(sp))
for k=0l,nk-1l do begin
    galaxy_flux[k,*]=10.d^(-0.4d*(sp.petrocounts[k]-sp.reddening[k]-17.d $
                                  +shiftband[k]))
    galaxy_invvar[k,*]=galaxy_flux[k,*]*0.4d*alog(10.d)* $
      sp.petrocountserr[k]
    galaxy_invvar[k,*]=1.d/(galaxy_invvar[k,*]^2+errband[k]^2)
endfor

k_fit_coeffs,galaxy_flux,galaxy_invvar,sp.z,coeff,ematrix=ematrix, $
  zvals=zvals, filterlist=filterlist, bmatrix=bmatrix, lambda=lambda
z=sp.z

k_write_ascii_table,coeff,outpath+'/coeff.'+outname+'.dat'
k_write_ascii_table,z,outpath+'/z.'+outname+'.dat'

if(keyword_set(savfile)) then $
  save,galaxy_flux,galaxy_invvar,z,coeff,ematrix,bmatrix,bflux,lambda,nt, $
  filename=savfile

if(keyword_set(outpts) then begin
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
