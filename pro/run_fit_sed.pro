;+
; NAME:
;   run_fit_sed
;
; PURPOSE:
;   Fit SED to broad-band colors and redshifts using k_fit_sed.
;
; CALLING SEQUENCE:
;   k_fit_sed, outname
;
; INPUTS:
;   outname  -  name to attach to output files
;
; OPTIONAL INPUTS:
;   spfile  - spectro file (defaults to '/data/sdss/spectro/spAll.fits')
;   photozplates - include the photoz plates in the analysis
;   chunksize -  read in the spfile in chunks of this size (default 10000l)
;   zlimits - limits to apply to redshifts fit to
;   nz - number of redshift chunks for picking out objects evenly in z [8]
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
;   k_fit_sed
;   k_write_ascii_table
;
; REVISION HISTORY:
;   05-Jan-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro run_fit_sed,outname,spfile=spfile,nophotozplates=nophotozplates,chunksize=chunksize,zlimits=zlimits,nz=nz,templatelist=templatelist,filtfile=filtfile,nl=nl,lambdalim=lambdalim,smoothtemplate=smoothtemplate,nt=nt,fraction=fraction,shiftband=shiftband,errband=errband,errlimit=errlimit,maglimit=maglimit,outpath=outpath, savfile=savfile, nk=nk,scale=scale, nsp=nsp,maxiter=maxiter

if(NOT keyword_set(nophotozplates)) then mustdo=[669,670,671,672] 

if(NOT keyword_set(maxiter)) then maxiter=10l
if(NOT keyword_set(outpath)) then outpath='.'
if(NOT keyword_set(chunksize)) then chunksize=10000l
if(NOT keyword_set(zlimits)) then zlimits=[0.,0.5]
if(NOT keyword_set(nz)) then nz=8l
if(NOT keyword_set(templatelist)) then $
  templatelist=['ssp_salp_z02.flux.0220.dat', $
                'ssp_salp_z02.flux.0180.dat', $
                'ssp_salp_z02.flux.0150.dat', $
                'ssp_salp_z02.flux.0120.dat', $
                'ssp_salp_z004.flux.0220.dat', $
                'ssp_salp_z004.flux.0180.dat', $
                'ssp_salp_z004.flux.0150.dat', $
                'ssp_salp_z004.flux.0120.dat', $
                'flat.flux.dat', $
                'cos1.flux.dat']
if(NOT keyword_set(nk)) then nk=5L
if(NOT keyword_set(nl)) then nl=500L
if(NOT keyword_set(lambdalim)) then lambdalim=[1000.,12000.]
if(NOT keyword_set(smoothtemplate)) then smoothtemplate=300.d
if(NOT keyword_set(subsmoothtemplate)) then subsmoothtemplate=150.d
if(NOT keyword_set(subsmoothlimits)) then subsmoothlimits=[3000.,5000.]
if(NOT keyword_set(nt)) then nt=4L
if(NOT keyword_set(fraction)) then fraction=1.
if(NOT keyword_set(spfile)) then spfile='/data/sdss/spectro/spAll.fits'
if(NOT keyword_set(shiftband)) then shiftband=dblarr(nk)
if(NOT keyword_set(errband)) then errband=dblarr(nk)
if(NOT keyword_set(errlimit)) then errlimit=dblarr(nk)+0.8d
if(NOT keyword_set(maglimit)) then maglimit=dblarr(nk)+22.5d
if(NOT keyword_set(scale)) then scale=1.d
if(NOT keyword_set(nozlim)) then nozlim=[60.,61.]

if(NOT keyword_set(filtfile)) then begin
    filtfile=getenv('KCORRECT_DIR')+'/data/etemplates/filterlist.'+outname $
      +'.dat'
    spawn,'cat '+filtfile+' | wc -l',nfilters
    nk=nfilters[0]
    filterlist=strarr(nk)
    openr,unit,getenv('KCORRECT_DIR')+'/data/etemplates/filterlist.'+outname+ $
      '.dat',/get_lun
    readf,unit,filterlist
    close,unit
    free_lun,unit
endif

columns=['z','petrocounts','petrocountserr','reddening','class', $
         'ra', 'dec', 'plate']

lambda=lambdalim[0]+dindgen(nl+1l)*(lambdalim[1]-lambdalim[0])/double(nl)

; Read the necessary columns from spAll.fits
openr,unit,spfile,/get_lun
mrd_hread,unit,hdrstr
mrd_hread,unit,hdrstr
close,unit
free_lun,unit
hdr=hdr2struct(hdrstr)
if (NOT keyword_set(nsp)) then nsp=hdr.naxis2
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

; Cut down the sample in redshift
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

indx=where(sp.z lt nozlim[0] or sp.z gt nozlim[1],count)
if(count gt 0) then $
  sp=sp[indx]

; Trim off *anything* with bad errors, magnitudes
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
sp=sp[goodindx]
help,sp

; Cut out weird colors (any 3-sigma points from the mean colors)
for i = 0l, nk-2l do begin
    color=sp.petrocounts[i]-sp.petrocounts[i+1l]
    result=moment(color)
    klog,result[0],sqrt(result[1])
    goodindx=where((color-result[0])^2/result[1] lt 9.)
    sp=sp[goodindx]
endfor

galaxy_flux=dblarr(nk,n_elements(sp))
galaxy_invvar=dblarr(nk,n_elements(sp))
for k=0l,nk-1l do begin
    galaxy_flux[k,*]=10.d^(-0.4d*(sp.petrocounts[k]-sp.reddening[k]-17.d $
                                  +shiftband[k]))
    galaxy_invvar[k,*]=galaxy_flux[k,*]*0.4d*alog(10.d)* $
      sp.petrocountserr[k]
    galaxy_invvar[k,*]=1.d/(galaxy_invvar[k,*]^2+errband[k]^2)
endfor

k_fit_sed,galaxy_flux,galaxy_invvar,sp.z,templatelist, $
  filterlist, coeff, ematrix, bmatrix, bflux, lambda, nt=nt, $
  model_flux=model_flux, /plotfluxes, smoothtemplate=smoothtemplate, $
  maxiter=maxiter, subsmoothtemplate=subsmoothtemplate, $
  subsmoothlimits=subsmoothlimits
z=sp.z

k_write_ascii_table,ematrix,outpath+'/ematrix.'+outname+'.dat'
k_write_ascii_table,bmatrix,outpath+'/bmatrix.'+outname+'.dat'
k_write_ascii_table,bflux,outpath+'/bflux.'+outname+'.dat'
k_write_ascii_table,lambda,outpath+'/lambda.'+outname+'.dat'
k_write_ascii_table,coeff,outpath+'/coeff.'+outname+'.dat'
k_write_ascii_table,z,outpath+'/z.'+outname+'.dat'

savfile=outname+'.sav'
save,galaxy_flux,galaxy_invvar,z,coeff,ematrix,bmatrix,bflux,lambda,nt, $
  filename=savfile

outpts=outname+'.pts'
out=fltarr(nt-1l,n_elements(sp.z))
for i=0, nt-2l do $
  out[i,*]=coeff[i+1l,*]/coeff[0,*]
openw,11,outpath+'/'+outpts
writeu,11,out
close,11
out=0d
    
end
;------------------------------------------------------------------------------
