;+
; NAME:
;   k_gmr_ktable
; PURPOSE:
;   Put together estimated K-corrections in ugriz bands for a 
;   set of g-r colors. 
; CALLING SEQUENCE:
;   k_gmr_ktable,filename
; INPUTS:
;   filename - file to create
; OPTIONAL INPUTS:
; OUTPUTS:
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   18-Jan-2003  Blanton (NYU)
;-
;------------------------------------------------------------------------------
pro k_gmr_ktable,ntypes=ntypes, version=version, vpath=vpath

; defaults
if(n_elements(ntypes) eq 0) then ntypes=10
if(n_elements(a3lo) eq 0) then a3lo=-0.35
if(n_elements(a3hi) eq 0) then a3hi=0.09
if(n_elements(nz) eq 0) then nz=101
if(n_elements(zlo) eq 0) then zlo=0.
if(n_elements(zhi) eq 0) then zhi=0.5

; initialize stuff
bands=['u','g','r','i','z']
nk=n_elements(bands)
mag=dblarr(nk)
magerr=dblarr(nk)+0.03
zz=0.
kcorrect,mag,magerr,zz,kcorrect,version=version,vpath=vpath, $
  rmatrix=rmatrix,ematrix=ematrix,zvals=zvals,bmatrix=bmatrix, $
  lambda=lambda,coeff=coeff
nt=n_elements(coeff)

; choose a range of coefficients
coeff=dblarr(nt,ntypes)
coeff[0,*]=1.
coeff[3,*]=a3lo+(a3hi-a3lo)*dindgen(ntypes)/double(ntypes-1)

; now reconstruct the magnitudes of the galaxy at a range of 
; redshifts
zz=zlo+(zhi-zlo)*dindgen(nz)/(double(nz-1)) 
mags=dblarr(nk,ntypes,nz)
for i=0L, nz-1L do begin
    k_reconstruct_maggies,coeff,replicate(zz[i],ntypes), maggies, $
      ematrix=ematrix,zvals=zvals,rmatrix=rmatrix,bmatrix=bmatrix, $
      lambda=lambda
    mags[*,*,i]=-2.5*alog10(maggies)
endfor
mags0=mags[*,*,0]
for i=0L, nz-1L do mags[*,*,i]=mags[*,*,i]-mags0
for i=0L, ntypes-1L do mags0[*,i]=mags0[*,i]-mags0[2,i]

; output magnitudes
magstr1={colors, bandminusr:dblarr(5)}
magstr=replicate(magstr1,ntypes)
magstr.bandminusr=mags0
pdata=ptr_new(magstr)
hdr=strarr(3)
hdr[0]='# Colors wrt r-band at z=0 of templates used for k_gmr_ktable.dat'
hdr[1]='# Created '+systime()
hdr[2]='kversion '+k_version()
yanny_write,'k_gmr_ktable.colors.dat',pdata,hdr=hdr

kstr1={kcorrect, redshift:0., kcorrect:dblarr(ntypes)}
kstr=replicate(kstr1,nz)
for i=0L, nk-1L do begin
    kstr.redshift=zz
    kstr.kcorrect=reform(mags[i,*,*],ntypes,nz)
    pdata=ptr_new(kstr)
    hdr=strarr(3)
    hdr[0]='# K-corrections for '+bands[i]+'-band for set of templates'
    hdr[1]='# Created '+systime()
    hdr[2]='kversion '+k_version()
    yanny_write,'k_gmr_ktable.'+bands[i]+'.dat',pdata,hdr=hdr
endfor

end
;------------------------------------------------------------------------------
