;+
; NAME:
;   lf_poisson_sample
; PURPOSE:
;   plot the expected distribution of absolute magnitudes 
;   and the actual distribution for a redshift slice
; INPUTS:
; OPTIONAL INPUTS:
; KEYWORDS:
; OUTPUTS:
; OPTIONAL OUTPUTS:
; BUGS:
; DEPENDENCIES:
; REVISION HISTORY:
;   2002-5-22  written - Blanton
;-
pro lf_poisson_sample,absmk,phi,p,q,sigabsmag,sample_absmmin,sample_absmmax, $
                      sample_zmin,sample_zmax,zvals,nzvals,ikcorrect, $
                      uniqkcorrect,mmin,mmax,zzlist,amlist,zzero=zzero, $
                      omega0=omega0,omegal0=omegal0,nsub=nsub

if(NOT keyword_set(nsub)) then nsub=50
if(n_elements(omega0) eq 0) then omega0=0.3
if(n_elements(omegal0) eq 0) then omegal0=0.7
dh=2.99792e+5/100.
pi=3.14159265358979D

ngals=n_elements(mmin)
nphi=n_elements(phi)

nborder=4
outphi=dblarr(nphi+nborder)
outphi[0:nphi-1]=(phi[0L:nphi-1L])
outphi[nphi+lindgen(nborder)]=phi[nphi-1]+0.*(dindgen(nborder)+1.)* $
  (phi[nphi-1]-phi[nphi-2])

absmk=dblarr(nphi+nborder)
absmk[0:nphi-1]=sample_absmmin+(sample_absmmax-sample_absmmin)* $
  (dindgen(nphi)+0.5)/double(nphi)
absmk[nphi+lindgen(nborder)]=absmk[nphi-1]+ $
  (dindgen(nborder)+1.)*(absmk[1]-absmk[0])

usenzvals=nzvals*20

sample_vmin=(dh*comdis2(sample_zmin,omega0,omegal0))^3
sample_vmax=(dh*comdis2(sample_zmax,omega0,omegal0))^3
zvals=sample_zmin+(sample_zmax-sample_zmin)* $
  (dindgen(usenzvals)+0.5)/double(usenzvals)
vvals=(dh*comdis2(zvals,omega0,omegal0))^3

vmap=sample_vmin+(sample_vmax-sample_vmin)* $
  (dindgen(usenzvals+1))/double(usenzvals)
zmap=dblarr(usenzvals+1)
for i=0L, usenzvals do begin
    indx=where(vvals lt vmap[i])
    j=indx[n_elements(indx)-1]
    if(j eq -1) then j=0
    if(j eq usenzvals-1) then j=usenzvals-2
    jp1=j+1
    sv=(vmap[i]-vvals[j])/(vvals[jp1]-vvals[j])
    zmap[i]=zvals[j]+sv*(zvals[jp1]-zvals[j])
endfor
zchances=dblarr(usenzvals+1)
for i=0, usenzvals do $
  zchances[i]=10.^(0.4*p*(zmap[i]-zzero))
zchances=zchances/max(zchances)

; choose redshift 

splog,'here'
gchances=outphi/max(outphi)
ntotal=0L
nin=0L
for i=0L, ngals-1L do begin
    if((i mod 1000) eq 0) then splog,i
; choose absolute magnitude 
    grandomc=randomu(seed,nsub)
    grandomi=long((nphi+nborder)*randomu(seed,nsub))
    zrandomc=randomu(seed,nsub)
    vvrandom=sample_vmin+randomu(seed,nsub)*(sample_vmax-sample_vmin)
    jvals=long(double(usenzvals)*(vvrandom-sample_vmin)/ $
               (sample_vmax-sample_vmin))
    jp1vals=jvals+1
    svvals=(vvrandom-vmap[jvals])/(vmap[jp1vals]-vmap[jvals])
    zzrandom=zmap[jvals]+svvals*(zmap[jp1vals]-zmap[jvals])
    szvals=(zzrandom-zmap[jvals])/(zmap[jp1vals]-zmap[jvals])
    zchances_here=zchances[jvals]+szvals*(zchances[jp1vals]-zchances[jvals])

    indx=where(grandomc lt gchances[grandomi] and $
               zrandomc lt zchances_here,count)
    if(count gt 0) then begin
        zzrandom=zzrandom[indx]
        amrandom=absmk[grandomi[indx]]+sigabsmag*randomn(seed,count)
        indx=where(amrandom gt sample_absmmin and amrandom lt sample_absmmax, $
                   count)
        if(count gt 0) then begin
            amrandom=amrandom[indx]
            zzrandom=zzrandom[indx]

; now we are left with things distributed with the right
; probability in the M-z plane, but have not yet chopped 
; it to the flux limits; we will have to keep the number 
; of these objects for the normalization
            ntotal=ntotal+count
            amrandom=amrandom+absmerr[i]*randomn(seed,count)- $
              q*(zzrandom-double(zzero))
            
; apply errors

; are we within the limits? get dm, K, then check
            dcrandom=dh*comdis2(zzrandom,omega0,omegal0,cdtable=cdtable, $
                                zmaxtable=zmaxtable)
            dlrandom=(1.+zzrandom)*dcrandom
            dmrandom=5.*alog10(dlrandom)+25.
            izzrandom=long(double(nzvals)*(zzrandom-sample_zmin)/ $
                           (sample_zmax-sample_zmin))
            krandom=kcorrect[izzrandom,ikcorrect[i]]
            mrandom=amrandom+dmrandom+krandom
            indx=where(mrandom gt mmin[i] and mrandom lt mmax[i],count)
            if(count gt 0) then begin
                if(nin eq 0L) then begin
                    amlist=amrandom[indx]
                    zzlist=zzrandom[indx]
                endif else begin
                    amlist=[amlist,amrandom[indx]]
                    zzlist=[zzlist,zzrandom[indx]]
                endelse
                nin=nin+count
            endif
        endif
    endif
endfor
splog,ngals,nin,ntotal,total(phi),total(outphi)

; calculate normalization times square degrees; ie. multiply this by
; the number of square degrees in the sample to get the correct
; normalization factor
vtot=4.*pi*(sample_vmax-sample_vmin)/3.
fullspheredeg2=4.*pi*(180./pi)^2
norm=fullspheredeg2*double(ntotal)*double(ngals)/double(n_elements(amlist)) $
  /vtot

end
