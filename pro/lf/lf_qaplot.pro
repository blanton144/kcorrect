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
pro lf_poisson_sample

; get names
cd,current=cwd
path=strsplit(cwd,'/',/extract)
sample=path[n_elements(path)-3]
letter=path[n_elements(path)-2]
post=path[n_elements(path)-1]
if(NOT keyword_set(samplename)) then samplename=sample+letter+post

; use a previous fit if it exists
restore,'npgauss_fit.'+samplename+'.sav'
dh=2.99792e+5/100.

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

sample_vmin=comdis2(sample_zmin,omega0,omegal0)^3
sample_vmax=comdis2(sample_zmax,omega0,omegal0)^3
zvals=sample_zmin+(sample_zmax-sample_zmin)* $
  (dindgen(nzvals)+0.5)/double(nzvals)
vvals=comdis2(zvals,omega0,omegal0)^3

; choose redshift 
nsub=100
grandomi_full=long(nphi*randomu(seed,ngals*nsub))
grandomc_full=randomu(seed,ngals*nsub)
vvrandom_full=sample_vmin+randomu(seed,ngals*nsub)*(sample_vmax-sample_vmin)
zzrandom_full=dblarr(ngals*nsub)
for i=0L, ngals*nsub-1L do begin
    indx=where(vvals lt vvrandom_full[i])
    j=indx[n_elements(indx)-1]
    if(j eq -1) then j=0
    if(j eq nzvals-1) then j=nzvals-2
    jp1=j+1
    sv=(vvrandom_full[i]-vvals[j])/(vvals[jp1]-vvals[j])
    zzrandom_full[i]=zvals[j]+sv*(zvals[jp1]-zvals[j])
endfor

gchances=outphi/max(outphi)
for i=0L, ngals-1L do begin
; choose absolute magnitude 
    usegals=i*nsub+lindgen(nsub)
    grandomc=grandomc_full[usegals]
    grandomi=grandomi_full[usegals]
    zzrandom=zzrandom_full[usegals]
    indx=where(grandomc lt gchances[grandomi],count)
    if(count gt 0) then begin
        zzrandom=zzrandom[indx]
        amrandom=absmk[grandomi[indx]]+sqrt(sigabsmag^2+useabsmerr[i]^2)* $
          randomn(seed,count)-q*(zzrandom-zzero)

; are we within the limits? get dm, K, then check
        dcrandom=dh*comdis2(zzrandom,omega0,omegal0,cdtable=cdtable, $
                            zmaxtable=zmaxtable)
        dlrandom=(1.+zzrandom)*dcrandom
        dmrandom=5.*alog10(dlrandom)+25.
        izzrandom=long(double(nzvals)*(zzrandom-sample_zmin)/ $
                       (sample_zmax-sample_zmin))
        krandom=uniqkcorrect[izzrandom,4]
        mrandom=amrandom+dmrandom+krandom
        indx=where(mrandom gt mmin[i] and mrandom lt mmax[i] and $
                   amrandom gt sample_absmmin and amrandom lt sample_absmmax, $
                   count)
        if(count gt 0) then begin
            if(n_elements(amlist) eq 0) then begin
                amlist=amrandom[indx]
                zzlist=zzrandom[indx]
            endif else begin
                amlist=[amlist,amrandom[indx]]
                zzlist=[zzlist,zzrandom[indx]]
            endelse
        endif
    endif
endfor
stop

end
