;+
; NAME:
;   lf_sample_npgauss
; PURPOSE:
;   create a sampling in redshift and absolute magnitude based on npgauss lf
; USAGE:
;   lf_sample_npgauss, ndesired, absmk, phi, p, q, sigabsmag, model_absmmin, $
;      model_absmmax, sample_zmin, sample_zmax [ , absmout=, zzout=, $
;      fpadjust=, omega0=, omegal0=, zzero=, nzvals=, rhozpars=, $
;      rhozfunc=, zvals=, sample_nzvals=]
; INPUTS:
;   ndesired - number of samples wanted 
;   absmk, phi, p, q, sigabsmag - description of npgauss lf
;   model_absmmin, model_absmmax - abs. mag. limits of model
;   sample_zmin, sample_zmax - redshift limits for sampling
; OPTIONAL INPUTS:
;   nzvals - number of redshift in kcorrect templates (default 50)
;   zvals - [nzvals] redshift grid
;   zzero - reference redshift for q and p evolution (default 0.1)
;   omega0 - matter density (default 0.3)
;   omegal0 - vacuum energy density (default 0.7)
;   rhozpars, rhozfunc - if we have fit out the radial density function
; OPTIONAL OUTPUTS:
;   absmout - [nout] abs. mag. of final galaxies
;   zzout - [nout] redshifts of final galaxies
;   fadjust - ratio of number in actual universion to from number
;             gotten in mock to number in actual (to compare redshift
;             histograms, say)
; REVISION HISTORY:
;   2002-5-22  written - Blanton
;-
pro lf_sample_npgauss,ndesired,absmk,phi,p,q,sigabsmag,model_absmmin, $
                      model_absmmax,sample_zmin,sample_zmax,absmout=absmout, $
                      zzout=zzout,fpadjust=fpadjust,omega0=omega0, $
                      omegal0=omegal0,zzero=zzero,nzvals=nzvals, $
                      rhozpars=rhozpars,rhozfunc=rhozfunc,zvals=zvals, $
                      sample_nzvals=sample_nzvals

if(NOT keyword_set(nzvals)) then nzvals=10000
if(n_elements(zzero) eq 0) then zzero=0.1
if(n_elements(omega0) eq 0) then omega0=0.3
if(n_elements(omegal0) eq 0) then omegal0=0.7
if(n_elements(sample_nzvals) eq 0) then sample_nzvals=10000
dh=2.99792e+5/100.
pi=3.14159265358979D+0

nphi=n_elements(absmk)

; set up gaussian choosing stuff
indx=where(absmk-4.*sigabsmag lt model_absmmax and $
           absmk+4.*sigabsmag gt model_absmmin)
gchances=dblarr(n_elements(phi))
gchances[indx]=phi[indx]*0.5* $
  (errorf((model_absmmax-absmk[indx])/(sqrt(2.)*sigabsmag))- $
   errorf((model_absmmin-absmk[indx])/(sqrt(2.)*sigabsmag)))
gchances=gchances/max(gchances)

; set up redshift choosing stuff
sample_dcmin=(dh*comdis(sample_zmin,double(omega0),double(omegal0)))
sample_dcmax=(dh*comdis(sample_zmax,double(omega0),double(omegal0)))
sample_vmin=sample_dcmin^3.
sample_vmax=sample_dcmax^3.
zvals=sample_zmin+(sample_zmax-sample_zmin)* $
  (dindgen(sample_nzvals+1))/double(sample_nzvals)
dcvals=(dh*comdis(zvals,double(omega0),double(omegal0)))
dcmap=sample_dcmin+(sample_dcmax-sample_dcmin)* $
  (dindgen(sample_nzvals+1))/double(sample_nzvals)
zmap=dblarr(sample_nzvals+1)
for i=0L, sample_nzvals do begin
    indx=where(dcvals lt dcmap[i])
    j=indx[n_elements(indx)-1]
    if(j le -1) then j=0
    if(j ge sample_nzvals-1) then j=sample_nzvals-2
    jp1=j+1
    sdc=(dcmap[i]-dcvals[j])/(dcvals[jp1]-dcvals[j])
    zmap[i]=zvals[j]+sdc*(zvals[jp1]-zvals[j])
endfor
zchances=dblarr(sample_nzvals+1)
for i=0, sample_nzvals do $
  zchances[i]=10.^(0.4*p*(zmap[i]-zzero))
if(n_elements(rhozpars) gt 0) then begin
    rhoz_int=1.+rhozfunc#rhozpars
    rhoz_zvals=sample_zmin+(sample_zmax-sample_zmin)* $
      (dindgen(n_elements(rhoz_int))+0.5)/double(n_elements(rhoz_int))
    rhoz_chances=interpol(rhoz_int,rhoz_zvals,zmap)
    zchances=zchances*rhoz_chances
endif
fpadjust=max(zchances)          ; multiply any distribution by this number to 
                                ; get actual distribution
zchances=zchances/fpadjust
fpadjust=1.                     ; in fact, not necessary since we return 
                                ; ntotal points, distributed according
                                ; to 10^(0.4*P*(z-zzero))

npicked=0
nleft=ndesired
absmout=dblarr(ndesired)
zzout=dblarr(ndesired)
while(nleft gt 0) do begin

; choose a gaussian for each object
    npick=nleft
    grandomc=randomu(seed,npick)
    grandomi=long(nphi*randomu(seed,npick))
    indx=where(grandomc lt gchances[grandomi],count)
    
    if(count gt 0) then begin
        npick=n_elements(indx)
        grandomi=grandomi[indx]
; assign a redshift
        vvrandom=sample_vmin+randomu(seed,npick)*(sample_vmax-sample_vmin)
        dcrandom=vvrandom^(1./3.)
        jvals=long(double(sample_nzvals)*(dcrandom-sample_dcmin)/ $
                   (sample_dcmax-sample_dcmin))
        jp1vals=jvals+1
        sdcvals=(dcrandom-dcmap[jvals])/(dcmap[jp1vals]-dcmap[jvals])
        zzrandom=zmap[jvals]+sdcvals*(zmap[jp1vals]-zmap[jvals])
        szvals=(zzrandom-zmap[jvals])/(zmap[jp1vals]-zmap[jvals])
        zchances_here=zchances[jvals]+ $
          szvals*(zchances[jp1vals]-zchances[jvals])
        zrandomc=randomu(seed,npick)

        indx=where(zrandomc lt zchances_here,count)
        if(count gt 0) then begin
            npick=n_elements(indx)
            grandomi=grandomi[indx]
            zzrandom=zzrandom[indx]
            
            amrandom=absmk[grandomi]+sigabsmag*randomn(seed,npick)
            outindx=where(amrandom lt model_absmmin or $
                          amrandom gt model_absmmax, $
                          outcount)
            while(outcount gt 0) do begin
                amrandom[outindx]= $
                  absmk[grandomi[outindx]]+sigabsmag*randomn(seed,outcount)
                outindx=where(amrandom lt model_absmmin or $
                              amrandom gt model_absmmax, $
                              outcount)
            end
            
            amrandom=amrandom-q*(zzrandom-zzero)
            absmout[npicked:npicked+npick-1L]=amrandom
            zzout[npicked:npicked+npick-1L]=zzrandom
            npicked=npicked+npick
            nleft=nleft-npick
        endif
    endif

    help,npicked
    help,nleft
end

end
