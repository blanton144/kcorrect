;+
; NAME:
;   lf_selfunc_npgauss
; PURPOSE:
;   calculate selection function at a given redshift using
;   an npgauss model
; INPUTS:
;   zz               [N] redshifts
;   mmin
;   mmax
;   sample_absmmin   minimum abs mag for the sample
;   sample_absmmax   maximum abs mag for the sample
;   phi
;   absmk
;   sigabsmag        
;   q
;   p
; OPTIONAL INPUTS:
;   zzero            zeropoint in redshift (default 0.1)
; KEYWORDS:
; OUTPUTS:
; OPTIONAL OUTPUTS:
; BUGS:
; DEPENDENCIES:
; REVISION HISTORY:
;   2002-5-22  written - Blanton
;-
;
function lf_selfunc_npgauss,zz,kcorrect,mmin,mmax,sample_absmmin, $
                            sample_absmmax,model_absmmin,model_absmmax, $
                            phi,absmk,sigabsmag,q,p, $
                            zzero=zzero,omega0=omega0,omegal0=omegal0, $
                            rhozpars=rhozpars,rhozfunc=rhozfunc, $
                            zvals=zvals,sample_zmin=sample_zmin, $
                            sample_zmax=sample_zmax,absmmin=absmmin, $
                            absmmax=absmmax

; defaults
if(n_elements(zzero) eq 0) then zzero=0.1
if(n_elements(omega0) eq 0) then omega0=0.3
if(n_elements(omegal0) eq 0) then omegal0=0.7

; settings
pi=3.14159265358979D
dh=2.99792D+5/100.D
ngals=n_elements(zz)
nphi=n_elements(phi)

; get distance moduli
dc=dh*comdis2(zz,omega0,omegal0,cdtable=cdtable,zmaxtable=zmaxtable)
dl=(1.+zz)*dc
dm=5.*alog10(dl)+25.

; get min and max *model* absolute mags
; note that given a flux limit mmax, at higher redshift 
absmmin=mmin-dm-kcorrect+q*(zz-zzero)
absmmax=mmax-dm-kcorrect+q*(zz-zzero)

; make sure it is not out of sample bounds
absmmin=(absmmin > (sample_absmmin+q*(zz-zzero)))
absmmin=(absmmin < (sample_absmmax+q*(zz-zzero)))
absmmax=(absmmax > (sample_absmmin+q*(zz-zzero)))
absmmax=(absmmax < (sample_absmmax+q*(zz-zzero)))

; make sure it is not out of model bounds
absmmin=(absmmin > model_absmmin)
absmmin=(absmmin < model_absmmax)
absmmax=(absmmax > model_absmmin)
absmmax=(absmmax < model_absmmax)

; calculate selection function for each object
numer=dblarr(ngals)
denom=dblarr(ngals)
for i=0L, nphi-1L do begin
    denom_mmindiff=model_absmmin-absmk[i]
    denom_mmaxdiff=model_absmmax-absmk[i]
    denom=denom+phi[i]*0.5* $
      (errorf(denom_mmaxdiff/(sqrt(2.)*sigabsmag))- $
       errorf(denom_mmindiff/(sqrt(2.)*sigabsmag)))
    numer_mmindiff=absmmin-absmk[i]
    numer_mmaxdiff=absmmax-absmk[i]
    numer=numer+phi[i]*0.5* $
      (errorf(numer_mmaxdiff/(sqrt(2.)*sigabsmag))- $
       errorf(numer_mmindiff/(sqrt(2.)*sigabsmag)))
endfor
selfunc=numer/denom
selfunc=(selfunc > 0.)

; adjust for number density evolution
selfunc=selfunc*10.^(0.4*p*(zz-zzero))

; if desired, take out the large scale structure
if(n_elements(rhozpars) gt 0) then begin
    rhoz_int=1.+rhozfunc#rhozpars 
    indx=where(rhoz_int le 0.,count)
    if(count gt 0) then rhoz_int[indx]=1.e-3
    nzvals=n_elements(zvals)
    rhoz=interpolate(rhoz_int,(double(nzvals))*(zz-sample_zmin)/ $
                     (sample_zmax-sample_zmin)-0.5)
    selfunc=selfunc*rhoz
endif

return,selfunc

end
