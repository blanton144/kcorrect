;+
; NAME:
;   lf_create_mock
; PURPOSE:
;   Given an evolving Schechter function luminosity function, creates
;   a catalog in the following way:
; 
;     0) Create the selection function for the most favorable
;        K-corrections
;     1) Sample a Gaussian random field on a large grid out past the
;        largest redshift, but following the selection function. 
;     2) Now choose luminosities and K-corrections for the galaxies. 
;     3) Apply uncertainties
;     4) Select the final sample
;
; INPUTS:
;    nsample - number of galaxies to draw
;    schechter - structure with {mstar, alpha, q, p}
;    mmin, mmax - flux limits
;    sample_absmmin, sample_absmmax - absolute mag limits
;    zvals_kcorrect, uniq_kcorrect - kcorrect parameters
; OPTIONAL INPUTS:
;    omega0 - matter density
;    omegal0 - vacuum energy density
; KEYWORDS:
; OUTPUTS:
; OPTIONAL OUTPUTS:
; BUGS:
; DEPENDENCIES:
; REVISION HISTORY:
;   2003-01-18  written - Blanton
;-
pro lf_create_mock,nsample,schechter,mmin,mmax,sample_absmmin,sample_absmmax, $
                   sample_zmin,sample_zmax, $
                   zvals_kcorrect,uniq_kcorrect,zzout=zzout, $
                   omega0=omega0,omegal0=omegal0

; settings
dh=2.99792e+5/100.
pi=!DPI

; defaults
if(n_elements(omega0) eq 0) then omega0=0.3
if(n_elements(omegal0) eq 0) then omegal0=0.7
if(n_elements(nzvals) eq 0) then nzvals=100


; create the selection function for ikcorrect=0 (assumed to be most
; favorable) 
zvals=sample_zmin+(sample_zmax-sample_zmin)*(dindgen(nzvals)+0.5)/ $
  double(nzvals)

; get the total expected in the survey volume
int_zvals=sample_zmin+(sample_zmax-sample_zmin)* $
  (dindgen(int_nzvals)+0.5)/double(int_nzvals)
dz=int_zvals[1]-int_zvals[0]
dvdz=dh^3*dcomvoldz(int_zvals,omega0,omegal0)*10.^(0.4*p*(int_zvals-zzero))
volume=total(dvdz*dz,/double)
totalphi=total(phi*0.5* $
               (errorf((model_absmmax-absmk)/(sqrt(2.)*sigabsmag))- $
                errorf((model_absmmin-absmk)/(sqrt(2.)*sigabsmag))),/double)
nexp=volume*area*totalphi
splog,nexp

; now get ntotal objects distributed in the volume
; according to the model fit
lf_sample_npgauss,ntotal,absmk,phi,p,q,sigabsmag,model_absmmin, $
  model_absmmax,sample_zmin,sample_zmax,absmout=absmout, $
  zzout=zzout,fpadjust=fpadjust,omega0=omega0, $
  omegal0=omegal0,zzero=zzero,rhozpars=rhozpars, $
  rhozfunc=rhozfunc,sample_nzvals=sample_nzvals

; assign random values of absmerr from the sample itself
ngals=n_elements(absmerr)
absmerrout=absmerr[long(randomu(seed,ntotal)*ngals)]

; adjust absolute magnitude according to errors
absmout=absmout+absmerrout*randomn(seed,ntotal)

; get a K-correction 
ikcorrectout=ikcorrect[long(randomu(seed,ntotal)*ngals)]
nkcorrect=n_elements(uniqkcorrect)/nzvals
kcorrectout=dblarr(ntotal)-100.
kcorrect_zvals=sample_zmin+(sample_zmax-sample_zmin)*(dindgen(nzvals)+0.5)/ $
  double(nzvals)
for kkcorrect=0L, nkcorrect-1L do begin
    indx=where(ikcorrectout eq kkcorrect,count)
    if(count gt 0) then begin
        iz=long((zzout[indx]-sample_zmin)/(sample_zmax-sample_zmin)* $
                double(nzvals))
        tlindx=where(iz lt 0,count)
        if(count gt 0) then iz[tlindx]=0
        thindx=where(iz gt nzvals-2,count)
        if(count gt 0) then iz[thindx]=nzvals-2
        izp1=iz+1
        sz=(zzout[indx]-kcorrect_zvals[iz])/(kcorrect_zvals[izp1]- $
                                             kcorrect_zvals[iz])
        kcorrectout[indx]=uniqkcorrect[iz,kkcorrect]+sz* $
          (uniqkcorrect[izp1,kkcorrect]-uniqkcorrect[iz,kkcorrect])
    endif
endfor

; now assign an observed apparent magnitude
dc=dh*comdis2(zzout,double(omega0),double(omegal0))
dl=(1.+zzout)*dc
dm=5.*alog10(dl)+25.
appmout=absmout+dm+kcorrectout

; check the flux limits
mminout=mmin
mmaxout=mmax
fgotout=fgot
indx=where(appmout gt mminout and appmout lt mmaxout,count)
if(count eq 0) then begin
    splog,'no galaxies survive flux cuts'
    absmout=-1
    zzout=-1
    appmout=-1
    absmerrout=-1
    ikcorrectout=-1
    kcorrectout=-1
    mminout=-1
    mmaxout=-1
    fgotout=-1
    return
endif
absmout=absmout[indx]
zzout=zzout[indx]
appmout=appmout[indx]
absmerrout=absmerrout[indx]
ikcorrectout=ikcorrectout[indx]
kcorrectout=kcorrectout[indx]
mminout=mminout[indx]
mmaxout=mmaxout[indx]
fgotout=fgotout[indx]

; check the observed absolute magnitude limits
indx=where(absmout gt sample_absmmin and absmout lt sample_absmmax,count)
if(count eq 0) then begin
    splog,'no galaxies survive observed absmag cuts'
    absmout=-1
    zzout=-1
    appmout=-1
    absmerrout=-1
    ikcorrectout=-1
    kcorrectout=-1
    mminout=-1
    mmaxout=-1
    fgotout=-1
    return
endif
absmout=absmout[indx]
zzout=zzout[indx]
appmout=appmout[indx]
absmerrout=absmerrout[indx]
ikcorrectout=ikcorrectout[indx]
kcorrectout=kcorrectout[indx]
mminout=mminout[indx]
mmaxout=mmaxout[indx]
fgotout=fgotout[indx]

; check fgot
indx=where(randomu(seed,n_elements(absmout)) lt fgotout,count)
if(count eq 0) then begin
    splog,'no galaxies survive fgot cuts'
    absmout=-1
    zzout=-1
    appmout=-1
    absmerrout=-1
    ikcorrectout=-1
    kcorrectout=-1
    mminout=-1
    mmaxout=-1
    fgotout=-1
    return
endif
absmout=absmout[indx]
zzout=zzout[indx]
appmout=appmout[indx]
absmerrout=absmerrout[indx]
ikcorrectout=ikcorrectout[indx]
kcorrectout=kcorrectout[indx]
mminout=mminout[indx]
mmaxout=mmaxout[indx]
fgotout=fgotout[indx]

; now set the weight of these (to recover actual number of objects)
fadjust=fpadjust*nexp/(double(ntotal))

end
