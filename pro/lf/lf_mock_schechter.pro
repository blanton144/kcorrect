;+
; NAME:
;   lf_mock_schechter
; PURPOSE:
;   create mock catalog given schechter luminosity function
; USAGE:
;   lf_mock_npgauss, ntotal, mmin, mmax, fgot, area, schechter, p, q, $
;      sigabsmag, sample_absmmin, sample_absmmax, model_absmmin, $
;      model_absmmax, sample_zmin,sample_zmax, ikcorrect, $
;      uniqkcorrect, absmerr [, absmout=, zzout=, appmout=, absmerrout=, $
;      ikcorrectout=, kcorrectout=, fadjust=, omega0=, omegal0=, $
;      nzvals=, zzero=, rhozpars=, rhozfunc=,zvals=,sample_nzvals= ]
; INPUTS:
;   ntotal - total number of galaxies in volume (not total output!)
;   mmin, mmax - [ntotal] flux limits of each galaxy
;   fgot - [ntotal] chances of observing each galaxy
;   area - total area of sample
;   schechter, p, q - description of schechter l.f. + evolution
;   sample_absmmin, sample_absmmax - abs. mag. limits of sample
;   model_absmmin, model_absmmax - abs. mag. limits to model
;   sample_zmin, sample_zmax - redshift limits of sample
;   ikcorrect - [ngals] available kcorrections to choose from
;   uniqkcorrect - [nzvals, nk] available k-correction functions
;   absmerr - [ngals] available uncertainties to choose from
; OPTIONAL INPUTS:
;   nzvals - number of redshift in kcorrect templates (default 50)
;   zvals - [nzvals] redshift grid
;   zzero - reference redshift for q and p evolution (default 0.1)
;   omega0 - matter density (default 0.3)
;   omegal0 - vacuum energy density (default 0.7)
;   rhozpars, rhozfunc - if we have fit out the radial density function
; OPTIONAL OUPUTS:
;   absmout - [nout] abs. mag. of final galaxies
;   zzout - [nout] redshifts of final galaxies
;   appmout - [nout] apparent mag of final galaxies
;   absmerrout - [nout] uncertainties on abs. mag. of final galaxies
;   ikcorrectout - [nout] which kcorrection functions for final galaxies
;   kcorrectout - [nout] actual kcorrections for final galaxies
;   fadjust - ratio of number in actual universion to from number
;             gotten in mock to number in actual (to compare redshift
;             histograms, say)
; REVISION HISTORY:
;   2002-5-22  written - Blanton
;-
pro lf_mock_schechter,ntotal,mmin,mmax,fgot,area,schechter,p,q, $
                      sample_absmmin,sample_absmmax,sample_zmin,sample_zmax, $
                      ikcorrect,uniqkcorrect,absmerr, absmout=absmout, $
                      zzout=zzout,appmout=appmout,absmerrout=absmerrout, $
                      ikcorrectout=ikcorrectout,kcorrectout=kcorrectout, $
                      fadjust=fadjust,omega0=omega0,omegal0=omegal0, $
                      nzvals=nzvals,zzero=zzero,rhozpars=rhozpars, $
                      rhozfunc=rhozfunc,zvals=zvals,sample_nzvals=sample_nzvals

; defaults
if(n_elements(omega0) eq 0) then omega0=0.3
if(n_elements(omegal0) eq 0) then omegal0=0.7
if(n_elements(nzvals) eq 0) then nzvals=50
if(n_elements(int_nzvals) eq 0) then int_nzvals=1500
if(n_elements(zzero) eq 0) then zzero=0.1
dh=2.99792e+5/100.
pi=3.14159265358979D+0

; get the total expected in the survey volume
int_zvals=sample_zmin+(sample_zmax-sample_zmin)* $
  (dindgen(int_nzvals)+0.5)/double(int_nzvals)
dz=int_zvals[1]-int_zvals[0]
dvdz=dh^3*dcomvoldz(int_zvals,omega0,omegal0)*10.^(0.4*p*(int_zvals-zzero))
volume=total(dvdz*dz,/double)
totalphi=lf_schechter_nden(schechter,limits=[sample_absmmin,sample_absmmax])
nexp=volume*area*totalphi
splog,nexp

; now get ntotal objects distributed in the volume
; according to the model fit
lf_sample_schechter,ntotal,schechter,p,q,sample_absmmin, $
  sample_absmmax,sample_zmin,sample_zmax,absmout=absmout, $
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
