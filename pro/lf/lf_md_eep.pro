;+
; NAME:
;   lf_md_eep
; PURPOSE:
;   calculate normalization of EEP fit from lf_eep
; USAGE:
;   lf_md_eep, marea, zz, sel, zvals, selfunc, fraclimits, n1, n1err [, $
;        j3=, /simple, weight=, omega0=, omegal0=]
; INPUTS:
;   marea - [Nm] completeness-weight area of each region with a
;           constant flux limit
;   zz - [Ngals] redshift of each galaxy
;   sel - [Ngals] selection function for each galaxy
;   zvals - [Nz] redshift grid 
;   selfunc - [Nz, Nm] selection function grid in each region 
;   fraclimits - fraction of area in each region
; OPTIONAL INPUTS:
;   j3 - J_3 value of Davis&Huchra (1983) for optimal weighting scheme
;        (default to 10000.)
;   weight - weight to apply to each galaxy (default to 1.)
;   omega0 - matter density (default 0.3)
;   omegal0 - vacuum energy density (default 0.7)
; OPTIONAL KEYWORDS:
;   /simple - don't do optimal weighting, do simple weighting
; OUTPUTS:
;   n1 - normalization of EEP fit
;   n1err - uncertainty in normalization
; REVISION HISTORY:
;   2002-5-22  written - Blanton
;-
;
pro lf_md_eep,marea,zz,sel,zvals,selfunc,fraclimits,n1,n1err,j3=j3, $
              simple=simple, weight=weight, omega0=omega0, omegal0=omegal0

; defaults
if(n_elements(weight) eq 0) then weight=fltarr(n_elements(zz))+1.
if(n_elements(omega0) eq 0) then omega0=0.3
if(n_elements(omegal0) eq 0) then omegal0=0.7
if(n_elements(j3) eq 0) then j3=10.^4 ; ie. make it pretty irrelevant

; settings
pi=3.14159265358979D
dh=2.99792D+5/100.D
corrfact=1.

; get volume of survey
nzvals=n_elements(zvals)
nm=n_elements(selfunc)/nzvals
dz=zvals[1]-zvals[0]
dvdz=dh^3*dcomvoldz(zvals,omega0,omegal0)
volume=total(marea)*total(dvdz*dz,/double)

; get initial guess at mean density
sum=total(weight/sel,/double)
sumsqr=total((weight/sel)^2,/double)
n1=sum/volume
n1err=sqrt(sumsqr)/volume
n1=n1*corrfact
n1err=n1err*corrfact
splog,'initial: '+string(n1,format='(e11.4)')+' +/- '+ $
  string(n1err,format='(e11.4)')+' per Mpc^3 per steradian'

if(NOT keyword_set(simple)) then begin
    n1old=0.
    n1err=0.

    while(abs(n1old-n1) gt 0.002*n1err) do begin
        
; get weighted volume for each special case
        wvall=fltarr(nm)
        for j=0L, nm-1L do begin
; get weighted volume 
            dwvdz=total(marea)*dvdz*selfunc[*,j]/(1.+n1*j3*selfunc[*,j])
            wvall[j]=total(dwvdz*dz,/double)
        endfor
        
; now average the weighted volumes according to the fractional area
; within each set of flux limits
        wv=float(total(wvall*fraclimits,/double)/total(fraclimits,/double))
        
; add weighted number of galaxies to get new n1
        contrib=1./((1.+n1*j3*sel)*wv)
        n1old=n1
        n1=total(weight*contrib,/double)
        n1err=sqrt(total((weight*contrib)^2/wv,/double))
        n1=n1*corrfact
        n1err=n1err*corrfact
        
; output
        splog,'iterate: '+string(n1,format='(e11.4)')+' +/- '+ $
          string(n1err,format='(e11.4)')+' per Mpc^3 per steradian'
    end
endif

end
