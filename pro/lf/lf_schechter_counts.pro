;+
; NAME:
;   lf_schechter_counts
; PURPOSE: 
;   Predicts magnitude counts from Schechter fn.
; USAGE:
;   counts=lf_schechter_counts(mmin, mmax, schechter, q, p, kcorrect, $
;             zvals [, omega0=, omegal0=, zzero=, integrand= ]
; INPUTS:
;   mmin, mmax - flux range to consider
;   schechter, q, p - description of evolving schechter function
;   kcorrect, zvals - grid of typical K-corrections
; OPTIONAL INPUTS:
;   zzero - reference redshift for q and p evolution (default 0.1)
;   omega0 - matter density (default 0.3)
;   omegal0 - vacuum energy density (default 0.7)
; OPTIONAL OUTPUTS:
;   integrand - redshift counts
; REVISION HISTORY:
;   2002-8-25  written - Blanton
;-
function lf_schechter_counts,mmin,mmax,schechter,q,p,kcorrect,zvals, $
                             omega0=omega0,omegal0=omegal0,zzero=zzero, $
                             integrand=integrand

; settings
pi=3.14159265358979D
dh=2.99792D+5/100.D

; defaults
if(n_elements(omega0) eq 0) then omega0=0.3
if(n_elements(omegal0) eq 0) then omegal0=0.7
if(n_elements(zzero) eq 0) then zzero=0.1
if(NOT keyword_set(ambin)) then ambin=0.02
nzvals=n_elements(zvals)

; get comoving volume and distance modulus
dc=dh*comdis2(zvals,omega0,omegal0) 
dvdz=dh^3*dcomvoldz(zvals,omega0,omegal0)
dl=(1.+zvals)*dc
dm=5.*alog10(dl)+25.

; get the density at each redshift
mlimits=dblarr(2,nzvals)
mlimits[0,*]=mmin-dm-kcorrect+q*(zvals-zzero)
mlimits[1,*]=mmax-dm-kcorrect+q*(zvals-zzero)
integrand=dblarr(nzvals)
for i=0, nzvals-1 do begin
    nbin=long((mlimits[1,i]-mlimits[0,i])/ambin)
    am=mlimits[0,i]+(mlimits[1,i]-mlimits[0,i])*(dindgen(nbin)+0.5)/ $
      double(nbin)
    dm=am[1]-am[0]
    sch=lf_schechter(am,schechter.phistar,schechter.mstar,schechter.alpha)
    integrand[i]=dvdz[i]*10.^(0.4*p*(zvals[i]-zzero))*total(dm*sch,/double)
endfor

; get the counts
dz=zvals[1]-zvals[0]
counts=total(integrand*dz,/double)

return,counts

end
