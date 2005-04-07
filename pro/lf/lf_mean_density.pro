;+
; NAME:
;   lf_mean_density
; PURPOSE:
;   calculate mean density given npgauss luminosity function
; USAGE:
;   lf_mean_density, zz, sf, fraclimits, iabsmerr, immin, immax, $
;     ikcorrect, uniqabsmerr, uniqmmin, uniqmmax, uniqkcorrect, $
;     sample_absmmin, sample_absmmax, model_absmmin, model_absmmax, $
;     sample_zmin, sample_zmax, phi, absmk, sigabsmag, q, p, n1, n1err $
;     [, zzero=, nzvals=, omega0=, omegal0=, j3=, /simple, rhozpars=, $
;     rhozfunc=, zvals= ]
; INPUTS:
;   zz - [N] redshifts
;   sf - [N] selection function for each galaxy
;   fraclimits - [Nm] fraction of area in each region with constant
;                flux limit 
;   iabsmerr - [N] index indicating which discrete value of error to
;              use from uniqabsmerr
;   immin - [N] index indicating which discrete value of bright flux
;           limit to use from uniqmmin
;   immax - [N] index indicating which discrete value of faint flux
;           limit to use from uniqmmax
;   ikcorrect - [N] index indicating which kcorrect template to
;               use from uniqkcorrect
;   uniqabsmerr - [Nerr] discrete values of error 
;   uniqmmin - [Nmin] discrete values of bright flux limit
;   uniqmmax - [Nmin] discrete values of faint flux limit
;   uniqkcorrect - [nzvals, Nk] kcorrections as a function of
;                  redshiift for each template
;   sample_absmmin, sample_absmmax - sample limits on absolute mag
;   sample_zmin, sample_zmax - sample limits on redshift
;   phi, absmk, sigabsmag, q, p - description of l.f.
; OPTIONAL INPUTS:
;   nzvals - number of redshift in kcorrect templates (default 50)
;   zzero - reference redshift for q and p evolution (default 0.1)
;   omega0 - matter density (default 0.3)
;   omegal0 - vacuum energy density (default 0.7)
;   rhozpars, rhozfunc - if we have fit out the radial density function
; OPTIONAL KEYWORDS:
;   /simple - do simple weighting 
; OUTPUTS:
;   n1, n1err - mean density and error
; OPTIONAL OUTPUTS:
;   zvals - redshift steps
; BUGS:
; DEPENDENCIES:
; REVISION HISTORY:
;   2002-5-22  written - Blanton
;-
;
pro lf_mean_density,zz,sf,fraclimits,iabsmerr,immin,immax,ikcorrect, $
                    uniqabsmerr, uniqmmin, uniqmmax, uniqkcorrect, $
                    sample_absmmin,sample_absmmax, $
                    model_absmmin,model_absmmax, $
                    sample_zmin, sample_zmax, $
                    phi,absmk,sigabsmag,q,p,n1,n1err, $
                    zzero=zzero,nzvals=nzvals,omega0=omega0, $
                    omegal0=omegal0,j3=j3,simple=simple, $
                    rhozpars=rhozpars,rhozfunc=rhozfunc,zvals=zvals

; defaults
if(n_elements(omega0) eq 0) then omega0=0.3
if(n_elements(omegal0) eq 0) then omegal0=0.7
if(n_elements(nzvals) eq 0) then nzvals=50
if(n_elements(j3) eq 0) then j3=10.^6 ; ie. make it pretty irrelevant

; settings
pi=3.14159265358979D
dh=2.99792D+5/100.D
zvals=sample_zmin+(sample_zmax-sample_zmin)* $
  (dindgen(nzvals)+0.5)/double(nzvals)
ngals=n_elements(zz)
nmmin=n_elements(uniqmmin)
nmmax=n_elements(uniqmmax)
nabsmerr=n_elements(uniqabsmerr)
nkcorrect=n_elements(uniqkcorrect)/nzvals

; get volume of survey
dz=zvals[1]-zvals[0]
dvdz=dh^3*dcomvoldz(zvals,omega0,omegal0)
volume=total(dvdz*dz,/double)

; find correction factor to correct for <rhoz> != 1
corrfact=1.
if(n_elements(rhozpars) gt 0) then begin
    rhoz_int=1.+rhozfunc#rhozpars
    rhoz_zvals=sample_zmin+(sample_zmax-sample_zmin)* $
      (dindgen(n_elements(rhoz_int))+0.5)/double(n_elements(rhoz_int))
    rhoz_dvdz=dh^3*dcomvoldz(rhoz_zvals,omega0,omegal0)
    rhoz_dz=rhoz_zvals[1]-rhoz_zvals[0]
    corrfact=total(rhoz_dvdz*rhoz_dz*rhoz_int,/double)/ $
      total(rhoz_dvdz*rhoz_dz,/double)
endif

; get initial guess at mean density
sum=total(1./sf,/double)
sumsqr=total(1./sf^2,/double)
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
        wvall=dblarr(nabsmerr,nmmin,nmmax,nkcorrect)
        for kabsmerr=0L, nabsmerr-1L do begin
            for kmmin=0L, nmmin-1L do begin
                for kmmax=0L, nmmax-1L do begin
                    for kkcorrect=0L, nkcorrect-1L do begin
; get selection function from limits
                        sfvals= $
                          lf_selfunc_npgauss(zvals,uniqkcorrect[*,kkcorrect], $
                                             uniqmmin[kmmin],uniqmmax[kmmax], $
                                             sample_absmmin,sample_absmmax, $
                                             model_absmmin,model_absmmax, $
                                             phi,absmk, $
                                             sqrt(sigabsmag^2+ $ 
                                                  uniqabsmerr[kabsmerr]), $
                                             q,p,zzero=zzero,omega0=omega0, $
                                             omegal0=omegal0, $
                                             rhozpars=rhozpars, $
                                             rhozfunc=rhozfunc, $
                                             zvals=zvals, $
                                             sample_zmin=sample_zmin, $
                                             sample_zmax=sample_zmax)
                        
; get weighted volume 
                        dwvdz=dvdz*sfvals/(1.+n1*j3*sfvals)
                        wvall[kabsmerr,kmmin,kmmax,kkcorrect]= $
                          total(dwvdz*dz,/double)
                    endfor
                endfor
            endfor
        endfor
        
; now average the weighted volumes according to the fractional area
; within each set of flux limits
        wv=dblarr(nabsmerr,nkcorrect)
        for kabsmerr=0L, nabsmerr-1L do begin
            for kkcorrect=0L, nkcorrect-1L do begin
                wv[kabsmerr,kkcorrect]= $
                  total(wvall[kabsmerr,*,*,kkcorrect]*fraclimits,/double)/ $
                  total(fraclimits,/double)
            endfor
        endfor
        
; add weighted number of galaxies to get new n1
        contrib=1./((1.+n1*j3*sf)*wv[iabsmerr,ikcorrect])
        n1old=n1
        n1=total(contrib,/double)
        n1err=sqrt(total(contrib/wv[iabsmerr,ikcorrect],/double))
        n1=n1*corrfact
        n1err=n1err*corrfact
        
; output
        splog,'iterate: '+string(n1,format='(e11.4)')+' +/- '+ $
          string(n1err,format='(e11.4)')+' per Mpc^3 per steradian'
    end
endif

; we have normalized to the model range in phi, so we have to make 
; sure that n1 applies to that range of phi
totalphi=total(phi*0.5* $
               (errorf((model_absmmax-absmk)/(sqrt(2.)*sigabsmag))- $
                errorf((model_absmmin-absmk)/(sqrt(2.)*sigabsmag))),/double)
n1=n1/totalphi

end
