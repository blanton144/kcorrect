;+
; NAME:
;   lf_vmax
; PURPOSE:
;   calculate luminosity function using the simple vmax method
; USAGE:
;   lf_vmax,absmag,kcorrect,ikcorrect,uniqkcorrect, $
;     nzvals,mmin,mmax,sample_absmmin,sample_absmmax, $
;     sample_zmin,sample_zmax, absmk, phi, phierr [, omega0=, $
;     omegal0=, vmax=, imeankcorrect= ]
; INPUTS:
;   absmag           [N] absolute magnitudes
;   kcorrect         [N] K-correction for each galaxy
;   ikcorrect        [N] index of K-correction function for each gal
;   uniqkcorrect     [nzvals,N_k] K(z) for each K-correction function
;   nzvals           number of redshift slices in each K-correction fn
;   mmin             [N] minimum apparent mag for each galaxy
;   mmax             [N] maximum apparent mag for each galaxy
;   sample_absmmin   observed absolute magnitude minimum of sample
;   sample_absmmax   observed absolute magnitude maximum of sample
;   sample_zmin      minimum redshift of sample
;   sample_zmax      maximum redshift of sample
; OPTIONAL INPUTS:
;   omega0           omega_matter to use (default: 0.3)
;   omegal0          omega_lambda to use (default: 0.7)
;   imeankcorrect    mean index for K-correction to use for
;                    uncertainties in 0-bins (default mean of given
;                    ikcorrect)
; OUTPUTS:
;   absmk            center of each abs mag bin
;   phi              density in each abs mag bin, in # per cubic Mpc per mag
;   phierr           error in phi
; OPTIONAL INPUTS/OUTPUTS:
;   vmax             vmax used for each object
; REVISION HISTORY:
;   2002-11-17  written - Blanton
;-
pro lf_vmax,absmag,ikcorrect,uniqkcorrect, $
            nzvalskcorrect,mmin,mmax,sample_absmmin,sample_absmmax, nabsmk, $
            sample_zmin,sample_zmax,area,absmk,phi,phierr,omega0=omega0, $
            omegal0=omegal0, vmax=vmax, imeankcorrect=imeankcorrect

; settings
pi=3.14159265358979D
dh=2.99792D+5/100.D

; defaults
if(n_elements(omega0) eq 0) then omega0=0.3
if(n_elements(omegal0) eq 0) then omegal0=0.7
if(n_elements(imeankcorrect) eq 0) then $
  imeankcorrect=long(0.5*(min(ikcorrect)+max(ikcorrect)))

; calculate vmax for each object (if not input)
if(n_elements(vmax) ne n_elements(absmag)) then begin
    lf_calc_vmax,absmag,ikcorrect,uniqkcorrect,nzvalskcorrect, $
      mmin,mmax,sample_zmin,sample_zmax,area,vmax,omega0=omega0,omegal0=omegal0
endif

; check for zeros 
indx_nzero=where(vmax gt 0.,count_nzero)
if(count_nzero ne n_elements(vmax)) then $
  splog, 'WARNING: '+string(n_elements(vmax)-count_nzero)+ $
  ' objects with vmax <= 0, ignoring'

; sum, while binning
absmlimits=sample_absmmin+(dindgen(nabsmk+1L))* $
  (sample_absmmax-sample_absmmin)/(double(nabsmk))
absmk=0.5*(absmlimits[0L:nabsmk-2L]+absmlimits[1L:nabsmk-1L])
phi=dblarr(nabsmk)
phierr=dblarr(nabsmk)
for i=0L, n_elements(absmk)-1L do begin
    indx_in_bin=where(absmlimits[i+0] le absmag[indx_nzero] and $
                      absmlimits[i+1] gt absmag[indx_nzero], count_in_bin)
    dabsm=absmlimits[i+1]-absmlimits[i+0]
    if(count_in_bin gt 0) then begin
        phi[i]=total(1./vmax[indx_nzero[indx_in_bin]],/double)/dabsm
        phierr[i]=sqrt(total(1./vmax[indx_nzero[indx_in_bin]]^2,/double))/dabsm
    endif else begin
        phi[i]=0./dabsm
        lf_calc_vmax,absmk[i],imeankcorrect,uniqkcorrect,nzvalskcorrect, $
          mmin[i],mmax[i],sample_zmin,sample_zmax,area,vmax_ul,omega0=omega0, $
          omegal0=omegal0
        phierr[i]=1./vmax_ul/dabsm
    endelse 
endfor

end
