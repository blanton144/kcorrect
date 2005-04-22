;+
; NAME:
;   galex_to_maggies
; PURPOSE:
;   convert GALEX catalog input to Galactic-extcintion corrected AB maggies 
; CALLING SEQUENCE:
;   galex_to_maggies,galex,maggies,ivar
; INPUTS:
;   galex - [N] GALEX "mcat" style input, with:
;            .ALPHA_J2000
;            .DELTA_J2000
;            .NUV_MAG
;            .NUV_MAGERR
;            .FUV_MAG
;            .FUV_MAGERR
; OUTPUTS:
;   maggies - [2, N] output in AB maggies in FUV and NUV filters
;   ivar - [2, N] inverse variance of maggies
; COMMENTS:
;   It also ALWAYS applies a minimum error of [0.05, 0.05] in [FUV, NUV]
;
;   If the GALEX structure has a .E_BV entry, we use that for the
;   GALEX extinction, but we call dust_getval() to read the SFD maps
;   otherwise.
; REVISION HISTORY:
;   07-Apr-2005  Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro galex_to_maggies, galex, maggies, ivar

common com_galex_to_maggies, nuv_leff, fuv_leff, extvoebv, nuv_extoextv, $
  fuv_extoextv

minfuverror=0.05
minnuverror=0.05

maggies=fltarr(2,n_elements(galex))
ivar=fltarr(2,n_elements(galex))

if(n_elements(nuv_leff) eq 0) then begin
    nuv_leff=k_lambda_eff(filterlist='galex_NUV.par')
    fuv_leff=k_lambda_eff(filterlist='galex_FUV.par')
    extvoebv=3.10
    nuv_extoextv=(ext_ccm(nuv_leff))[0]
    fuv_extoextv=(ext_ccm(fuv_leff))[0]
endif

iebv=tag_indx(galex[0], 'e_bv')
if(iebv eq -1) then begin
    glactc, galex.alpha_j2000, galex.delta_j2000, 2000., gl, gb, 1, /deg
    ebv=dust_getval(gl,gb,/interp,/noloop)
endif else begin
    ebv=galex.e_bv
endelse
nuv_extinction=ebv*nuv_extoextv*extvoebv
fuv_extinction=ebv*fuv_extoextv*extvoebv

maggies[0,*]=0.
ivar[0,*]=0.
igood=where(galex.fuv_mag ne -999. and galex.fuv_mag ne -99., ngood)
if(ngood gt 0) then begin
    maggies[0,igood]=10.^(-0.4*(galex[igood].fuv_mag-fuv_extinction[igood]))
    ivar[0,igood]= $
      1./((galex[igood].fuv_magerr^2+minfuverror^2)*(0.4*alog(10.)* $
                                                     maggies[0,igood])^2)
endif 

maggies[1,*]=0.
ivar[1,*]=0.
igood=where(galex.nuv_mag ne -999. and galex.fuv_mag ne -99., ngood)
if(ngood gt 0) then begin
    maggies[1,igood]=10.^(-0.4*(galex[igood].nuv_mag-nuv_extinction[igood]))
    ivar[1,igood]= $
      1./((galex[igood].nuv_magerr^2+minnuverror^2)*(0.4*alog(10.)* $
                                                     maggies[1,igood])^2)
endif

end
