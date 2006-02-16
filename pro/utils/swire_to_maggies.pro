;+
; NAME:
;   swire_to_maggies
; PURPOSE:
;   convert SWIRE input to Galactic-extinction corrected AB maggies 
; CALLING SEQUENCE:
;   swire_to_maggies, swire, maggies, ivar
; INPUTS:
;   swire - [N] SWIRE style input, with:
;               .RA
;               .DEC
;               .FLUX_AP_36
;               .UNCF_AP_36
;               .FLUX_AP_45
;               .UNCF_AP_45
;               .FLUX_AP_58
;               .UNCF_AP_58
;               .FLUX_AP_80
;               .UNCF_AP_80
;               .FLUX_AP_24
;               .UNCF_AP_24
;               .AP_M_U
;               .MSIG_U
;               .AP_M_G
;               .MSIG_G
;               .AP_M_R
;               .MSIG_R
;               .AP_M_I
;               .MSIG_I
;               .AP_M_Z
;               .MSIG_Z
; OUTPUTS:
;   maggies - [10, N] output in AB maggies in
;             ugriz[3.6][4.5][5.8][8.0][24] bands
;   ivar - [10, N] inverse variance of maggies
; COMMENTS:
;   It ALWAYS applies a minimum error of 0.02 (in magnitude) in all bands
; REVISION HISTORY:
;   07-Apr-2005  Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro swire_to_maggies, swire, maggies, ivar

red_fac=[5.08248, 3.73964, 2.71230, 2.05665, 1.45819, 0., 0., 0., 0., 0.]
abfix=[0., 0.0, 0., 0., 0., 0., 0. ,0., 0. ,0.]

maggies=fltarr(10, n_elements(swire))
ivar=fltarr(10, n_elements(swire))

glactc, swire.ra, swire.dec, 2000., gl, gb, 1, /deg
ebv=dust_getval(gl,gb, /interp, /noloop)

i=0
ftag=tag_indx(swire[0], 'ap_m_u')
utag=tag_indx(swire[0], 'msig_u')
ig=where(swire.(utag) ne -99. AND swire.(utag) ne 0. AND $
         swire.(ftag) ne -99. AND swire.(ftag) ne 0., ng)
if(ng gt 0) then begin
    maggies[i,ig]= 10.^(-0.4*swire[ig].(ftag))*10.^(0.4*ebv[ig]*red_fac[i])
    ivar[i,ig]= 1./(0.4*alog(10.)*maggies[i,ig]*swire[ig].(utag))^2
endif

i=1
ftag=tag_indx(swire[0], 'ap_m_g')
utag=tag_indx(swire[0], 'msig_g')
ig=where(swire.(utag) ne -99. AND swire.(utag) ne 0. AND $
         swire.(ftag) ne -99. AND swire.(ftag) ne 0., ng)
if(ng gt 0) then begin
    maggies[i,ig]= 10.^(-0.4*swire[ig].(ftag))*10.^(0.4*ebv[ig]*red_fac[i])
    ivar[i,ig]= 1./(0.4*alog(10.)*maggies[i,ig]*swire[ig].(utag))^2
endif

i=2
ftag=tag_indx(swire[0], 'ap_m_r')
utag=tag_indx(swire[0], 'msig_r')
ig=where(swire.(utag) ne -99. AND swire.(utag) ne 0. AND $
         swire.(ftag) ne -99. AND swire.(ftag) ne 0., ng)
if(ng gt 0) then begin
    maggies[i,ig]= 10.^(-0.4*swire[ig].(ftag))*10.^(0.4*ebv[ig]*red_fac[i])
    ivar[i,ig]= 1./(0.4*alog(10.)*maggies[i,ig]*swire[ig].(utag))^2
endif

i=3
ftag=tag_indx(swire[0], 'ap_m_i')
utag=tag_indx(swire[0], 'msig_i')
ig=where(swire.(utag) ne -99. AND swire.(utag) ne 0. AND $
         swire.(ftag) ne -99. AND swire.(ftag) ne 0., ng)
if(ng gt 0) then begin
    maggies[i,ig]= 10.^(-0.4*swire[ig].(ftag))*10.^(0.4*ebv[ig]*red_fac[i])
    ivar[i,ig]= 1./(0.4*alog(10.)*maggies[i,ig]*swire[ig].(utag))^2
endif

i=4
ftag=tag_indx(swire[0], 'ap_m_z')
utag=tag_indx(swire[0], 'msig_z')
ig=where(swire.(utag) ne -99. AND swire.(utag) ne 0. AND $
         swire.(ftag) ne -99. AND swire.(ftag) ne 0., ng)
if(ng gt 0) then begin
    maggies[i,ig]= 10.^(-0.4*swire[ig].(ftag))*10.^(0.4*ebv[ig]*red_fac[i])
    ivar[i,ig]= 1./(0.4*alog(10.)*maggies[i,ig]*swire[ig].(utag))^2
endif

i=5
ftag=tag_indx(swire[0], 'flux_ap5_36')
utag=tag_indx(swire[0], 'uncf_ap5_36')
ig=where(swire.(utag) ne -99., ng)
if(ng gt 0) then begin
    maggies[i,ig]= swire[ig].(ftag)*1.e-6/3631.*10.^(0.4*ebv[ig]*red_fac[i])
    ivar[i,ig]= 1./(swire[ig].(utag)*1.e-6/3631.)^2
endif

i=6
ftag=tag_indx(swire[0], 'flux_ap5_45')
utag=tag_indx(swire[0], 'uncf_ap5_45')
ig=where(swire.(utag) ne -99., ng)
if(ng gt 0) then begin
    maggies[i,ig]= swire[ig].(ftag)*1.e-6/3631.*10.^(0.4*ebv[ig]*red_fac[i])
    ivar[i,ig]= 1./(swire[ig].(utag)*1.e-6/3631.)^2
endif

i=7
ftag=tag_indx(swire[0], 'flux_ap5_58')
utag=tag_indx(swire[0], 'uncf_ap5_58')
ig=where(swire.(utag) ne -99., ng)
if(ng gt 0) then begin
    maggies[i,ig]= swire[ig].(ftag)*1.e-6/3631.*10.^(0.4*ebv[ig]*red_fac[i])
    ivar[i,ig]= 1./(swire[ig].(utag)*1.e-6/3631.)^2
endif

i=8
ftag=tag_indx(swire[0], 'flux_ap5_80')
utag=tag_indx(swire[0], 'uncf_ap5_80')
ig=where(swire.(utag) ne -99., ng)
if(ng gt 0) then begin
    maggies[i,ig]= swire[ig].(ftag)*1.e-6/3631.*10.^(0.4*ebv[ig]*red_fac[i])
    ivar[i,ig]= 1./(swire[ig].(utag)*1.e-6/3631.)^2
endif

i=9
ftag=tag_indx(swire[0], 'flux_ap5_24')
utag=tag_indx(swire[0], 'uncf_ap5_24')
ig=where(swire.(utag) ne -99., ng)
if(ng gt 0) then begin
    maggies[i,ig]= swire[ig].(ftag)*1.e-6/3631.*10.^(0.4*ebv[ig]*red_fac[i])
    ivar[i,ig]= 1./(swire[ig].(utag)*1.e-6/3631.)^2
endif

k_minerror, maggies, ivar, replicate(0.02, 10)
k_abfix, maggies, ivar, aboff=abfix

end
