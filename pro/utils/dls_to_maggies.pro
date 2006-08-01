;+
; NAME:
;   dls_to_maggies
; PURPOSE:
;   convert DLS input to Galactic-extinction corrected AB maggies 
; CALLING SEQUENCE:
;   dls_to_maggies, dls, maggies, ivar
; INPUTS:
;   dls - [N] DLS file (for ?=BVRz):
;               .RA
;               .DEC
;               .MAG_AUTO? 
;               .MAGERR_AUTO?
; OUTPUTS:
;   maggies - [3, N] output in AB maggies in BRI
;   ivar - [3, N] inverse variance of maggies
; COMMENTS:
;   It ALWAYS applies a minimum error of [0.02, 0.02, 0.02, 0.02] in
;   BVRz.
;
;   The v4_0 version incorrectly applied SFD extinction, which is
;   already applied to the DEEP magnitudes.
; REVISION HISTORY:
;   07-Apr-2005  Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro dls_to_maggies, dls, maggies, ivar

euler,dls.ra,dls.dec,ll,bb,1
red_fac=[1.32, 0.99, 0.81, 0.472]*3.24
extinction= red_fac # dust_getval(ll, bb, /interp, /noloop)

maggies=fltarr(4, n_elements(dls))
ivar=fltarr(4, n_elements(dls))

i=0
maggies[i,*]=10.^(-0.4*(dls.mag_autob-extinction[i,*]+ $
                        k_vega2ab(filterlist='bessell_B.par', /kur)))
ivar[i,*]=1./(0.02^2+(dls.magerr_autob)^2)/ $
  (0.4*alog(10.))^2/maggies[i,*]^2

i=1
maggies[i,*]=10.^(-0.4*(dls.mag_autov-extinction[i,*]+ $
                        k_vega2ab(filterlist='bessell_V.par', /kur)))
ivar[i,*]=1./(0.02^2+(dls.magerr_autov)^2)/ $
  (0.4*alog(10.))^2/maggies[i,*]^2

i=2
maggies[i,*]=10.^(-0.4*(dls.mag_autor-extinction[i,*]+ $
                        k_vega2ab(filterlist='bessell_R.par', /kur)))
ivar[i,*]=1./(0.02^2+(dls.magerr_autor)^2)/ $
  (0.4*alog(10.))^2/maggies[i,*]^2

i=3
maggies[i,*]=10.^(-0.4*(dls.mag_autoz-extinction[i,*]+ $
                        k_vega2ab(filterlist='sdss_z0.par', /kur)))
ivar[i,*]=1./(0.02^2+(dls.magerr_autoz)^2)/ $
  (0.4*alog(10.))^2/maggies[i,*]^2


end
