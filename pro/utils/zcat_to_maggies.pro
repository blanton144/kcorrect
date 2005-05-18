;+
; NAME:
;   zcat_to_maggies
; PURPOSE:
;   convert DEEP zcat input to Galactic-extinction corrected AB maggies 
; CALLING SEQUENCE:
;   zcat_to_maggies, zcat, maggies, ivar
; INPUTS:
;   zcat - [N] DEEP style input, with:
;               .MAGB
;               .MAGR
;               .MAGI
;               .SFD_EBV
; OUTPUTS:
;   maggies - [3, N] output in AB maggies in BRI
;   ivar - [3, N] inverse variance of maggies
; COMMENTS:
;   Calculates an approximate error based on the flux.
; 
;   It ALWAYS applies a minimum error of [0.02, 0.02, 0.02] in
;   BRI.
; REVISION HISTORY:
;   07-Apr-2005  Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro zcat_to_maggies, zcat, maggies, ivar

red_fac=[4.32, 2.63, 1.96]

maggies=fltarr(3, n_elements(zcat))
ivar=fltarr(3, n_elements(zcat))

mbase=24.0
sigmbase=0.05

i=0
extinction=zcat.sfd_ebv*red_fac[i]
maggies[i,*]=10.^(-0.4*(zcat.magb-extinction))
ivar[i,*]=1./(0.02^2+(sigmbase*10.^(-0.4*(mbase-zcat.magb)))^2)/ $
  (0.4*alog(10.))^2/maggies[i,*]^2

i=1
extinction=zcat.sfd_ebv*red_fac[i]
maggies[i,*]=10.^(-0.4*(zcat.magr-extinction))
ivar[i,*]=1./(0.02^2+(sigmbase*10.^(-0.4*(mbase-zcat.magr)))^2)/ $
  (0.4*alog(10.))^2/maggies[i,*]^2

i=2
extinction=zcat.sfd_ebv*red_fac[i]
maggies[i,*]=10.^(-0.4*(zcat.magi-extinction))
ivar[i,*]=1./(0.02^2+(sigmbase*10.^(-0.4*(mbase-zcat.magi)))^2)/ $
  (0.4*alog(10.))^2/maggies[i,*]^2

end
