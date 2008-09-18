;+
; NAME:
;   deep_to_maggies
; PURPOSE:
;   convert DEEP zcat input to Galactic-extinction corrected AB maggies 
; CALLING SEQUENCE:
;   deep_to_maggies, zcat, maggies, ivar
; INPUTS:
;   zcat - [N] DEEP style input, with:
;               .MAGB
;               .MAGR
;               .MAGI
; OUTPUTS:
;   maggies - [3, N] output in AB maggies in BRI
;   ivar - [3, N] inverse variance of maggies
; COMMENTS:
;   Calculates an approximate error based on the flux.
; 
;   It ALWAYS applies a minimum error of [0.02, 0.02, 0.02] in
;   BRI.
;
;   The v4_0 version incorrectly applied SFD extinction, which is
;   already applied to the DEEP magnitudes.
; REVISION HISTORY:
;   07-Apr-2005  Mike Blanton, NYU
;   22-July-2008 Guangtun Zhu, NYU, Initiate ivar with 1E-32
;   17-Sep-2008 Guangtun Zhu, NYU, Initiate ivar with 0.
;-
;------------------------------------------------------------------------------
pro deep_to_maggies, zcat, maggies, ivar

maggies=fltarr(3, n_elements(zcat))
ivar=fltarr(3, n_elements(zcat))

mbase=24.0
sigmbase=0.05

i=0
maggies[i,*]=10.^(-0.4*(zcat.magb))
ivar[i,*]=1./(0.02^2+(sigmbase*10.^(-0.4*(mbase-zcat.magb)))^2)/ $
  (0.4*alog(10.))^2/maggies[i,*]^2

i=1
maggies[i,*]=10.^(-0.4*(zcat.magr))
ivar[i,*]=1./(0.02^2+(sigmbase*10.^(-0.4*(mbase-zcat.magr)))^2)/ $
  (0.4*alog(10.))^2/maggies[i,*]^2

i=2
maggies[i,*]=10.^(-0.4*(zcat.magi))
ivar[i,*]=1./(0.02^2+(sigmbase*10.^(-0.4*(mbase-zcat.magi)))^2)/ $
  (0.4*alog(10.))^2/maggies[i,*]^2

end
