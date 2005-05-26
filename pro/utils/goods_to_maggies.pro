;+
; NAME:
;   goods_to_maggies
; PURPOSE:
;   convert GOODS catalog input to Galactic-extcintion corrected AB maggies 
; CALLING SEQUENCE:
;   goods_to_maggies,goods,maggies,ivar
; INPUTS:
;   goods - [N] GOODS style input
;               .RA (J2000 degrees)
;               .DEC (J2000 degrees)
;               .BMAG_MAGAUTO
;               .BMAGERR_MAGAUTO
;               .VMAG_MAGAUTO
;               .VMAGERR_MAGAUTO
;               .IMAG_MAGAUTO
;               .IMAGERR_MAGAUTO
;               .ZMAG_MAGAUTO
;               .ZMAGERR_MAGAUTO
;               .JMAG_MAGAUTO
;               .JMAGERR_MAGAUTO
;               .HMAG_MAGAUTO
;               .HMAGERR_MAGAUTO
;               .KMAG_MAGAUTO
;               .KMAGERR_MAGAUTO
; OUTPUTS:
;   maggies - [7, N] output in AB maggies in BVizJHK
;   ivar - [7, N] inverse variance of maggies
; COMMENTS:
;   It ALWAYS applies a minimum error of [0.02, 0.02, 0.02] in
;   all bandpasses
;
;   Requires you to have the dust maps so that dust_getval can find
;   them. (If somebody wants me to set "default" columns in the
;   goods structure that this code looks for, let me know).
; REVISION HISTORY:
;   07-Apr-2005  Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro goods_to_maggies, goods, maggies, ivar

minerrors=replicate(0.02, 7)
dfactors=[4.32, 3.32, 2.00, 1.54, 0.90, 0.58, 0.37]
names=['B','V','I','Z','J','H','K']

;; get Galactic extinction
euler,goods.ra,goods.dec,ll,bb,1
ebv=dust_getval(ll, bb, /interp, /noloop)

maggies=fltarr(7,n_elements(goods))
ivar=fltarr(7,n_elements(goods))
for iband=0L, 6L do begin
    itag_m=tag_indx(goods[0], names[iband]+'mag_magauto')
    itag_msig=tag_indx(goods[0], names[iband]+'magerr_magauto')
    indx=where(goods.(itag_m) gt 0. AND $
               goods.(itag_msig) gt 0., count)
    if(count gt 0) then begin
        maggies[iband, indx]=10.^(-0.4*(goods[indx].(itag_m) $
                                        -ebv[indx]*dfactors[iband]))
        ivar[iband, indx]= $
          1./(0.4*alog(10.)*maggies[iband,indx]* $
              goods[indx].(itag_msig))^2
    endif
endfor

k_minerror, maggies, ivar, minerrors

end
