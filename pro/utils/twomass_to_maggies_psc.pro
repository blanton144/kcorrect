;+
; NAME:
;   twomass_to_maggies_psc
; PURPOSE:
;   convert 2MASS catalog input to Galactic-extcintion corrected AB maggies 
; CALLING SEQUENCE:
;   twomass_to_maggies,twomass,maggies,ivar
; INPUTS:
;   twomass - [N] 2MASS style input, with:
;               .RA (J2000 degrees)
;               .DECL (J2000 degrees)
;               .J_M_EXT
;               .J_MSIG_EXT
;               .J_FLG_EXT
;               .H_M_EXT
;               .H_MSIG_EXT
;               .H_FLG_EXT
;               .K_M_EXT
;               .K_MSIG_EXT
;               .K_FLG_EXT
; OUTPUTS:
;   maggies - [2, N] output in AB maggies in FUV and NUV filters
;   ivar - [2, N] inverse variance of maggies
; COMMENTS:
;   Uses k_vega2ab() to get correction from 2MASS Vega magnitudes to
;   2MASS AB magnitudes. As of 2005-04-18, this yields [0.91, 1.39,
;   1.85] for AB-Vega in JHKs respectively.
;
;   It ALWAYS applies a minimum error of [0.03, 0.03, 0.03] in
;   JHKs
;
;   Requires you to have the dust maps so that dust_getval can find
;   them. (If somebody wants me to set "default" columns in the
;   twomass structure that this code looks for, let me know).
; 
;   **This routine is identical to TWOMASS_TO_MAGGIES** but it works
;   with the PSC.
;
; REVISION HISTORY:
;   07-Apr-2005  Mike Blanton, NYU
;   03-Feb-2010  J. Moustakas, UCSD
;-
;------------------------------------------------------------------------------
pro twomass_to_maggies_psc, twomass, maggies, ivar

common com_g2m, vega2ab

minerrors=[0.03, 0.03, 0.03]
red_fac=[0.902, 0.576, 0.367]
names=['j','h','k']

;; get AB corrections
filterlist=['twomass_J.par', 'twomass_H.par', 'twomass_Ks.par']
if(n_elements(vega2ab) eq 0) then $
  vega2ab=k_vega2ab(filterlist=filterlist, /kurucz)

;; get Galactic extinction
euler,twomass.ra,twomass.dec,ll,bb,1
ebv=dust_getval(ll, bb, /interp, /noloop)

maggies=fltarr(3,n_elements(twomass))
ivar=fltarr(3,n_elements(twomass))
for iband=0L, 2L do begin
; http://www.ipac.caltech.edu/2mass/releases/second/doc/ancillary/pscformat.html
    itag_m=tag_indx(twomass[0], names[iband]+'_m') ; "default" magnitude
    itag_msig=tag_indx(twomass[0], names[iband]+'_msigcom') ; "total" magnitude error
    indx=where(twomass.(itag_m) gt 0. AND $
               twomass.(itag_msig) gt 0., count)
    if(count gt 0) then begin
        maggies[iband, indx]=10.^(-0.4*(twomass[indx].(itag_m) $
                                        -ebv[indx]*red_fac[iband] $
                                        +vega2ab[iband]))
        ivar[iband, indx]= $
          1./((0.4*alog(10.)*maggies[iband,indx])^2* $
              (twomass[indx].(itag_msig)^2+minerrors[iband]^2))
    endif
endfor

end
