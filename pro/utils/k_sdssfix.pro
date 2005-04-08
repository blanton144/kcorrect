;+
; NAME:
;   k_sdssfix
; PURPOSE:
;   Take SDSS database asinh magnitudes and errors and "fixes" them
; CALLING SEQUENCE:
;   k_sdssfix, mags, mags_err, maggies, maggies_ivar [, /standard, aboff=]
; INPUTS:
;   mags - input asinh magnitudes (luptitudes)
;   mags_err - input uncertainties in mags
; OUTPUTS:
;   maggies - best AB maggies to use 
;   maggies_ivar - best inverse variance to use
; KEYWORDS:
;   /standard - assume standard magnitudes, not luptitudes
; COMMENTS:
;   This converts from SDSS database asinh magnitudes to AB maggies. 
;
;   It "fixes" errors in the sense that for "bad" measurements or
;   errors you assign values which are not absurd. Not necessary for
;   Princeton-style input (which sets ivar to zero as
;   appropriate). 
;  
;   Also adds errors in quadrature:
;     sigma(ugriz) = [0.05, 0.02, 0.02, 0.02, 0.03]
;   to account for calibration uncertainties.
; BUGS:
;   Needs better tracking of Eisenstein numbers
; REVISION HISTORY:
;   07-Feb-2002  Written by Mike Blanton, NYU
;   03-Jun-2003  Updated to v3_0 by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_sdssfix, mags, mags_err, maggies, maggies_ivar, standard=standard

bvalues=[1.4D-10, 0.9D-10, 1.2D-10, 1.8D-10, 7.4D-10]
errband=[0.05,0.02,0.02,0.02,0.03]

if((size(mags))[0] eq 1) then begin
    nk=5
    ngalaxy=1
endif else begin
    nk=(size(mags,/dimens))[0]
    ngalaxy=(size(mags,/dimens))[01]
endelse

if(nk ne 5) then begin
  klog, 'k_sdssfix for 5-band SDSS observations only! ignoring ...'
  return 
endif

mags_ivar=reform(k_sdss_err2ivar(mags_err),nk,ngalaxy)

indx=where(mags eq -9999 and mags_ivar ne 0.,count)
if(count gt 0) then begin
    klog, string(count)+' cases of good error but bad value in k_sdssfix!'
    mags_ivar[indx]=0.
endif

maggies=dblarr(nk,ngalaxy)
maggies_ivar=dblarr(nk,ngalaxy)
for k=0L, nk-1L do begin
    indx=where(mags_ivar[k,*] ne 0.,count)
    if(count gt 0) then begin
        if(NOT keyword_set(standard)) then begin
            err=1./sqrt(mags_ivar[k,indx])
            maggies[k,indx]=k_lups2maggies(mags[k, indx], err, $
                                           maggies_err=merr, $
                                           bvalues=bvalues[k])
            maggies_ivar[k,indx]=1./(merr^2)
        endif else begin
            maggies[k,indx]=(10.D)^(-(0.4D)*(mags[k,indx]))
            maggies_ivar[k,indx]= $
              mags_ivar[k,indx]/(0.4*alog(10.)*maggies[k,indx])^2.
        endelse
    endif
endfor
indx=where(mags_ivar eq 0.,count)
if(count gt 0) then begin
   maggies[indx]=0.
   maggies_ivar[indx]=0.
endif

k_minerror, maggies, maggies_ivar, errband
k_abfix, maggies, maggies_ivar, aboff=aboff

end
;------------------------------------------------------------------------------
