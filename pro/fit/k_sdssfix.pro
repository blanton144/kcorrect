;+
; NAME:
;   k_sdssfix
; PURPOSE:
;   Take a set of SDSS pipeline magnitudes and errors and "fixes" them, in the
;   sense that for "bad" measurements or errors you assign values
;   which are not absurd. Not necessary for Princeton-style input
;   (which sets ivar to zero as appropriate). 
; CALLING SEQUENCE:
;   k_sdssfix, mags, mags_stddev, maggies, maggies_ivar 
; INPUTS:
;   mags - input magnitudes
;   mags_stddev - input uncertainties in magnitude
; OPTIONAL INPUTS:
; OUTPUTS:
;   maggies - best maggies to use
;   maggies_ivar - best inverse variance to use
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   07-Feb-2002  Written by Mike Blanton, NYU
;   03-Jun-2003  Updated to v3_0 by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function k_sdss_err2ivar, err, verbose=verbose

errband=[0.05,0.02,0.02,0.02,0.03]

ivar=fltarr(5,n_elements(err)/5L)+1.E

error_indx=where(err eq -9999, error_count) 
if(error_count gt 0) then ivar[error_indx]=0.E

error_indx=where(err eq -1000, error_count) 
if(error_count gt 0) then ivar[error_indx]=0.E

error_indx=where(err eq 0, error_count) 
if(error_count gt 0) then ivar[error_indx]=0.E

if(keyword_set(verbose)) then begin
    unusual_indx=where(ivar eq 1. and err lt 0,unusual_count)
    if(unusual_count gt 0) then begin
        for i=0L,unusual_count-1L do $
          splog,'Unusual ERR value of '+string(err[unusual_indx[i]])+ $
          '; setting IVAR to zero'
        ivar[unusual_indx]=0.E
    endif
endif

q=err gt 0
for i=0L, 4L do $
  ivar[i,*]=ivar[i,*]*q/((err[i,*]^2+errband[i]^2)*q+(1-q))

return, ivar

end
;
pro k_sdssfix, mags, mags_stddev, maggies, maggies_ivar

if((size(mags))[0] eq 1) then begin
    nk=1
    ngalaxy=n_elements(mags)
endif else begin
    nk=(size(mags,/dimens))[0]
    ngalaxy=(size(mags,/dimens))[01]
endelse

if(nk ne 5) then begin
  klog, 'k_sdssfix for 5-band SDSS observations only! ignoring ...'
  return 
endif

mags_ivar=reform(k_sdss_err2ivar(mags_stddev),nk,ngalaxy)

indx=where(mags eq -9999 and mags_ivar ne 0.,count)
if(count gt 0) then begin
    klog, string(count)+' cases of good error but bad value in k_sdssfix!'
    mags_ivar[indx]=0.
endif

maggies=dblarr(nk,ngalaxy)
indx=where(mags_ivar ne 0.,count)
if(count gt 0) then begin
   maggies[indx]=(10.D)^(-(0.4D)*mags)
   maggies_ivar[indx]=mags_ivar/(0.4*alog(10.)*maggies)^2.
endif
indx=where(mags_ivar eq 0.,count)
if(count gt 0) then begin
   maggies[indx]=0.
   maggies_ivar[indx]=0.
endif


end
;------------------------------------------------------------------------------