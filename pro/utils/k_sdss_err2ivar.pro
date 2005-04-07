;+
; NAME:
;   k_sdss_err2ivar
; PURPOSE:
;   convert SDSS database magnitude error values to inverse variances
; CALLING SEQUENCE:
;   k_sdss_err2ivar, err [, /verbose]
; INPUTS:
;   err - [5, N] error values
; KEYWORDS:
;   /verbose - loud about possible mistakes
; COMMENTS:
;   This "fixes" errors in the sense that for "bad" measurements or
;   errors you assign values which are not absurd. Not necessary for
;   Princeton-style input (which sets ivar to zero as
;   appropriate). 
; REVISION HISTORY:
;   07-Feb-2002  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function k_sdss_err2ivar, err, verbose=verbose

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

for i=0L, 4L do begin
    q=err[i,*] gt 0.
    ivar[i,*]=ivar[i,*]*q/((err[i,*]^2)*q+(1-q))
endfor

return, ivar

end
;------------------------------------------------------------------------------
