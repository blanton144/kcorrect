;+
; NAME:
;   kcorrect_so_ext
; PURPOSE:
;   returns appropriate dynamic library extension given arch
; CALLING SEQUENCE:
;   so_ext= kcorrect_so_ext()
; COMMENTS:
;   necessary to deal with non-standard .dylib extension on darwin
; REVISION HISTORY:
;   20-Feb-2004  Written by M. Blanton, NYU
;-
;------------------------------------------------------------------------------
function kcorrect_so_ext

common com_kcorrect_so_ext, so_ext

if(NOT keyword_set(so_ext)) then begin
    spawn,'uname',uname
    uname=uname[0]
    so_ext='so'
    if(uname eq 'HP-UX') then so_ext='sl'
    if(uname eq 'Darwin') then so_ext='dylib'
endif

return, so_ext
end
