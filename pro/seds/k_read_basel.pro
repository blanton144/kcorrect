;+
; NAME:
;   k_read_basel
; PURPOSE:
;   read in spectrum from a Basel spectrum file
; CALLING SEQUENCE:
;   k_read_basel, lambda, flux, filename [, teff=, logg=, mh=, vturb=, xh= ]
; INPUTS:
;   filename - file containing Basel format spectrum
; OPTIONAL KEYWORDS:
;   /silent - shut up
; OUTPUTS:
;   lambda - wavelengths at pixel centers
;   flux - flux at each pixel
;   teff, logg, mh, vturb, xh - basel parameters in file
; COMMENTS:
;   1221 elements *hard-coded*
; REVISION HISTORY:
;   17-Jan-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_read_basel,lambda,flux,filename,teff=teff,logg=logg, mh=mh,vturb=vturb, $
                 xh=xh, silent=silent

; Need at least 3 parameters
if (N_params() LT 3) then begin
    print, 'Syntax - k_read_basel, lambda, flux, filename, [teff=,logg=, mh=,vturb=, xh=]'
    return
endif

lambda=fltarr(1221)

tmpflux=fltarr(1221)
tmpmodelno=0L
tmpteff=0L
tmplogg=0.d
tmpmh=0.d
tmpvturb=0.d
tmpxh=0.d

openr,unit,filename,/get_lun
readf,unit,lambda
nunits=0L
while(NOT eof(unit)) do begin
    readf,unit,tmpmodelno,tmpteff,tmplogg,tmpmh,tmpvturb,tmpxh
    readf,unit,tmpflux
    nunits=nunits+1L
end
if(NOT keyword_set(silent)) then begin
    outstr=strtrim(string(nunits),2)+' block(s) of spectra'
    klog,outstr
endif
close,unit
free_lun,unit

flux=fltarr(1221,nunits)
modelno=lonarr(nunits)
teff=lonarr(nunits)
logg=fltarr(nunits)
mh=fltarr(nunits)
vturb=fltarr(nunits)
xh=fltarr(nunits)

openr,unit,filename,/get_lun
readf,unit,lambda
nunits=0L
while(NOT eof(unit)) do begin
    readf,unit,tmpmodelno,tmpteff,tmplogg,tmpmh,tmpvturb,tmpxh
    readf,unit,tmpflux
    modelno[nunits]=tmpmodelno
    teff[nunits]=tmpteff
    logg[nunits]=tmplogg
    mh[nunits]=tmpmh
    vturb[nunits]=tmpvturb
    xh[nunits]=tmpxh
    flux[*,nunits]=tmpflux
    nunits=nunits+1L
end
close,unit
free_lun,unit

end
;------------------------------------------------------------------------------

