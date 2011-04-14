;+
; NAME:
;   wise_to_maggies
; PURPOSE:
;   convert WISE catalog input to Galactic-extinction corrected AB maggies 
; CALLING SEQUENCE:
;   wise_to_maggies,wise,mgy,ivar
; INPUTS:
;   wise - [N] structure with:
;               .RA
;               .DEC
;               .W1MPRO
;               .W1SIGMPRO
;               .W2MPRO
;               .W2SIGMPRO
;               .W3MPRO
;               .W3SIGMPRO
;               .W4MPRO
;               .W4SIGMPRO
; OUTPUTS:
;   mgy - [4, N] output in AB maggies in FUV and NUV filters
;   ivar - [4, N] inverse variance of maggies
; COMMENTS:
;   CCM Galactic extinction curve assumed (any reason that might be
;     right?)
;   No use of apertures yet.
; REVISION HISTORY:
;   14-Apr-2011  Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro wise_to_maggies, wise, mgy, ivar, nodust=nodust, mpro=mpro

filterlist= ['wise_w'+['1','2','3','4']+'.par']

; Handle Galactic reddening
if(NOT keyword_set(nodust)) then begin
    leff=k_lambda_eff(filterlist=filterlist)
    extvoebv=3.10
    extoextv=(ext_ccm(leff))
    
    glactc, wise.ra, wise.dec, 2000., gl, gb, 1, /deg
    ebv=dust_getval(gl,gb,/interp,/noloop)

    extinction=(extoextv#ebv)*extvoebv
endif else begin
    extinction= fltarr(4, n_elements(wise))
endelse

v2ab= k_vega2ab(filterlist=filterlist, /kurucz)

mgy= fltarr(4, n_elements(wise))
ivar= fltarr(4, n_elements(wise))

if(keyword_set(mpro)) then begin
    itag= lonarr(4)
    itag[0]= tag_indx(wise[0], 'w1mpro')
    itag[1]= tag_indx(wise[0], 'w2mpro')
    itag[2]= tag_indx(wise[0], 'w3mpro')
    itag[3]= tag_indx(wise[0], 'w4mpro')
    isigtag= lonarr(4)
    isigtag[0]= tag_indx(wise[0], 'w1sigmpro')
    isigtag[1]= tag_indx(wise[0], 'w2sigmpro')
    isigtag[2]= tag_indx(wise[0], 'w3sigmpro')
    isigtag[3]= tag_indx(wise[0], 'w4sigmpro')
endif else begin
    itag= lonarr(4)
    itag[0]= tag_indx(wise[0], 'w1mag_8')
    itag[1]= tag_indx(wise[0], 'w2mag_8')
    itag[2]= tag_indx(wise[0], 'w3mag_8')
    itag[3]= tag_indx(wise[0], 'w4mag_8')
    isigtag= lonarr(4)
    isigtag[0]= tag_indx(wise[0], 'w1sigm_8')
    isigtag[1]= tag_indx(wise[0], 'w2sigm_8')
    isigtag[2]= tag_indx(wise[0], 'w3sigm_8')
    isigtag[3]= tag_indx(wise[0], 'w4sigm_8')
endelse


minerr=0.01
for ib= 0L, 3L do begin
    mgy[ib,*]= 10.^(-0.4*(wise.(itag[ib])-extinction[ib,*]))
    igood= where(wise.(isigtag[ib]) ne 0., ngood)
    if(ngood gt 0) then $
      ivar[ib,igood]= 1./((0.4*alog(10.)*mgy[ib,igood]^2)* $
                          (wise[igood].(isigtag[ib])^2+minerr^2))
endfor

end

